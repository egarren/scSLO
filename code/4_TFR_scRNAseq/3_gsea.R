rm(list=ls())
###Pathway analysis
library(org.Mm.eg.db)
library(tidyverse)
library(RDAVIDWebService)
library(Seurat)
library(cowplot)
library(fgsea)
library(clusterProfiler)
library(enrichplot)
library(pathview)
library(topGO)
library(scde)
library(biomaRt)
library(GO.db)
library(DBI)
library(msigdbr)
library(DOSE)
library(ggpubr)
library(AnnotationHub)
library(SPIA)

#load DEs
dir.create("./gsea")
setwd("./gsea")
load("../analyzed.RData")
m_df = msigdbr(species = "Mus musculus")#, category = "C7") #H = hallmarks, C2=curated, C5=GO, C7=immune
pathways<-m_df %>% split(x = .$gene_symbol, f = .$gs_name) 
msig.df<-m_df %>% dplyr::select(gs_name,entrez_gene)
ah = AnnotationHub()
ms_ens <- query(ah, c("Mus musculus", "EnsDb"))
ms_ens <- ms_ens[["AH89211"]]
annotations_ahb <- ensembldb::genes(ms_ens, return.type = "data.frame")

#custom genesets
load("custom.gse.RData")
pathways<-append(pathways,gse.list) #,CTD
gse.entrez2<-list(gse.entrez) %>%
  map_df(enframe, name = "gs_name", value="entrez_gene") %>% 
  unnest
gse.entrez2$entrez_gene<-as.numeric(as.character(gsub("///.*","",gse.entrez2$entrez_gene)))
gse.entrez2<-gse.entrez2[!is.na(gse.entrez2$entrez_gene),]
msig.df<-bind_rows(msig.df,gse.entrez2) #,CTD.entrez
CTD.entrez<-list(CTD.entrez) %>%
  map_df(enframe, name = "gs_name", value="entrez_gene") %>% 
  unnest
CTD.entrez$entrez_gene<-as.numeric(as.character(gsub("///.*","",CTD.entrez$entrez_gene)))
CTD.entrez<-CTD.entrez[!is.na(CTD.entrez$entrez_gene),]



for(i in ls(pattern=".DE")){
  DE.df<-get(i)
  DE.df<-DE.df %>% rownames_to_column(var = "gene") %>% 
    dplyr::filter(!(grepl("Rps", gene) | grepl("Rpl", gene)| grepl("mt.", gene)| grepl("H2.", gene))) 
  if(sum(DE.df$avg_log2FC)!=0){
    ## FSGEA
    df<-DE.df
    df<-df[df$p_val_adj<0.05,] #select for sig genes, abs(df$avg_log2FC)>0.2&
    ranks <- df$avg_log2FC
    names(ranks) <- df$gene
    png(paste0(i,".ranks.png"),width=4,height=4,units="in",res=200)
    barplot(sort(ranks, decreasing = T))
    dev.off()
    for(h in c("pathways","CTD")){
      path<-get(h)
      fgseaRes <- fgsea(path, ranks, minSize=15, maxSize = 500, nperm=1000)
      if(!is.null(fgseaRes)){if(nrow(fgseaRes)>1){
        #enrichment plot
        head(fgseaRes[order(padj, -abs(NES)), ], n=15)
        plotEnrichment(path[[fgseaRes[order(padj, -abs(NES)), ]$pathway[1]]], ranks)+
          ggtitle(fgseaRes[order(padj, -abs(NES)), ]$pathway[1])+
          theme(plot.title = element_text(size=4),axis.title=element_text(size=5),axis.text=element_text(size=5))#plot top pathway enrichment
        ggsave2(paste0(i,".",h,".topenrichment.png"),width=2, height=2,device="png")
        write.csv(fgseaRes[order(padj, -abs(NES)), ][,1:7],file=paste0(i,".gsea.results.csv"))
      }}
    }
    
    ##clusterProfiler (https://bioconductor.statistik.tu-dortmund.de/packages/3.6/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html)
    df<-DE.df
    df<-left_join(df,annotations_ahb, by=c("gene"="symbol"))    
    sigGenes <- as.character(na.exclude(df$entrezid[df$p_val_adj < 0.05])) #select sig genes, abs(df$avg_log2FC) 
    all_genes<-as.character(df$entrezid)
    df2<-df[!is.na(df$entrezid),]
    df2<-df2[which(duplicated(df2$entrezid) == F),]
    ranks <- df2$avg_log2FC
    names(ranks) <- df2$entrezid
    geneList <- sort(ranks, decreasing = TRUE)
    #topGO
    ggo.table<-groupGO(gene=sigGenes,OrgDb=org.Mm.eg.db,ont="BP",level=4,readable=T) #BP (biological process), MF (molecular function), CC (cellular componnent)
    if(!is.null(ggo.table)){if(dim(ggo.table)[1]!=0){
      barplot(ggo.table, drop=TRUE, showCategory=6)+ggtitle(i)
      ggsave2(paste0(i,".gGO.bar.png"),width=8, height=2,device="png")
    }}
    for(j in c("eGO","gsea","kk","kk2","mkk","mkk2","davidKEGG","davidBP","egmt","egmt2","ctd","ctd2")){
      if(j=="eGO"){tab<-try(enrichGO(gene=sigGenes,universe=all_genes,OrgDb=org.Mm.eg.db,ont="ALL", readable=T))}
      if(j=="gsea"){tab<-try(gseGO(geneList,OrgDb=org.Mm.eg.db,ont="ALL",pvalueCutoff = 0.05))}
      if(j=="kk"){tab <- try(enrichKEGG(gene = sigGenes,universe=all_genes,organism = 'mmu'))}#search_kegg_organism('mmu', by='kegg_code')}
      if(j=="kk2"){tab<-try(gseKEGG(geneList,organism="mmu",pvalueCutoff=0.05))} #KEGG gsea}
      if(j=="mkk"){tab<-try(enrichMKEGG(sigGenes,universe=all_genes,organism="mmu"))}
      if(j=="mkk2"){tab<-try(gseMKEGG(geneList,organism="mmu",pvalueCutoff=0.05)) }
      if(j=="davidKEGG"){tab<-try(enrichDAVID(sigGenes,universe=all_genes,annotation="KEGG_PATHWAY",david.user="***")) }
      if(j=="davidBP"){tab<-try(enrichDAVID(sigGenes,universe=all_genes,annotation="GOTERM_BP_FAT",david.user="***"))}
      if(j=="egmt"){tab<-try(enricher(sigGenes,universe=all_genes,TERM2GENE = msig.df))}
      if(j=="egmt2"){tab<-try(GSEA(geneList,TERM2GENE = msig.df,pvalueCutoff=0.05))}
      if(j=="ctd"){tab<-try(enricher(sigGenes,universe=all_genes,TERM2GENE = CTD.entrez))}
      if(j=="ctd2"){tab<-try(GSEA(geneList,TERM2GENE = CTD.entrez,pvalueCutoff=0.05))}
      if(!is.null(tab)&!is(tab,"try-error")){if(nrow(tab)>=1){ 
        if(j %in% c("eGO","gsea")){tab<-setReadable(tab,OrgDb=org.Mm.eg.db)}else{tab<-setReadable(tab,OrgDb=org.Mm.eg.db,keyType="ENTREZID")}
        write.csv(tab,file=paste0(i,".",j,".tab.csv"))
        if(j %in% c("eGO","kk","mkk","davidKEGG")){
          barplot(tab, drop=TRUE, showCategory=6)+ggtitle(i)
          ggsave2(paste0(i,".",j,".bar.png"),width=8, height=2,device="png")
        }
        clusterProfiler::dotplot(tab, showCategory=6)+ggtitle(i)+  theme(legend.direction = "vertical", legend.box = "horizontal")
        ggsave2(paste0(i,".",j,".dot.png"),width=9, height=5,device="png")
        clusterProfiler::dotplot(tab, showCategory=6,font.size=7)+ggtitle(i)+  
          guides(size=F)+scale_color_continuous(name="P-val",guide=guide_colorbar(reverse=TRUE))
        ggsave2(paste0(i,".",j,".dot2.png"),width=4.25, height=2,device="png")
        clusterProfiler::cnetplot(tab, categorySize="pvalue", foldChange=geneList)+
          ggtitle(i)+scale_color_gradient2(name="LogFC")
        ggsave2(paste0(i,".",j,".cnet.png"),width=6, height=4,device="png")
        clusterProfiler::cnetplot(tab, categorySize="pvalue", foldChange=geneList,node_label="gene")+
          guides(size=F)+scale_color_gradient2(name="LogFC")
        ggsave2(paste0(i,".",j,".cnet2.png"),width=4.5, height=3.5,device="png")
           if(j %in% c("gsea","kk2","egmt2","ctd2")){
          clusterProfiler::gseaplot(tab,geneSetID=tab$ID[1],title=tab$Description[1])
          ggsave2(paste0(i,".",j,".plot.png"),width=4, height=5,device="png")
          for(m in 1:10){
            clusterProfiler::gseaplot(tab,geneSetID=tab$ID[m],title=tab$Description[m])
            ggsave2(paste0(i,".",j,".",m,".plot.png"),width=3, height=5,device="png")
          }
        }
        if(j %in% c("kk","kk2","davidKEGG")){
          pathview(gene.data = geneList, pathway.id = tab$ID[2], species = "mmu", out.suffix=paste0(i,".",j,".pathview"))
        }
      }}
    }
  }
}



#custom plot
tab<-try(enrichDAVID(sigGenes,universe=all_genes,annotation="GOTERM_BP_FAT",david.user="***"))
df<-as.data.frame(tab)
df$lp<- -log10(df$p.adjust)
df<-df[order(df$lp,decreasing=T),]
df<-head(df,n=15)
ggplot(df,aes(x=reorder(Description,lp),y=lp))+geom_col()+coord_flip()+theme_classic()+
  theme(text=element_text(family="Arial"))+
  geom_hline(yintercept=1,linetype=2,size=0.2)+labs(x=NULL,y=NULL)
ggsave2("tfr_gsea_table.png",width=4, height=3,device="png")

