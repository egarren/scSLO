rm(list=ls())
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
library(ggrepel)

#load DEs
dir.create("./gsea")
setwd("./gsea")
load("../temp2.DEseq.RData")
m_df = msigdbr(species = "Mus musculus")#, category = "C7") #H = hallmarks, C2=curated, C5=GO, C7=immune
pathways<-m_df %>% split(x = .$gene_symbol, f = .$gs_name) 
msig.df<-m_df %>% dplyr::select(gs_name,entrez_gene)
ah = AnnotationHub()
ms_ens <- query(ah, c("Mus musculus", "EnsDb"))
ms_ens <- ms_ens[["AH89211"]]
annotations_ahb <- ensembldb::genes(ms_ens, return.type = "data.frame")

#GSEA
for(i in ls(pattern="res_")){
  DE.df<-get(i)
  DE.df<-data.frame(ensembl=rownames(DE.df),log2FC=DE.df$log2FoldChange,p_val_adj=DE.df$padj)
  DE.df<-left_join(DE.df,unique(tx2gene[,c("gene_id","symbol")]),by=c("ensembl"="gene_id"),keep=F)
  DE.df<-DE.df[!is.na(DE.df$p_val_adj),c("log2FC","p_val_adj","symbol")]
  DE.df$gene<-DE.df$symbol
  DE.df$avg_log2FC = -DE.df$log2FC
  if(sum(DE.df$avg_log2FC)!=0){
    ## FSGEA
    df<-DE.df
    df<-df[df$p_val_adj<0.05,] 
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
        topUp <- head(fgseaRes %>% dplyr::filter(ES > 0) %>% dplyr::top_n(5, wt=-padj),n=5)
        topDown <- head(fgseaRes %>% dplyr::filter(ES < 0) %>% top_n(5, wt=-padj),n=5)
        topPathways <- bind_rows(topUp, topDown) %>% arrange(-ES)
        png(paste0(i,".",h,".gsea.table.png"),width=10,height=4,units="in",res=200)
        plotGseaTable(path[topPathways$pathway], ranks, fgseaRes, gseaParam = 0.5)
        dev.off()
        write.csv(fgseaRes[order(padj, -abs(NES)), ][,1:7],file=paste0(i,".gsea.results.csv"))
        write.csv(topPathways[,1:7],file=paste0(i,".gsea.table.csv"))
      }}
    }
    
    df<-DE.df
    df<-left_join(df,annotations_ahb, by=c("gene"="symbol"))    
    sigGenes <- as.character(na.exclude(df$entrezid[df$p_val_adj < 0.05 & df$avg_log2FC >1])) 
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
      ggsave2(paste0(i,".gGO.bar.png"),width=4, height=2,device="png")
    }}
    for(j in c("eGO","gsea","kk","kk2","mkk","mkk2","davidKEGG","davidBP","egmt","egmt2","ctd","ctd2")){
      if(j=="eGO"){tab<-try(enrichGO(gene=sigGenes,universe=all_genes,OrgDb=org.Mm.eg.db,ont="ALL", readable=T))}
      if(j=="gsea"){tab<-try(gseGO(geneList,OrgDb=org.Mm.eg.db,ont="ALL",pvalueCutoff = 0.05))}
      if(j=="kk"){tab <- try(enrichKEGG(gene = sigGenes,universe=all_genes,organism = 'mmu'))}#search_kegg_organism('mmu', by='kegg_code')}
      if(j=="kk2"){tab<-try(gseKEGG(geneList,organism="mmu",pvalueCutoff=0.05))} #KEGG gsea}
      if(j=="mkk"){tab<-try(enrichMKEGG(sigGenes,universe=all_genes,organism="mmu"))}
      if(j=="mkk2"){tab<-try(gseMKEGG(geneList,organism="mmu",pvalueCutoff=0.05)) }
      if(j=="davidKEGG"){tab<-try(enrichDAVID(sigGenes,universe=all_genes,annotation="KEGG_PATHWAY",david.user="XXX")) }
      if(j=="davidBP"){tab<-try(enrichDAVID(sigGenes,universe=all_genes,annotation="GOTERM_BP_FAT",david.user="XXX"))}
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
          for(m in 1:pmin(length(tab$ID[m]),5)){
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

#custom GSEA
load("analyzed.RData")
Idents(obj)<-"seurat_clusters"
DimPlot(obj,label=T)
i="res_genotype_TFH" 
df<-get(i)
df<-data.frame(ensembl=rownames(df),log2FC=df$log2FoldChange,p_val_adj=df$padj)
df<-left_join(df,unique(tx2gene[,c("gene_id","symbol")]),by=c("ensembl"="gene_id"),keep=F)
df<-df[!is.na(df$p_val_adj),c("log2FC","p_val_adj","symbol")]
df$gene<-df$symbol
df$avg_log2FC = -df$log2FC
df<-df[df$p_val_adj<0.05,] 
rownames(df)<-df$gene
ranks <- df$avg_log2FC
names(ranks) <- df$gene
sostdc1_sig<-unique(obj.markers$gene[obj.markers$cluster %in% c("3","7") &
                                   obj.markers$avg_log2FC>0.7 & obj.markers$p_val_adj <0.05])
tcm_sig<-unique(obj.markers$gene[obj.markers$cluster %in% c("5") &
                                  obj.markers$avg_log2FC>0.7 & obj.markers$p_val_adj <0.05])
sig_list<-list(sostdc1_sig,tcm_sig)
names(sig_list)<-c("sostdc1","tcm")
res<-fgsea(pathways=sig_list, ranks,eps=0)
plotEnrichment(sostdc1_sig, ranks,ticksSize = 0.05)+
  theme(plot.title = element_text(size=5),axis.title=element_blank(),axis.text=element_text(size=5))#plot top pathway enrichment
ggsave2(paste0("../output/",i,".sostdc1.gsea.png"),width=1.5, height=1.5,device="png")
plotEnrichment(tcm_sig, ranks,ticksSize = 0.05)+
  theme(plot.title = element_text(size=5),axis.title=element_blank(),axis.text=element_text(size=5))#plot top pathway enrichment
ggsave2(paste0("../output/",i,".tcm.gsea.png"),width=1.5, height=1.5,device="png")

#custom plot
i="res_genotype_TFH" 
df<-get(i)
df<-data.frame(ensembl=rownames(df),log2FC=df$log2FoldChange,p_val_adj=df$padj)
df<-left_join(df,unique(tx2gene[,c("gene_id","symbol")]),by=c("ensembl"="gene_id"),keep=F)
df<-df[!is.na(df$p_val_adj),c("log2FC","p_val_adj","symbol")]
df$gene<-df$symbol
df$avg_log2FC = -df$log2FC
df<-df %>%
  dplyr::filter(!(grepl("Rps", gene) | grepl("Rpl", gene)| grepl("mt.", gene)| grepl("H2.", gene))) 
df<-left_join(df,annotations_ahb, by=c("gene"="symbol"))    
sigGenes <- as.character(na.exclude(df$entrezid[df$p_val_adj < 0.05 & df$avg_log2FC >1]))
all_genes<-as.character(df$entrezid)
tab<-try(enrichDAVID(sigGenes,universe=all_genes,annotation="GOTERM_BP_FAT",david.user="XXX"))
tab<-as.data.frame(tab)
tab$lp<- -log10(tab$p.adjust)
tab<-tab[order(tab$lp,decreasing=T),]
tab<-head(tab,n=20)
ggplot(tab,aes(x=reorder(Description,lp),y=lp))+geom_col()+coord_flip()+theme_classic()+
  theme(text=element_text(family="Arial"))+
  geom_hline(yintercept=1,linetype=2,size=0.2)+labs(x=NULL,y=NULL)
ggsave2("../output/tfh_gsea_table2.png",width=6, height=4,device="png")
write.csv(tab,"../output/tfh_gsea_tab.csv")

#Volcano
gene.list2<-sostdc1_sig
i="res_genotype_TFH"
res_df<-get(i)
res_df<- res_df[order(res_df$pvalue),]
volcano<-data.frame(gene=rownames(res_df),log2FC=res_df$log2FoldChange,p_val_adj=res_df$padj)
volcano<-left_join(volcano,unique(tx2gene[,c("gene_id","symbol")]),by=c("gene"="gene_id"),keep=F)
volcano<-volcano[!is.na(volcano$p_val_adj),c("log2FC","p_val_adj","symbol")]
volcano$gene<-volcano$symbol
volcano$avg_log2FC = -volcano$log2FC
thresh_p_val_adj <- 0.05
thresh_lfc <-1
plt_df<- volcano %>% #rownames_to_column(var = "gene") %>% 
  dplyr::filter(!(grepl("Rps", gene) | grepl("Rpl", gene)| grepl("mt.", gene)|grepl("Igh", gene)| grepl("H2.", gene))) %>%
  mutate(up_in = ifelse(p_val_adj >= thresh_p_val_adj, "NS", ifelse(avg_log2FC > thresh_lfc, "UP", ifelse(avg_log2FC < -thresh_lfc, "DOWN", "NS"))),
         rank_pval = rank(p_val_adj), rank_lfc_inc = rank(avg_log2FC), rank_lfc_dec = rank(-abs(avg_log2FC)),
         lp = -log10(p_val_adj)) %>% 
  arrange(-abs(avg_log2FC))
plt_df <- plt_df %>% mutate(gene_label = ifelse(up_in != "NS" & gene %in% gene.list2 &
                                                  (rank_pval < 50 | rank_lfc_inc < 50 | rank_lfc_dec < 50), gene, NA))
plt_df <- plt_df %>% mutate(label = ifelse(gene %in% gene.list2, "sig", "notsig"))
plt_df<-plt_df[order(plt_df$label),]
plt_df$lp <- pmin(plt_df$lp, 10)
plt_df$avg_log2FC <- pmin(plt_df$avg_log2FC, 5)
plt_df$avg_log2FC <- pmax(plt_df$avg_log2FC, -5)
ggplot(plt_df, aes(x = avg_log2FC, y = lp, color = label, label = gene_label)) +
  geom_point(size=0.5) +
  geom_text_repel(color="black",size = 2, segment.size = 0.2,seed = 1, show.legend=T, box.padding=0.3,max.overlaps=20) +
  geom_hline(yintercept=1.3,linetype=2,size=0.2)+
  geom_vline(xintercept=1,linetype=2,size=0.2)+
  geom_vline(xintercept= -1,linetype=2,size=0.2)+
  scale_size_manual(values=c(0.5,0.7))+
  scale_color_manual(values = c("sig" = "black", "notsig" = "grey80")) +
  scale_x_continuous(limits = c(-max(plt_df$avg_log2FC), max(plt_df$avg_log2FC)), expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(limits = c(0, max(plt_df$lp)), expand = expansion(mult = c(0.05, 0.01))) +
  labs(x=NULL,y=NULL) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text=element_text(size=5),
        panel.background = element_blank(),legend.position="none",#legend.title=element_text(size=8),
        plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 10, color = "black"))
ggsave2(paste0("../output/volcano_gsea_",i,".png"),width=3, height=2.5,device="png")

