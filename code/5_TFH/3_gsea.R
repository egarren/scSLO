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
# library(VISION)
library(msigdbr)
# library(msigdb)
# library(goseq)
# library(nicethings)
library(DOSE)
library(ggpubr)
library(AnnotationHub)
library(SPIA)

#load DEs
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
rm(list=ls(pattern=".txt.DE"))

i="tfh_pd1.DE"
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
      ggsave2(paste0(i,".gGO.bar.png"),width=4, height=2,device="png")
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


#custom GSEA
Idents(obj)<-"seurat_clusters"
DimPlot(obj,label=T)
i="tfh_pd1.DE" 
df<-get(i)#load DE of interest
df<-df[df$p_val_adj<0.01,] #select for sig genes
df$gene<-rownames(df)
ranks <- df$avg_log2FC
names(ranks) <- df$gene
sostdc1_sig<-unique(obj.markers$gene[obj.markers$cluster %in% c("3","7") &
                                   obj.markers$avg_log2FC>1 & obj.markers$p_val_adj <0.01])
tcm_sig<-unique(obj.markers$gene[obj.markers$cluster %in% c("5") &
                                  obj.markers$avg_log2FC>1 & obj.markers$p_val_adj <0.01])
sig_list<-list(sostdc1_sig,tcm_sig)
names(sig_list)<-c("sostdc1","tcm")
res<-fgsea(pathways=sig_list, ranks,eps=0)
plotEnrichment(sostdc1_sig, ranks,ticksSize = 0.05)+
  theme(plot.title = element_text(size=5),axis.title=element_blank(),axis.text=element_text(size=5))#plot top pathway enrichment
ggsave2(paste0(i,".sostdc1.gsea.png"),width=1.5, height=1.5,device="png")
plotEnrichment(tcm_sig, ranks,ticksSize = 0.05)+
  theme(plot.title = element_text(size=5),axis.title=element_blank(),axis.text=element_text(size=5))#plot top pathway enrichment
ggsave2(paste0(i,".tcm.gsea.png"),width=1.5, height=1.5,device="png")



#custom plot
tab<-try(enrichDAVID(sigGenes,universe=all_genes,annotation="GOTERM_BP_FAT",david.user="***"))
df<-as.data.frame(tab)
df$lp<- -log10(df$p.adjust)
df<-df[order(df$lp,decreasing=T),]
df<-head(df,n=20)
ggplot(df,aes(x=reorder(Description,lp),y=lp))+geom_col()+coord_flip()+theme_classic()+
  theme(text=element_text(family="Arial"))+
  geom_hline(yintercept=1,linetype=2,size=0.2)+labs(x=NULL,y=NULL)
ggsave2("tfh_gsea_table2.png",width=6, height=4,device="png")


#FCvsFC
m564_pd1.DE$gene<-rownames(m564_pd1.DE)
B6_pd1.DE$gene<-rownames(B6_pd1.DE)
df<-left_join(m564_pd1.DE[,c("avg_log2FC","gene")],B6_pd1.DE[,c("avg_log2FC","gene")],by="gene",keep=F)
colnames(df)<-c("m564FC","gene","b6FC")
rownames(df)<-df$gene
df$diff= abs(df$m564FC-df$b6FC)
df$sum = df$m564FC + df$b6FC
df<-df %>% mutate(rank_diff = rank(-diff), rank_sum = rank(sum), rank_sum_dec = rank(-sum)) 
gene.list<-df$gene[df$rank_diff<10 | df$rank_sum <5 |df$rank_sum_dec<10]
gene.list<-c(gene.list,"Sostdc1")
df$label<-"NS"
df$label[df$m564FC>1 & df$b6FC>1]<-"up"
df$label[df$m564FC< -1 & df$b6FC< -1]<-"down"
df$label[df$diff>1]<-"diff"
df<-df[order(-df$rank_diff),]
df$m564FC <- pmin(df$m564FC, 3)
df$b6FC <- pmin(df$b6FC, 3)


#custom volcano GSEA
genes.meta<-getBM(attributes=c("ensembl_gene_id", "mgi_symbol", "start_position", "end_position", "chromosome_name", 
                               "percentage_gene_gc_content", "external_gene_name", "gene_biotype","go_id","name_1006"),filters=
                    "mgi_symbol",values=list(rownames(obj@assays[["RNA"]]@meta.features)),
                  mart=useMart("ensembl", dataset = "mmusculus_gene_ensembl",host="ensembl.org"),useCache=F) #useast.


for(i in c("M2810","M15594","M13559","M11785","M15156","M2071","M12750","M18127","M14645","M11809","M12269","M11809","M18190","M12919",
           "GO:0042110","GO:0045333","GO:0008380","GO:0009060","GO:0009127","GO:0009168",
           "GO:0006754","GO:0030217","GO:0007159","GO:0035456","GO:0050900",
           "GO:0035458","GO:0030098","GO:0006119","GO:0009152")){
  if(grepl("GO:",i)){gene.list2<-genes.meta$mgi_symbol[genes.meta$go_id==i]}else{
    gene.list2<-unique(m_df$gene_symbol[m_df$gs_id==i]) #M18810
  }
  if(length(gene.list2)>5){
    j="tfh_pd1.DE"
    res<-get(j)
    thresh_p_val_adj <- 1e-10
    thresh_lfc <-0.3
    plt_df<- res %>% rownames_to_column(var = "gene") %>% 
      mutate(up_in = ifelse(p_val_adj >= thresh_p_val_adj, "NS", ifelse(avg_log2FC > thresh_lfc, "FDC", ifelse(avg_log2FC < -thresh_lfc, "TRC", "NS"))),
             rank_pval = rank(p_val_adj), rank_lfc_inc = rank(avg_log2FC), rank_lfc_dec = rank(-abs(avg_log2FC)),
             lp = -log10(p_val_adj)) %>% 
      arrange(-abs(avg_log2FC))
    plt_df <- plt_df %>% mutate(gene_label = ifelse(up_in != "NS" & gene %in% gene.list2 &
                                                      (rank_pval < 100 | rank_lfc_inc < 100 | rank_lfc_dec < 100), gene, NA))
    plt_df <- plt_df %>% mutate(label = ifelse(gene %in% gene.list2, "sig", "notsig"))
    plt_df<-plt_df[order(plt_df$label),]
    plt_df$lp <- pmin(plt_df$lp, 150)
    plt_df$avg_log2FC <- pmin(plt_df$avg_log2FC, 3)
    ggplot(plt_df, aes(x = avg_log2FC, y = lp, color = label, label = gene_label,size=label)) +
      geom_point() +
      geom_text_repel(color="black",size = 2, segment.size = 0.2,seed = 1, show.legend=T, box.padding=0.3,max.overlaps=20) +
      geom_hline(yintercept=10,linetype=2,size=0.2)+
      geom_vline(xintercept=0.3,linetype=2,size=0.2)+
      geom_vline(xintercept= -0.3,linetype=2,size=0.2)+
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
    ggsave2(paste0(j,"_volcano_gsea",i,".png"),width=2, height=1.7,device="png")
    
    gene.list<-df$gene[df$gene %in% gene.list2 & (df$rank_diff<100 | df$rank_sum <100 |df$rank_sum_dec<100)]
    df <- df %>% mutate(label = ifelse(gene %in% gene.list2, "sig", "notsig")) #genes to color
    df<-df[order(df$label),]
    df<- df %>% mutate(gene_label = ifelse(gene %in% gene.list , gene, NA))
    ggplot(df, aes(b6FC,m564FC,colour=label,size=label,label=gene_label)) + 
      geom_point()+#geom_point(fill=NA,colour=alpha("black",0.5),pch=21,size=3) +
      theme_classic()+
      scale_size_manual(values=c(0.5,0.7))+
      geom_text_repel(color="black",size = 2, segment.size = 0.2,seed = 1, show.legend=T, box.padding=0.3,max.overlaps=20) +
      theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position="none",
            axis.ticks=element_blank(),axis.line=element_blank())+ #axis.text=element_blank(),
      geom_hline(yintercept=0,color=alpha("black",0.5))+geom_vline(xintercept=0,color=alpha("black",0.5))+#geom_abline(intercept = 0, slope = 1)+
      scale_color_manual(values = c("sig" ="black", "notsig" = "grey70"))
    ggsave2(paste0(i,"_FCvsFC.png"),width=3, height=3,device="png")
    
    obj3 <- AddModuleScore(object = obj2,features = list(gene.list2),name = 'temp')
    FeaturePlot(obj3, "temp1") +
      NoLegend()+ NoAxes()+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
                                 plot.title=element_blank()) & #NoLegend()+
      scale_color_gradientn(colors = viridis(n = 10, direction = -1))
    ggsave2(paste0(i,"_umap.png"),width=3.5, height=3.5,device="png")
    
  }
}
save.image("gsea.RData")


FeaturePlot(obj3, "temp1") & scale_color_gradientn(colors = viridis(n = 10, direction = -1))
ggsave2("enrichment_legend.png",width=3.5, height=3.5,device="png")



