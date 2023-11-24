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


m_df = msigdbr(species = "Homo sapiens")#, category = "C7") #H = hallmarks, C2=curated, C5=GO, C7=immune
pathways<-m_df %>% split(x = .$gene_symbol, f = .$gs_name) 
for(i in c(ls(pattern="DE"))){
  df<-get(i)#load DE of interest
  df<-df[df$p_val_adj<0.01,] #select for sig genes
  df$gene<-rownames(df)
  ranks <- df$avg_log2FC
  names(ranks) <- df$gene
  png(paste0(i,".ranks.png"),width=4,height=4,units="in",res=200)
  barplot(sort(ranks, decreasing = T))
  dev.off()
  fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize = 500, nperm=1000)
  if(!is.null(fgseaRes)){if(dim(fgseaRes)[1]!=0){
    #enrichment plot
    head(fgseaRes[order(padj, -abs(NES)), ], n=15)
    plotEnrichment(pathways[[fgseaRes[order(padj, -abs(NES)), ]$pathway[1]]], ranks)+
      ggtitle(fgseaRes[order(padj, -abs(NES)), ]$pathway[1])+
      theme(plot.title = element_text(size=5),axis.title=element_text(size=5),axis.text=element_text(size=5))#plot top pathway enrichment
    ggsave2(paste0(i,".topenrichment.png"),width=3, height=3,device="png")
    #gsea table
    topUp <- head(fgseaRes %>% dplyr::filter(ES > 0) %>% top_n(10, wt=-padj),n=10)
    topDown <- head(fgseaRes %>% dplyr::filter(ES < 0) %>% top_n(10, wt=-padj),n=10)
    topPathways <- bind_rows(topUp, topDown) %>% arrange(-ES)
    png(paste0(i,".gsea.table.png"),width=20,height=7,units="in",res=200)
    plotGseaTable(pathways[topPathways$pathway], ranks, fgseaRes, gseaParam = 0.5)
    dev.off()
    write.csv(fgseaRes[order(padj, -abs(NES)), ][,1:7],file=paste0(i,".gsea.results.csv"))
    write.csv(topPathways[,1:7],file=paste0(i,".gsea.table.csv"))
  }}
}

#custom GSEA
i="pdl1.DE" 
df<-get(i)#load DE of interest
df$gene<-rownames(df)
ranks <- df$avg_log2FC
names(ranks) <- df$gene
asc_sig<-unique(obj.markers$gene[obj.markers$cluster %in% c("ASC_1","PC","ASC_2","ASC_3") &
                                   obj.markers$avg_log2FC>1 & obj.markers$p_val_adj <0.01])
gc_sig<-unique(obj.markers$gene[obj.markers$cluster %in% c("GC_1","GC_2","LZ","DZ") &
                                   obj.markers$avg_log2FC>1 & obj.markers$p_val_adj <0.01])
sig_list<-list(asc_sig,gc_sig,asc2_sig)
names(sig_list)<-c("ASC","GC","ASC_2")
res<-fgsea(pathways=sig_list, ranks,eps=0,nPermSimple=10000)
plotEnrichment(asc_sig, ranks,ticksSize = 0.05)+
  theme(plot.title = element_text(size=5),axis.title=element_blank(),axis.text=element_text(size=5))#plot top pathway enrichment
ggsave2(paste0(i,".asc.gsea.png"),width=1.5, height=1.5,device="png")
plotEnrichment(gc_sig, ranks,ticksSize = 0.05)+
  theme(plot.title = element_text(size=5),axis.title=element_blank(),axis.text=element_text(size=5))#plot top pathway enrichment
ggsave2(paste0(i,".gc.gsea.png"),width=1.5, height=1.5,device="png")



##clusterProfiler (https://bioconductor.statistik.tu-dortmund.de/packages/3.6/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html)
m_df = msigdbr(species = "Homo sapiens")#, category = "C7") #H = hallmarks, C2=curated, C5=GO, C7=immune
msig.df<-m_df %>% dplyr::select(gs_name,entrez_gene)
for(i in c("Naive_pdl1.DE","pd1.DE","pdl1.DE")){ #ls(pattern="DE")
  df<-get(i)
  df$gene<-rownames(df)
  genes.meta<-bitr(df$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") #pull in entrezID, could also use bitr_kegg, see: keytypes(org.Hs.eg.db)
  df<-cbind(df,genes.meta[match(df$gene, genes.meta$SYMBOL),])
  all_genes<-na.exclude(as.character(df$ENTREZID))
  sigGenes <- na.exclude(df$ENTREZID[df$p_val_adj < 0.01 &df$avg_log2FC > 1 ]) #select sig genes, abs(df$avg_log2FC) 
  df2<-df[df$p_val_adj<0.01 & !is.na(df$ENTREZID),]
  ranks <- df2$avg_log2FC
  names(ranks) <- df2$ENTREZID
  geneList <- sort(ranks, decreasing = TRUE)
  #topGO
  ggo.table<-groupGO(gene=sigGenes,OrgDb="org.Hs.eg.db",ont="MF",level=4,readable=T) #BP (biological process), MF (molecular function), CC (cellular componnent)
  if(!is.null(ggo.table)){if(dim(ggo.table)[1]!=0){
    write.csv(ggo.table,file=paste0(i,".gGO.table.csv"))
    barplot(ggo.table, drop=TRUE, showCategory=12)+ggtitle(i)
    ggsave2(paste0(i,".gGO.bar.png"),width=8, height=4,device="png")
  }}
  ego.table<-enrichGO(gene=sigGenes,OrgDb="org.Hs.eg.db",ont="MF", readable=T)
  if(!is.null(ego.table)){ if(dim(ego.table)[1]!=0){
    write.csv(ego.table,file=paste0(i,".eGO.table.csv"))
    barplot(ego.table, drop=TRUE, showCategory=12)+ggtitle(i)
    ggsave2(paste0(i,".eGO.bar.png"),width=7, height=4,device="png")
    clusterProfiler::dotplot(ego.table)+ggtitle(i)
    ggsave2(paste0(i,".eGO.dot.png"),width=7, height=4,device="png")
    x2<-enrichplot::pairwise_termsim(ego.table)
    emapplot(x2)+ggtitle(i)
    ggsave2(paste0(i,".eGO.emap.png"),width=10, height=8,device="png")
    clusterProfiler::cnetplot(ego.table, categorySize="pvalue", foldChange=geneList)+ggtitle(i)
    ggsave2(paste0(i,".eGO.cnet.png"),width=10, height=8,device="png")
    png(paste0(i,".eGO.graph.png"),width=20,height=15,units="in",res=200)
    clusterProfiler::plotGOgraph(ego.table)
    dev.off()
  }}
  #gsea
  gsea.table<-gseGO(geneList,OrgDb="org.Hs.eg.db",ont="MF",pvalueCutoff = 0.2)
  if(!is.null(gsea.table)){if(dim(gsea.table)[1]!=0){
    gsea.table<-setReadable(gsea.table,OrgDb="org.Hs.eg.db")
    write.csv(gsea.table,file=paste0(i,".gsea.table.csv"))
    clusterProfiler::dotplot(gsea.table)+ggtitle(i)
    ggsave2(paste0(i,".gsea.dot.png"),width=7, height=4,device="png")
    x2<-enrichplot::pairwise_termsim(gsea.table)
    emapplot(x2)+ggtitle(i)
    ggsave2(paste0(i,".gsea.emap.png"),width=10, height=8,device="png")
    clusterProfiler::cnetplot(gsea.table, categorySize="pvalue", foldChange=geneList)+ggtitle(i)
    ggsave2(paste0(i,".gsea.cnet.png"),width=10, height=8,device="png")
    png(paste0(i,".gsea.graph.png"),width=20,height=15,units="in",res=200)
    clusterProfiler::plotGOgraph(gsea.table)
    dev.off()
    clusterProfiler::gseaplot(gsea.table,geneSetID=gsea.table$ID[1],title=gsea.table$Description[1])
    ggsave2(paste0(i,".gsea.plot.png"),width=6, height=9,device="png")
  }}
  #KEGG
  kk.table <- enrichKEGG(gene = sigGenes,organism = 'hsa')#search_kegg_organism('mmu', by='kegg_code')
  if(!is.null(kk.table)){ if(dim(kk.table)[1]!=0){
    kk.table <- setReadable(kk.table,OrgDb="org.Hs.eg.db",keyType="ENTREZID")
    write.csv(kk.table,file=paste0(i,".kegg.table.csv"))
    barplot(kk.table, drop=TRUE, showCategory=12)+ggtitle(i)
    ggsave2(paste0(i,".kk.bar.png"),width=7, height=4,device="png")
    clusterProfiler::dotplot(kk.table)+ggtitle(i)
    ggsave2(paste0(i,".kk.dot.png"),width=7, height=4,device="png")
    x2<-enrichplot::pairwise_termsim(kk.table)
    emapplot(x2)+ggtitle(i)
    ggsave2(paste0(i,".kk.emap.png"),width=10, height=8,device="png")
    clusterProfiler::cnetplot(kk.table, categorySize="pvalue", foldChange=geneList)+ggtitle(i)
    ggsave2(paste0(i,".kk.cnet.png"),width=10, height=8,device="png")
    pathview(gene.data = geneList, pathway.id = kk.table$ID[1], species = "hsa", out.suffix=paste0(i,".KEGGpath"))
  }}
  kk2.table<-gseKEGG(geneList,organism="hsa",pvalueCutoff=0.2) #KEGG gsea
  if(!is.null(kk2.table)){if(dim(kk2.table)[1]!=0){
    kk2.table<-setReadable(kk2.table,OrgDb="org.Hs.eg.db",keyType="ENTREZID") #KEGG gsea
    write.csv(kk2.table,file=paste0(i,".kk2.table.csv"))
    clusterProfiler::dotplot(kk2.table)+ggtitle(i)
    ggsave2(paste0(i,".kk2.dot.png"),width=7, height=4,device="png")
    emapplot(enrichplot::pairwise_termsim(kk2.table))+ggtitle(i)
    ggsave2(paste0(i,".kk2.emap.png"),width=10, height=8,device="png")
    clusterProfiler::cnetplot(kk2.table, categorySize="pvalue", foldChange=geneList)+ggtitle(i)
    ggsave2(paste0(i,".kk2.cnet.png"),width=10, height=8,device="png")
    clusterProfiler::gseaplot(kk2.table,geneSetID=kk2.table$ID[1],title=kk2.table$Description[1])
    ggsave2(paste0(i,".kk2.plot.png"),width=6, height=9,device="png")
    pathview(gene.data = geneList, pathway.id = kk2.table$ID[1], species = "hsa", out.suffix=paste0(i,".KEGG.gsea.path"))
  }}
  mkk.table<-enrichMKEGG(sigGenes,organism="hsa") #KEGG Module
  if(!is.null(mkk.table)){if(dim(mkk.table)[1]!=0){ 
    mkk.table<-setReadable(mkk.table,OrgDb="org.Hs.eg.db",keyType="ENTREZID") #KEGG Module
    write.csv(mkk.table,file=paste0(i,".mkk.table.csv"))
    barplot(mkk.table, drop=TRUE, showCategory=12)+ggtitle(i)
    ggsave2(paste0(i,".mkk.bar.png"),width=7, height=4,device="png")
    clusterProfiler::dotplot(mkk.table)+ggtitle(i)
    ggsave2(paste0(i,".mkk.dot.png"),width=7, height=4,device="png")
    emapplot(enrichplot::pairwise_termsim(mkk.table))+ggtitle(i)
    ggsave2(paste0(i,".mkk.emap.png"),width=10, height=8,device="png")
    clusterProfiler::cnetplot(mkk.table, categorySize="pvalue", foldChange=geneList)+ggtitle(i)
    ggsave2(paste0(i,".mkk.cnet.png"),width=10, height=8,device="png")
  }}
  mkk2.table<-gseMKEGG(geneList,organism="hsa",pvalueCutoff=0.2) #KEGG Module gsea
  if(!is.null(mkk2.table)){if(dim(mkk2.table)[1]!=0){
    mkk2.table<-setReadable(mkk2.table,OrgDb="org.Hs.eg.db",keyType="ENTREZID") #KEGG Module gsea
    write.csv(mkk2.table,file=paste0(i,".mkk2.table.csv"))
    clusterProfiler::dotplot(mkk2.table)+ggtitle(i)
    ggsave2(paste0(i,".mkk2.dot.png"),width=7, height=4,device="png")
    emapplot(enrichplot::pairwise_termsim(mkk2.table))+ggtitle(i)
    ggsave2(paste0(i,".mkk2.emap.png"),width=10, height=8,device="png")
    clusterProfiler::cnetplot(mkk2.table, categorySize="pvalue", foldChange=geneList)+ggtitle(i)
    ggsave2(paste0(i,".mkk2.cnet.png"),width=10, height=8,device="png")
    clusterProfiler::gseaplot(mkk2.table,geneSetID=mkk2.table$ID[1],title=mkk2.table$Description[1])
    ggsave2(paste0(i,".mkk2.plot.png"),width=6, height=9,device="png")
  }}
  david.KEGG.table<-enrichDAVID(sigGenes,annotation="KEGG_PATHWAY",david.user="elliot_akama-garren@hms.harvard.edu") 
  if(!is.null(david.KEGG.table)){if(dim(david.KEGG.table)[1]!=0){
    david.KEGG.table<-setReadable(david.KEGG.table,OrgDb="org.Hs.eg.db",keyType="ENTREZID") 
    write.csv(david.KEGG.table,file=paste0(i,".david.kegg.table.csv"))
    barplot(david.KEGG.table, drop=TRUE, showCategory=12)+ggtitle(i)
    ggsave2(paste0(i,".david.KEGG.bar.png"),width=7, height=4,device="png")
    clusterProfiler::dotplot(david.KEGG.table)+ggtitle(i)
    ggsave2(paste0(i,".david.KEGG.dot.png"),width=7, height=4,device="png")
    emapplot(enrichplot::pairwise_termsim(david.KEGG.table))+ggtitle(i)
    ggsave2(paste0(i,".david.KEGG.emap.png"),width=10, height=8,device="png")
    clusterProfiler::cnetplot(david.KEGG.table, categorySize="pvalue", foldChange=geneList)+ggtitle(i)
    ggsave2(paste0(i,".david.KEGG.cnet.png"),width=10, height=8,device="png")
    pathview(gene.data = geneList, pathway.id = david.KEGG.table$ID[1], species = "hsa", out.suffix=paste0(i,".KEGGpath"))
  }}
  david.BP.table<-enrichDAVID(sigGenes,annotation="GOTERM_BP_FAT",david.user="elliot_akama-garren@hms.harvard.edu")
  if(!is.null(david.BP.table)){if(dim(david.BP.table)[1]!=0){
    david.BP.table<-setReadable(david.BP.table,OrgDb="org.Hs.eg.db",keyType="ENTREZID")
    write.csv(david.BP.table,file=paste0(i,".david.BP.table.csv"))
    barplot(david.BP.table, drop=TRUE, showCategory=12)+ggtitle(i)
    ggsave2(paste0(i,".david.BP.bar.png"),width=7, height=4,device="png")
    clusterProfiler::dotplot(david.BP.table)+ggtitle(i)
    ggsave2(paste0(i,".david.BP.dot.png"),width=7, height=4,device="png")
    emapplot(enrichplot::pairwise_termsim(david.BP.table))+ggtitle(i)
    ggsave2(paste0(i,".david.BP.emap.png"),width=10, height=8,device="png")
    clusterProfiler::cnetplot(david.BP.table, categorySize="pvalue", foldChange=geneList)+ggtitle(i)
    ggsave2(paste0(i,".david.BP.cnet.png"),width=10, height=8,device="png")
  }}
  #MsigDb
  egmt.table<-enricher(sigGenes,TERM2GENE = msig.df)
  if(!is.null(egmt.table)){if(dim(egmt.table)[1]!=0){
    egmt.table<-setReadable(egmt.table,OrgDb="org.Hs.eg.db",keyType="ENTREZID")
    write.csv(egmt.table,file=paste0(i,".egmt.table.csv"))
    barplot(egmt.table, drop=TRUE, showCategory=12)+ggtitle(i)
    ggsave2(paste0(i,".egmt.bar.png"),width=14, height=4,device="png")
    clusterProfiler::dotplot(egmt.table)+ggtitle(i)
    ggsave2(paste0(i,".egmt.dot.png"),width=14, height=4,device="png")
    emapplot(enrichplot::pairwise_termsim(egmt.table))+ggtitle(i)
    ggsave2(paste0(i,".egmt.emap.png"),width=10, height=8,device="png")
    clusterProfiler::cnetplot(egmt.table, categorySize="pvalue", foldChange=geneList)+ggtitle(i)
    ggsave2(paste0(i,".egmt.cnet.png"),width=10, height=8,device="png")
  }}
  egmt2.table<-GSEA(geneList,TERM2GENE = msig.df,pvalueCutoff=0.2)
  if(!is.null(egmt2.table)){if(dim(egmt2.table)[1]!=0){
    egmt2.table<-setReadable(egmt2.table,OrgDb="org.Hs.eg.db",keyType="ENTREZID")
    write.csv(egmt2.table,file=paste0(i,".egmt2.table.csv"))
    clusterProfiler::dotplot(egmt2.table)+ggtitle(i)
    ggsave2(paste0(i,".egmt2.dot.png"),width=14, height=4,device="png")
    emapplot(enrichplot::pairwise_termsim(egmt2.table))+ggtitle(i)
    ggsave2(paste0(i,".egmt2.emap.png"),width=10, height=8,device="png")
    clusterProfiler::cnetplot(egmt2.table, categorySize="pvalue", foldChange=geneList)+ggtitle(i)
    ggsave2(paste0(i,".egmt2.cnet.png"),width=10, height=8,device="png")
    clusterProfiler::gseaplot(egmt2.table,geneSetID=egmt2.table$ID[1],title=egmt2.table$Description[1])
    ggsave2(paste0(i,".egmt2.plot.png"),width=6, height=9,device="png")
  }}
}


#custom plot
tab<-gseKEGG(geneList,organism="hsa",pvalueCutoff=0.2) #KEGG gsea
df<-as.data.frame(tab)
df$lp<- -log10(df$p.adjust)
df<-df[order(df$lp,decreasing=T),]
df<-head(df,n=10)
ggplot(df,aes(x=reorder(Description,lp),y=lp))+geom_col()+coord_flip()+theme_classic()+
  theme(text=element_text(family="Arial"))+
  geom_hline(yintercept=1,linetype=2,size=0.2)+labs(x=NULL,y=NULL)
ggsave2("pdl1_gsea_table.png",width=5, height=3,device="png")



#custom volcano GSEA
obj2<-obj
genes.meta<-getBM(attributes=c("ensembl_gene_id", "mgi_symbol", "start_position", "end_position", "chromosome_name", 
                               "percentage_gene_gc_content", "external_gene_name", "gene_biotype","go_id","name_1006"),filters=
                    "mgi_symbol",values=list(rownames(obj@assays[["RNA"]]@meta.features)),
                  mart=useMart("ensembl", dataset = "mmusculus_gene_ensembl",host="ensembl.org"),useCache=F) #useast.

m_df = msigdbr(species = "Mus musculus")#, category = "C7") #H = hallmarks, C2=curated, C5=GO, C7=immune
pathways<-m_df %>% split(x = .$gene_symbol, f = .$gs_name)
msig.df<-m_df %>% dplyr::select(gs_name,entrez_gene)
m_df[m_df$gs_id=="M15399",]

for(i in c("M15399","M15816","M16940","M24299",
           "M2810","M14645","M40427","M13926","M16933",
           "GO:0072676","GO:0006090","GO:1903037","GO:0072678","GO:0006096","GO:0042110",
           "GO:0007159","GO:0040017","GO:0050852","GO:2000404","GO:0061615")){
  if(grepl("GO:",i)){gene.list2<-genes.meta$mgi_symbol[genes.meta$go_id==i]}else{
    gene.list2<-unique(m_df$gene_symbol[m_df$gs_id==i]) #M18810
  }
  if(length(gene.list2)>5){
    j="tfr_pd1.DE"
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
      geom_text_repel(color="black",size = 2, segment.size = 0.2,seed = 1, show.legend=T, box.padding=0.5) +
      geom_hline(yintercept=10,linetype=2,size=0.2)+
      geom_vline(xintercept=0.3,linetype=2,size=0.2)+
      geom_vline(xintercept= -0.3,linetype=2,size=0.2)+
      scale_size_manual(values=c(0.5,1))+
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
    ggsave2(paste0(j,"_volcano_gsea",i,".png"),width=1.7, height=1.5,device="png")
    
    gene.list<-df$gene[df$gene %in% gene.list2 & (df$rank_diff<100 | df$rank_sum <100 |df$rank_sum_dec<100)]
    df <- df %>% mutate(label = ifelse(gene %in% gene.list2, "sig", "notsig")) #genes to color
    df<-df[order(df$label),]
    df<- df %>% mutate(gene_label = ifelse(gene %in% gene.list , gene, NA))
    ggplot(df, aes(b6FC,m564FC,colour=label,size=label,label=gene_label)) + 
      geom_point()+#geom_point(fill=NA,colour=alpha("black",0.5),pch=21,size=3) +
      theme_classic()+
      scale_size_manual(values=c(0.5,1))+
      geom_text_repel(color="black",size = 2, segment.size = 0.2,seed = 1, show.legend=T, box.padding=0.5) +
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




