###Pathway analysis
library(org.Hs.eg.db)
library(tidyverse)
library(RDAVIDWebService)
library(Seurat)
library(cowplot)
library(fgsea)
library(clusterProfiler)
library(pathview)
library(topGO)
library(scde)
library(biomaRt)
library(GO.db)
library(DBI)
library(VISION)
library(msigdbr)
library(DOSE)

m_df = msigdbr(species = "Homo sapiens")#, category = "C7") #H = hallmarks, C2=curated, C5=GO, C7=immune
pathways<-m_df %>% split(x = .$gene_symbol, f = .$gs_name) 
for(i in c(ls(pattern="DE"),ls(pattern=".markers"),ls(pattern=".response"))){
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
    topUp <- head(fgseaRes %>% filter(ES > 0) %>% top_n(10, wt=-padj),n=10)
    topDown <- head(fgseaRes %>% filter(ES < 0) %>% top_n(10, wt=-padj),n=10)
    topPathways <- bind_rows(topUp, topDown) %>% arrange(-ES)
    png(paste0(i,".gsea.table.png"),width=20,height=7,units="in",res=200)
    plotGseaTable(pathways[topPathways$pathway], ranks, fgseaRes, gseaParam = 0.5)
    dev.off()
    write.csv(fgseaRes[order(padj, -abs(NES)), ][,1:7],file=paste0(i,".gsea.results.csv"))
    write.csv(topPathways[,1:7],file=paste0(i,".gsea.table.csv"))
  }}
}

##clusterProfiler (https://bioconductor.statistik.tu-dortmund.de/packages/3.6/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html)
m_df = msigdbr(species = "Homo sapiens")#, category = "C7") #H = hallmarks, C2=curated, C5=GO, C7=immune
msig.df<-m_df %>% dplyr::select(gs_name,entrez_gene)
for(i in c(ls(pattern="DE"),ls(pattern=".markers"),ls(pattern=".response"))){
  df<-get(i)
  df$gene<-rownames(df)
  genes.meta<-bitr(df$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") #pull in entrezID, could also use bitr_kegg, see: keytypes(org.Hs.eg.db)
  df<-cbind(df,genes.meta[match(df$gene, genes.meta$SYMBOL),])
  sigGenes <- na.exclude(df$ENTREZID[df$p_val_adj < 0.01 &df$avg_log2FC > 1 ]) #select sig genes, abs(df$avg_log2FC) 
  df2<-df[df$p_val_adj<0.01 & !is.na(df$ENTREZID),]
  ranks <- df2$avg_log2FC
  names(ranks) <- df2$ENTREZID
  geneList <- sort(ranks, decreasing = TRUE)
  #topGO
  ggo.table<-groupGO(gene=sigGenes,OrgDb=org.Hs.eg.db,ont="MF",level=4,readable=T) #BP (biological process), MF (molecular function), CC (cellular componnent)
  if(!is.null(ggo.table)){if(dim(ggo.table)[1]!=0){
    write.csv(ggo.table,file=paste0(i,".gGO.table.csv"))
    barplot(ggo.table, drop=TRUE, showCategory=12)+ggtitle(i)
    ggsave2(paste0(i,".gGO.bar.png"),width=8, height=4,device="png")
  }}
  ego.table<-enrichGO(gene=sigGenes,OrgDb=org.Hs.eg.db,ont="MF", readable=T)
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
  gsea.table<-gseGO(geneList,OrgDb=org.Hs.eg.db,ont="MF",pvalueCutoff = 0.2)
  if(!is.null(gsea.table)){if(dim(gsea.table)[1]!=0){
    gsea.table<-setReadable(gsea.table,OrgDb=org.Hs.eg.db)
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
    kk.table <- setReadable(kk.table,OrgDb=org.Hs.eg.db,keyType="ENTREZID")
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
    kk2.table<-setReadable(kk2.table,OrgDb=org.Hs.eg.db,keyType="ENTREZID") #KEGG gsea
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
    mkk.table<-setReadable(mkk.table,OrgDb=org.Hs.eg.db,keyType="ENTREZID") #KEGG Module
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
    mkk2.table<-setReadable(mkk2.table,OrgDb=org.Hs.eg.db,keyType="ENTREZID") #KEGG Module gsea
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
  david.KEGG.table<-enrichDAVID(sigGenes,annotation="KEGG_PATHWAY",david.user="***") 
  if(!is.null(david.KEGG.table)){if(dim(david.KEGG.table)[1]!=0){
    david.KEGG.table<-setReadable(david.KEGG.table,OrgDb=org.Hs.eg.db,keyType="ENTREZID") 
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
  david.BP.table<-enrichDAVID(sigGenes,annotation="GOTERM_BP_FAT",david.user="***")
  if(!is.null(david.BP.table)){if(dim(david.BP.table)[1]!=0){
    david.BP.table<-setReadable(david.BP.table,OrgDb=org.Hs.eg.db,keyType="ENTREZID")
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
    egmt.table<-setReadable(egmt.table,OrgDb=org.Hs.eg.db,keyType="ENTREZID")
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
    egmt2.table<-setReadable(egmt2.table,OrgDb=org.Hs.eg.db,keyType="ENTREZID")
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
df<-as.data.frame(david.BP.table)
df$lp<- -log10(df$p.adjust)
df<-df[order(df$lp,decreasing=T),]
# df<-df[!grepl("process|resorption|lipid|neuron|myelination|axon",df$Description),]
df<-df[grepl("T cell|B cell|adhesion|endocytosis|activation|integrin-mediated|activation",df$Description),]
df<-df[!grepl("positive|plasma",df$Description),]
ggplot(df,aes(x=reorder(Description,lp),y=lp))+geom_col()+coord_flip()+theme_classic()+
  theme(text=element_text(family="Arial"))+
  geom_hline(yintercept=1,linetype=2,size=0.2)+labs(x=NULL,y=NULL)
ggsave2("david_kegg_table.png",width=4, height=3,device="png")

gene.list2<-unique(unlist(str_split(df$geneID[df$Description=="regulation of T cell activation"],pattern="/")))
gene.list2<-unique(unlist(str_split(df$geneID[grepl("T cell",df$Description)],pattern="/")))
gene.list2<-unique(m_df$gene_symbol[m_df$gs_id=="M16253"]) #M18810
gene.list2

#Cluster GO dotplot
load("SLO_all.RData")
Idents(SLO_all.combined)<-"my.clusters"
clust.mark<- FindAllMarkers(object = SLO_all.combined, test.use = "MAST",only.pos = TRUE)
genes.meta<-bitr(clust.mark$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") #pull in entrezID, could also use bitr_kegg, see: keytypes(org.Hs.eg.db)
clust.mark<-cbind(clust.mark,genes.meta[match(clust.mark$gene, genes.meta$SYMBOL),])
clust.mark<-clust.mark[!is.na(clust.mark$ENTREZID),]
clust.mark[] <- lapply(clust.mark, as.character)
clust.mark.genes<-clust.mark %>% split(x = .$ENTREZID, f = .$cluster) 
ck<-compareCluster(geneCluster=clust.mark.genes,fun="enrichKEGG",organism="hsa") 
clusterProfiler::dotplot(ck)
ggsave2("cluster.KEGG.dot.png",width=10, height=8,device="png")
c.ggo<-compareCluster(geneCluster=clust.mark.genes,fun="groupGO",OrgDb=org.Hs.eg.db,ont="BP",level=4,readable=T) 
clusterProfiler::dotplot(c.ggo)
ggsave2("cluster.gGO.dot.png",width=10, height=8,device="png")
c.ego<-compareCluster(geneCluster=clust.mark.genes,fun="enrichGO",OrgDb=org.Hs.eg.db,ont="BP", readable=T)
clusterProfiler::dotplot(c.ego)
ggsave2("cluster.eGO.dot.png",width=10, height=8,device="png")
c.msigdb<-compareCluster(geneCluster=clust.mark.genes,fun="enricher",TERM2GENE = msig.df)
clusterProfiler::dotplot(c.msigdb)
ggsave2("cluster.msigdb.dot.png",width=14, height=8,device="png")

