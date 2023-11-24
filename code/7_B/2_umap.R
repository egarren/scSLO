rm(list=ls())
library(Seurat)
library(cowplot)
library(ggplot2)
library(plyr)
library(dplyr)
library(biomaRt)
library(plotly)
library(scales)
library(EnhancedVolcano)
library(data.table)
library(ggpubr)
library(limma)
library(VennDiagram)
library(viridis)
library(pheatmap)
library(phylotools) 
library(ggforce)
library(harmony)
library(readxl)

obj<-B_all


##QC and Filter
#add gene metadata
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
obj@assays[["RNA"]]@meta.features$original_hgnc<-rownames(obj@assays[["RNA"]]@meta.features)
obj <- subset(obj, features=rownames(obj[!(grepl("NA..",obj@assays[["RNA"]]@meta.features$hgnc_symbol)),]))
obj <- subset(obj, features=rownames(obj[!(grepl("X..",obj@assays[["RNA"]]@meta.features$hgnc_symbol)),]))
genes.meta<-getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "external_gene_name", "gene_biotype","go_id","name_1006"),filters=
                    "hgnc_symbol",values=list(rownames(obj@assays[["RNA"]]@meta.features)), mart=human,useCache=F) #useast.
m <- match(obj@assays[["RNA"]]@meta.features$original_hgnc, genes.meta$hgnc_symbol)
obj@assays[["RNA"]]@meta.features<-cbind(obj@assays[["RNA"]]@meta.features,genes.meta[m,])
rownames(obj@assays[["RNA"]]@meta.features) <- make.names(obj@assays[["RNA"]]@meta.features$hgnc_symbol,unique=T)
rownames(obj@assays[["RNA"]]@data) <- make.names(obj@assays[["RNA"]]@meta.features$hgnc_symbol,unique=T)
rownames(obj@assays[["RNA"]]@counts) <- make.names(obj@assays[["RNA"]]@meta.features$hgnc_symbol,unique=T)
save.image("temp.merged2.RData")

# #TCR and Ig identification and filtering
ig_list <- c("IG_C_pseudogene", "IG_D_gene", "IG_D_pseudogene", "IG_J_gene", "IG_LV_gene", 
             "IG_pseudogene", "IG_V_gene", "IG_V_pseudogene")#"IG_C_gene", 
tr_list <-c("TR_V_gene", "TR_V_pseudogene", "TR_D_gene", "TR_J_gene", "TR_J_pseudogene", "TR_C_gene")
obj <- subset(obj, features=rownames(obj[!(obj@assays[["RNA"]]@meta.features$gene_biotype %in% ig_list),]))
obj <- subset(obj, features=rownames(obj[!(obj@assays[["RNA"]]@meta.features$gene_biotype %in% tr_list),]))
obj <- subset(obj, features=rownames(obj[!(is.na(obj@assays[["RNA"]]@meta.features$hgnc_symbol)),]))
# mitochondrial DNA percentage
obj <- PercentageFeatureSet(obj, pattern = "^MT", col.name = "percent.mt")
obj <- PercentageFeatureSet(obj, pattern = "^RPL", col.name = "percent.Rpl")
obj <- PercentageFeatureSet(obj, pattern = "^RPS", col.name = "percent.Rps")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
obj$CC.Difference <- obj$S.Score - obj$G2M.Score
#HSP regression
HSP_genes<-genes.meta[genes.meta$go_id=="GO:0034605",]$mgi_symbol
# HSP_genes<-genes.meta[genes.meta$go_id=="GO:0034605",]$hgnc_symbol
obj <- AddModuleScore(object = obj,features = list(HSP_genes),name = 'HSP.score')
Idents(obj)<-"GSE"
VlnPlot(object = obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",
                                     "percent.Rpl","percent.Rps","HSP.score1"),pt.size=0, ncol = 3)
ggsave2("vln.QC.png",device="png")
save.image("temp.preQC.RData")
#filter
dim(obj)
hist(log10(obj@meta.data$percent.Rpl),breaks=100)
hist(obj@meta.data$percent.Rps,breaks=100)
hist(obj@meta.data$percent.mt[obj@meta.data$percent.mt<15],breaks=100)
obj<- subset(x = obj, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & nCount_RNA>300 & 
                  percent.mt >  -Inf & percent.mt < 6 & percent.Rpl < 20 & percent.Rps < 20 ) 
dim(obj)
rm(list=setdiff(ls(), c("obj", "metadata","genes.meta")))
save.image("temp.QCd.RData")

#Harmony integration
obj <- obj %>% Seurat::NormalizeData(verbose = FALSE) %>%
  ScaleData(vars.to.regress=c("nCount_RNA", "percent.mt","HSP.score1","percent.Rpl","percent.Rps","nFeature_RNA"),verbose = FALSE) %>%
  RunPCA(npcs = 20, verbose = FALSE)
colnames(obj@meta.data)
as.data.frame(table(obj@meta.data$Sample))
obj.combined<-RunHarmony(obj,c("Species"),plot_convergence=T)
obj.combined <- RunUMAP(obj.combined, reduction = "harmony", dims = 1:20) #change based on elbow plot
obj.combined <- FindNeighbors(obj.combined, reduction = "harmony", dims = 1:20) #change based on elbow plot
obj.combined <- FindClusters(obj.combined, resolution = 0.2) #adjust resolution (bigger=more clusters), initially used 0.2
DimPlot(obj.combined, reduction = "umap", split.by="GSE")
DimPlot(obj.combined, reduction = "umap")
for (i in c("Species","Organ","Tx","seurat_clusters","GSE","Sample","Chemistry")){
  Idents(obj.combined)<-i
  DimPlot(obj.combined, reduction = "umap",pt.size=0.1)+NoLegend()
  ggsave2(paste0(i,".umap.png"),width=6, height=5,device="png")
}
for (i in c("nCount_RNA","nFeature_RNA","percent.mt","percent.Rps","percent.Rpl","CC.Difference","HSP.score1","S.Score","G2M.Score")){
  FeaturePlot(obj.combined, features= i,split.by = "Species",pt.size=0.1, order=T)
  ggsave2(paste0(i,".umap.png"),width=10, height=5,device="png")
}


rm(list=setdiff(ls(), c("obj.combined", "metadata","genes.meta")))
save.image("merged.RData")


## Cluster Analysis
Idents(obj.combined)<-"seurat_clusters"
DimPlot(obj.combined, reduction = "umap",pt.size=0.1,label=T)+NoLegend()
ggsave2("umap.label.png",width=6, height=5,device="png")
#https://link.springer.com/chapter/10.1007/978-3-030-86016-5_5/tables/1
FeaturePlot(object = obj.combined, features = c("IGHD","CD21","ITGAX","TLR7","TRAF5"))
FeaturePlot(object = obj.combined, features = "CD93")#, cols = c("grey", "blue"), reduction = "umap",pt.size=0.001, order=T)
FeaturePlot(object = obj.combined, features = "CD274",split.by="Species")#, cols = c("grey", "blue"), reduction = "umap",pt.size=0.001, order=T)
gene.list<-c("IGHD","SELL","FCMR","CD27","TNFRSF13B","FCRL5","MKI67","STMN1","BCL2A1","JCHAIN","XBP1","PRDM1","S1PR1",
             "MYC","IRF4","BATF","PAX5","FCER2","CXCR4","CD86","CD83","SDC1","CD38","AICDA","PLD4","MZB1","CD1C","SOX4",
             "ITGAX","TBX21","CR2","IL4R","CCR6","MS4A1","ICOSLG","IL21R","CXCR5","IGHM","NOTCH2","BCL6",
             "ITGAM","GPR138","TNFRSF13B","CD27","CXCR3","ITGB2","VIM","PLAC8","CD38","APOE")
FeaturePlot(object = obj.combined, features = gene.list, 
               cols = c("grey", "blue"), reduction = "umap",pt.size=0.001, order=T)
ggsave2("umap.goi.png",width=10, height=10,device="png")
VlnPlot(obj.combined, features = "CD1D",pt.size = 0.001)
# VlnPlot(obj.combined, features = "Cd274",pt.size = 0.001)
VlnPlot(obj.combined, features = "CD274",pt.size = 0.001)
ggsave2("vln.cd274.png",width=4, height=3,device="png")

##cluster DE
obj.markers <- FindAllMarkers(object = obj.combined, test.use = "MAST",max.cells.per.ident=1000)
write.csv(obj.markers,"cluster.markers.csv")
write.csv(obj.markers %>% group_by(cluster) %>% top_n(20, avg_log2FC),"cluster.markers.top20.csv")
top6 <- obj.markers %>% group_by(cluster) %>% top_n(6, avg_log2FC)
DoHeatmap(object = obj.combined, features = top6$gene, label = TRUE)  #slim.col.label to TRUE prints cluster IDs instead of cells
ggsave2("cluster.heatmap.png",width=10, height=10,device="png")
obj.combined@meta.data$my.clusters <- Idents(obj.combined)  # Store cluster identities in object@meta.data$my.clusters
save.image("temp.analyzed1.RData")

#rename clusters
new.cluster.ids <- c("Naive","GC_1","MBC_1","DZ","ASC_1","Activated","PC","LZ","MBC_2","GC_2","ASC_2","MZ_1","ASC_3","MZ_2",
                     "","","","","")
names(new.cluster.ids) <- levels(obj.combined)
obj.combined <- RenameIdents(obj.combined, new.cluster.ids)
obj.combined@meta.data$my.clusters2 <- Idents(obj.combined) 
cluster.select <- c("Naive","GC_1","MBC_1","DZ","ASC_1","Activated","PC","LZ","MBC_2","GC_2","ASC_2","MZ_1","ASC_3","MZ_2")
obj.combined<-subset(obj.combined,idents=cluster.select)
obj.combined@meta.data$my.clusters2  <- factor(obj.combined@meta.data$my.clusters2,levels = cluster.select)
DimPlot(obj.combined, reduction = "umap",pt.size=0.1,label=T)+NoLegend()
ggsave2("umap.label2.png",width=5, height=5,device="png")
dim(obj.combined)
save.image("temp.analyzed2.RData")

obj2<-obj.combined
obj2[["SCT"]]<-NULL
obj2[["integrated"]]<-NULL
obj.sce <- as.SingleCellExperiment(obj2)
save(obj.sce,file="sce.filtered.RData")
show_col(hue_pal()(14))


#export mouse data for RNA velocity
obj2<-obj.combined
#mgi to human gene name
obj2@assays[["RNA"]]@meta.features$hgnc_symbol<-rownames(obj2@assays[["RNA"]]@meta.features)
hgnc_list<-getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = list(rownames(obj2@assays[["RNA"]]@meta.features)) , 
                  mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl",host="dec2021.archive.ensembl.org"), attributesL = c("mgi_symbol"), 
                  martL = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host="dec2021.archive.ensembl.org"), uniqueRows=T)
m <- match(obj2@assays[["RNA"]]@meta.features$hgnc_symbol, hgnc_list$HGNC.symbol)
obj2@assays[["RNA"]]@meta.features<-cbind(obj2@assays[["RNA"]]@meta.features,hgnc_list[m,])
obj2@assays[["RNA"]]@meta.features$mgi_symbol<-obj2@assays[["RNA"]]@meta.features$MGI.symbol
rownames(obj2@assays[["RNA"]]@meta.features) <- make.names(obj2@assays[["RNA"]]@meta.features$mgi_symbol,unique=T)
rownames(obj2@assays[["RNA"]]@data) <- make.names(obj2@assays[["RNA"]]@meta.features$mgi_symbol,unique=T)
rownames(obj2@assays[["RNA"]]@counts) <- make.names(obj2@assays[["RNA"]]@meta.features$mgi_symbol,unique=T)
obj2 <- subset(obj2, features=rownames(obj2[!(is.na(obj2@assays[["RNA"]]@meta.features$mgi_symbol)),]))
obj2 <- subset(obj2, features=rownames(obj2[obj2@assays[["RNA"]]@meta.features$mgi_symbol != "",]))
obj.sce <- as.SingleCellExperiment(obj2,assay="RNA")
save(obj.sce,file="sce.ms.RData")

#export mouse data for RNA velocity
Idents(obj.combined)<-"GSE"
obj2<-subset(obj.combined,idents=c("GSE203132","GSE192762"))
# save metadata table:
obj2$barcode <- colnames(obj2)
obj2$UMAP_1 <- obj2@reductions$umap@cell.embeddings[,1]
obj2$UMAP_2 <- obj2@reductions$umap@cell.embeddings[,2]
write.csv(obj2@meta.data, file='./velocity/metadata.csv', quote=F, row.names=F)
# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(obj2, assay='RNA', slot='counts')
writeMM(counts_matrix, file='./velocity/counts.mtx')
# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(obj2@reductions$pca@cell.embeddings, file='./velocity/pca.csv', quote=F, row.names=F)
# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='./velocity/gene_names.csv',
  quote=F,row.names=F,col.names=F
)






##cluster DE2
obj.markers <- FindAllMarkers(object = obj.combined, test.use = "MAST",max.cells.per.ident=1000)
write.csv(obj.markers,"cluster.markers2.csv")
write.csv(obj.markers %>% group_by(cluster) %>% top_n(20, avg_log2FC),"cluster.markers.top20.csv")
top6 <- obj.markers %>% group_by(cluster) %>% top_n(6, avg_log2FC)
DoHeatmap(object = obj.combined, features = top6$gene, label = TRUE)  #slim.col.label to TRUE prints cluster IDs instead of cells
ggsave2("cluster.heatmap2.png",width=10, height=10,device="png")
obj.combined@meta.data$my.clusters <- Idents(obj.combined)  # Store cluster identities in object@meta.data$my.clusters
save.image("temp.analyzed3.RData")


##cluster DE
obj<-obj.combined
Idents(obj) <- "my.clusters2" #setting idents to condition metadata
lz_dz.DE <- FindMarkers(obj, ident.1 = "LZ", ident.2 = "DZ", min.pct=0,logfc.threshold = -Inf,test.use = "MAST",max.cells.per.ident=10000)
write.csv(lz_dz.DE ,"lz_dz_DE.csv")

asc2_dz.DE <- FindMarkers(obj, ident.1 = "ASC_2", ident.2 = "DZ", min.pct=0,logfc.threshold = -Inf,test.use = "MAST",max.cells.per.ident=10000)
write.csv(asc2_dz.DE ,"asc2_dz_DE.csv")

new.cluster.ids <- c("Naive","GC","MBC","DZ","ASC","Activated","ASC","LZ","MBC","GC","ASC","MZ","ASC","MZ")
names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj, new.cluster.ids)
obj@meta.data$my.clusters3 <- Idents(obj) 
asc_gc.DE <- FindMarkers(obj, ident.1 = "ASC", ident.2 = "GC", min.pct=0,logfc.threshold = -Inf,test.use = "MAST",max.cells.per.ident=10000)
write.csv(asc_gc.DE ,"asc_gc.csv")

#PDL1 DE
test<-FetchData(obj,vars=c("CD274"))
hist(test$CD274,breaks=100)
obj@meta.data$my.clusters4<-"pdl1_lo"
obj@meta.data$my.clusters4[test$CD274==0]<-"pdl1_neg"
obj@meta.data$my.clusters4[test$CD274>0.5]<-"pdl1_hi"
Idents(obj) <- "my.clusters4" #setting idents to condition metadata
pdl1.DE <- FindMarkers(obj, ident.1 = "pdl1_hi", ident.2 = "pdl1_neg", min.pct=0,logfc.threshold = -Inf,test.use = "MAST",max.cells.per.ident=10000)
write.csv(pdl1.DE ,"pdl1_vs.csv")



Idents(obj)<-"my.clusters3"
Naive<-subset(obj,idents="Naive")
Idents(Naive) <- "my.clusters4" #setting idents to condition metadata
Naive_pdl1.DE <- FindMarkers(Naive, ident.1 = "pdl1_hi", ident.2 = "pdl1_neg", min.pct=0,logfc.threshold = -Inf,test.use = "MAST")
write.csv(Naive_pdl1.DE ,"Naive_pdl1.csv")
GC<-subset(obj,idents="GC")
Idents(GC) <- "my.clusters4" #setting idents to condition metadata
GC_pdl1.DE <- FindMarkers(GC, ident.1 = "pdl1_hi", ident.2 = "pdl1_neg", min.pct=0,logfc.threshold = -Inf,test.use = "MAST")
write.csv(GC_pdl1.DE ,"GC_pdl1.csv")
ASC<-subset(obj,idents="ASC")
Idents(ASC) <- "my.clusters4" #setting idents to condition metadata
ASC_pdl1.DE <- FindMarkers(ASC, ident.1 = "pdl1_hi", ident.2 = "pdl1_neg", min.pct=0,logfc.threshold = -Inf,test.use = "MAST")
write.csv(ASC_pdl1.DE ,"ASC_pdl1.csv")

save.image("temp.analyzed1.RData")

#PD1 DE
test<-FetchData(obj,vars=c("PDCD1"))
hist(test$PDCD1,breaks=100)
obj@meta.data$my.clusters4<-"pdl1_mid"
obj@meta.data$my.clusters4[test$PDCD1==0]<-"pd1_neg"
obj@meta.data$my.clusters4[test$PDCD1>0.5]<-"pd1_hi"
Idents(obj) <- "my.clusters4" #setting idents to condition metadata
pd1.DE <- FindMarkers(obj, ident.1 = "pd1_hi", ident.2 = "pd1_neg", min.pct=0,logfc.threshold = -Inf,test.use = "MAST")
write.csv(pd1.DE,"pd1_DE.csv")

FeaturePlot(obj, "PDCD1") & scale_color_gradientn(colors = plasma(n = 10, direction = -1))
ggsave2("pd1_umap.png",width=4, height=4,device="png")

Idents(obj)<-"my.clusters3"
Naive<-subset(obj,idents="Naive")
Idents(Naive) <- "my.clusters4" #setting idents to condition metadata
Naive_pd1.DE <- FindMarkers(Naive, ident.1 = "pd1_hi", ident.2 = "pd1_neg", min.pct=0,logfc.threshold = -Inf,test.use = "MAST")
write.csv(Naive_pd1.DE ,"Naive_pd1.csv")
GC<-subset(obj,idents="GC")
Idents(GC) <- "my.clusters4" #setting idents to condition metadata
GC_pd1.DE <- FindMarkers(GC, ident.1 = "pd1_hi", ident.2 = "pd1_neg", min.pct=0,logfc.threshold = -Inf,test.use = "MAST")
write.csv(GC_pd1.DE ,"GC_pd1.csv")
ASC<-subset(obj,idents="ASC")
Idents(ASC) <- "my.clusters4" #setting idents to condition metadata
ASC_pd1.DE <- FindMarkers(ASC, ident.1 = "pd1_hi", ident.2 = "pd1_neg", min.pct=0,logfc.threshold = -Inf,test.use = "MAST")
write.csv(ASC_pd1.DE ,"ASC_pd1.csv")



brewer.pal(n = 2, name = "Dark2")

save.image("analyzed.RData")
