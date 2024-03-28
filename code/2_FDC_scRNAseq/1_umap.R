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

setwd("/scFDC")

load("temp.merged.RData")

df<-SLO_all@meta.data
df<-metadata
unique(df$Organ[df$Species=="Mus"])
unique(df$Organ[df$Species=="Hs"])

table(df$Organ2[df$Species=="Mus"])
table(df$Organ2[df$Species=="Hs"])

unique(df$Tx2[df$Species=="Mus"])
unique(df$Tx2[df$Species=="Hs"])

length(unique(df$GSE[df$Species=="Mus"]))
length(unique(df$GSE[df$Species=="Hs"]))

length(unique(df$Sample[df$Species=="Mus"]))
length(unique(df$Sample[df$Species=="Hs"]))

##QC and Filter
#add gene metadata
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
SLO_all@assays[["RNA"]]@meta.features$original_hgnc<-rownames(SLO_all@assays[["RNA"]]@meta.features)
SLO_all <- subset(SLO_all, features=rownames(SLO_all[!(grepl("NA..",SLO_all@assays[["RNA"]]@meta.features$hgnc_symbol)),]))
genes.meta<-getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "external_gene_name", "gene_biotype","go_id","name_1006"),filters=
                    "hgnc_symbol",values=list(rownames(SLO_all@assays[["RNA"]]@meta.features)), mart=human,useCache=F) #useast.
m <- match(SLO_all@assays[["RNA"]]@meta.features$original_hgnc, genes.meta$hgnc_symbol)
SLO_all@assays[["RNA"]]@meta.features<-cbind(SLO_all@assays[["RNA"]]@meta.features,genes.meta[m,])
rownames(SLO_all@assays[["RNA"]]@meta.features) <- make.names(SLO_all@assays[["RNA"]]@meta.features$hgnc_symbol,unique=T)
rownames(SLO_all@assays[["RNA"]]@data) <- make.names(SLO_all@assays[["RNA"]]@meta.features$hgnc_symbol,unique=T)
rownames(SLO_all@assays[["RNA"]]@counts) <- make.names(SLO_all@assays[["RNA"]]@meta.features$hgnc_symbol,unique=T)

# #TCR and Ig identification and filtering
ig_list <- c("IG_C_gene", "IG_C_pseudogene", "IG_D_gene", "IG_D_pseudogene", "IG_J_gene", "IG_LV_gene", 
             "IG_pseudogene", "IG_V_gene", "IG_V_pseudogene")
tr_list <-c("TR_V_gene", "TR_V_pseudogene", "TR_D_gene", "TR_J_gene", "TR_J_pseudogene", "TR_C_gene")
SLO_all <- subset(SLO_all, features=rownames(SLO_all[!(SLO_all@assays[["RNA"]]@meta.features$gene_biotype %in% ig_list),]))
SLO_all <- subset(SLO_all, features=rownames(SLO_all[!(SLO_all@assays[["RNA"]]@meta.features$gene_biotype %in% tr_list),]))
SLO_all <- subset(SLO_all, features=rownames(SLO_all[!(is.na(SLO_all@assays[["RNA"]]@meta.features$hgnc_symbol)),]))
# mitochondrial DNA percentage
SLO_all <- PercentageFeatureSet(SLO_all, pattern = "^MT", col.name = "percent.mt")
SLO_all <- PercentageFeatureSet(SLO_all, pattern = "^RPL", col.name = "percent.Rpl")
SLO_all <- PercentageFeatureSet(SLO_all, pattern = "^RPS", col.name = "percent.Rps")
#cell cycle regression
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
SLO_all <- CellCycleScoring(SLO_all, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
SLO_all$CC.Difference <- SLO_all$S.Score - SLO_all$G2M.Score
#HSP regression
HSP_genes<-genes.meta[genes.meta$go_id=="GO:0034605",]$hgnc_symbol
SLO_all <- AddModuleScore(object = SLO_all,features = list(HSP_genes),name = 'HSP.score')
Idents(SLO_all)<-"GSE"
VlnPlot(object = SLO_all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",
                                     "percent.Rpl","percent.Rps","HSP.score1"),pt.size=0, ncol = 3)
ggsave2("vln.QC.png",device="png")
save.image("temp.preQC.RData")
#filter
dim(SLO_all)
save<-SLO_all
hist(log10(SLO_all@meta.data$S.Score),breaks=100)
SLO_all<- subset(x = SLO_all, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & nCount_RNA>300 & 
                  percent.mt >  -Inf & percent.mt < 7 & percent.Rpl < 20 & percent.Rps < 20 &
                 S.Score <0.5 & G2M.Score<0.5) 
dim(SLO_all)
SLO_all.pp.sce <- as.SingleCellExperiment(SLO_all)
save(SLO_all.pp.sce,file="pp.sce.filtered.RData")
rm(list=setdiff(ls(), c("SLO_all", "metadata","genes.meta")))
save.image("temp.QCd.RData")

#Harmony integration
SLO_all <- SLO_all %>% Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(vars.to.regress=c("nCount_RNA", "percent.mt","HSP.score1","percent.Rpl","percent.Rps"),verbose = FALSE) %>% 
  RunPCA(pc.genes = SLO_all@var.genes, npcs = 20, verbose = FALSE)
SLO_all.combined<-RunHarmony(SLO_all,c("GSE","Organ","Species"),plot_convergence=T)
SLO_all.combined <- RunUMAP(SLO_all.combined, reduction = "harmony", dims = 1:20) #change based on elbow plot
SLO_all.combined <- FindNeighbors(SLO_all.combined, reduction = "harmony", dims = 1:20) #change based on elbow plot
SLO_all.combined <- FindClusters(SLO_all.combined, resolution = 0.15) #adjust resolution (bigger=more clusters), initially used 0.2
DimPlot(SLO_all.combined, reduction = "umap")
for (i in c("Species","Organ2","Tx2","seurat_clusters","GSE","Sample","Chemistry")){
  Idents(SLO_all.combined)<-i
  DimPlot(SLO_all.combined, reduction = "umap",pt.size=0.1)
  ggsave2(paste0(i,".umap.png"),width=6, height=5,device="png")
}
for (i in c("nCount_RNA","nFeature_RNA","percent.mt","percent.Rps","percent.Rpl","CC.Difference","HSP.score1","S.Score","G2M.Score")){
  FeaturePlot(SLO_all.combined, features= i,split.by = "Species",pt.size=0.1, order=T)
  ggsave2(paste0(i,".umap.png"),width=10, height=5,device="png")
}
rm(list=setdiff(ls(), c("SLO_all.combined", "metadata","genes.meta")))
save.image("merged.RData")

#Add cyster metadata
Cyster.metadata<-read.csv("/Cyster.metadata.csv", header=T) #define path
cyster_cellID<-gsub("SRR6976738","1",SLO_all.combined@meta.data$cellID)
cyster_cellID<-gsub("SRR6976741","2",cyster_cellID)
cyster_cellID<-substr(cyster_cellID,1,nchar(cyster_cellID)-4)
SLO_all.combined@meta.data$Cyster.cellID<-cyster_cellID
SLO_all.combined@meta.data<-left_join(x = SLO_all.combined@meta.data, y = Cyster.metadata[,c("Barcode","CC.label")], 
                                      by = c("Cyster.cellID"="Barcode"),keep=F) #add Cyster metadata (clusters)
rownames(SLO_all.combined@meta.data)<-SLO_all.combined@meta.data$cellID
save.image("merged2.RData")


## Cluster Analysis
Idents(SLO_all.combined)<-"seurat_clusters"
DimPlot(SLO_all.combined, reduction = "umap",group.by="CC.label",pt.size=0.1,
        cells=SLO_all.combined@meta.data$cellID[!is.na(SLO_all.combined@meta.data$CC.label)])
ggsave2("cyster.umap.png",width=6, height=5,device="png")
DimPlot(SLO_all.combined, reduction = "umap",pt.size=0.1,label=T)+NoLegend()
ggsave2("umap.label.png",width=6, height=5,device="png")
FeaturePlot(object = SLO_all.combined, features = "CD68")#, cols = c("grey", "blue"), reduction = "umap",pt.size=0.001, order=T)
FeaturePlot(object = SLO_all.combined, features = "CD274")#, cols = c("grey", "blue"), reduction = "umap",pt.size=0.001, order=T)
gene.list<-c("MADCAM1","TNFSF11","CXCL13","CR2","CXCL12","CCL19","CCL21","PDPN","CXCL9","CXCL10",
             "CH25H","IL7","LEPR","NR4A1","LEPR","CD34","TNFSF13B","INMT","CD274",
             "CD4","CD8A","CD3E","CD19","PTPRC","PECAM1","ACTA2","ITGAX","ITGAM","CD74",
             "CXCR5","SDC1","MARCO","PTX3","SIGLEC1","PROX1","DCN","VWF","TNFSF13B","VCAM1","MYH11",
             "APOE","AGT","FBN1","PTX3","ATF3","CD5L","EPCAM")
FeaturePlot(object = SLO_all.combined, features = gene.list, 
               cols = c("grey", "blue"), reduction = "umap",pt.size=0.001, order=T)
ggsave2("umap.goi.png",width=10, height=20,device="png")
VlnPlot(SLO_all.combined, features = "CD274",pt.size = 0.001, combine=F,raster=F)
ggsave2("vln.cd274.png",width=4, height=3,device="png")

##cluster DE
# SLO_all.markers <- FindAllMarkers(object = SLO_all.combined, test.use = "MAST")
SLO_all.markers <- FindAllMarkers(object = SLO_all.combined, test.use = "MAST",max.cells.per.ident=1000)
write.csv(SLO_all.markers,"cluster.markers.csv")
write.csv(SLO_all.markers %>% group_by(cluster) %>% top_n(20, avg_log2FC),"cluster.markers.top20.csv")
top6 <- SLO_all.markers %>% group_by(cluster) %>% top_n(6, avg_log2FC)
DoHeatmap(object = SLO_all.combined, features = top6$gene, label = TRUE)  #slim.col.label to TRUE prints cluster IDs instead of cells
ggsave2("cluster.heatmap.png",width=10, height=10,device="png")
SLO_all.combined@meta.data$my.clusters <- Idents(SLO_all.combined)  # Store cluster identities in object@meta.data$my.clusters
save.image("temp.analyzed1.RData")

#rename clusters
new.cluster.ids <- c("MRC","Adventitial","BEC","FDC_TRC","B","LEC","Myeloid","TRC","T","Epithelial","Periventricular","EPSTI1","","","","")
names(new.cluster.ids) <- levels(SLO_all.combined)
SLO_all.combined <- RenameIdents(SLO_all.combined, new.cluster.ids)
SLO_all.combined@meta.data$my.clusters2 <- Idents(SLO_all.combined) 
cluster.select <- c("MRC","Adventitial","BEC","FDC_TRC","B","LEC","Myeloid","TRC","T","Epithelial","Periventricular","EPSTI1")
SLO_all.combined<-subset(SLO_all.combined,idents=cluster.select)
SLO_all.combined@meta.data$my.clusters2  <- factor(SLO_all.combined@meta.data$my.clusters2,levels = cluster.select)
DimPlot(SLO_all.combined, reduction = "umap",pt.size=0.1,label=T)+NoLegend()
ggsave2("umap.label2.png",width=5, height=5,device="png")
dim(SLO_all.combined)
save.image("temp.analyzed2.RData")


##############################

###Subset 
k="new1"
# Idents(obj)<-"seurat_clusters"
obj<-SLO_all.combined
obj<-subset(obj,idents=c("MRC","Adventitial","FDC_TRC","TRC","EPSTI1"))
obj<-subset(obj,PTPRC <1 & (PDPN>0 | CXCL12>0 |CCL19>0|CCL21>0|CR2>0))
# obj<-subset(obj,idents=c(0:1))
# obj<-subset(obj,idents=c(0,2:8))

#Harmony integration
obj <- obj %>% Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(vars.to.regress=c("nCount_RNA", "percent.mt","HSP.score1","percent.Rpl","percent.Rps"),verbose = FALSE) %>% 
  RunPCA(pc.genes = obj@var.genes, npcs = 20, verbose = FALSE)
obj<-RunHarmony(obj,c("GSE","Organ2","Species","Sample","Chemistry"),plot_convergence=T)
ElbowPlot(obj)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:10) #change based on elbow plot
obj <- RunTSNE(obj, reduction = "harmony", dims = 1:10) #change based on elbow plot
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:10) #change based on elbow plot
obj <- FindClusters(obj, resolution = 0.2) #adjust resolution (bigger=more clusters), initially used 0.2
# DimPlot(obj, reduction = "tsne")
obj2<-obj
DimPlot(obj2, reduction = "umap")
FeaturePlot(object = obj2, features = "CD19", reduction = "umap",pt.size=0.001, order=T)
FeaturePlot(object = obj2, features = c("CR2","CD274"), reduction = "umap",pt.size=0.001, order=T)
VlnPlot(obj2, features = "CD274",pt.size = 0.001, combine=F,raster=F)
for (i in c("Species","Organ2","Tx2","seurat_clusters","GSE","Sample","Chemistry")){
  Idents(obj)<-i
  DimPlot(obj, reduction = "umap",pt.size=0.1)
  ggsave2(paste0(i,".umap_subset",k,".png"),width=6, height=5,device="png")
}
for (i in c("nCount_RNA","nFeature_RNA","percent.mt","percent.Rps","percent.Rpl","CC.Difference","HSP.score1","S.Score","G2M.Score")){
  FeaturePlot(obj, features= i,pt.size=0.1, order=T)
  ggsave2(paste0(i,".umap_subset",k,".png"),width=5, height=5,device="png")
}
save.image(paste0("merged_subset",k,".RData"))

## Cluster Analysis
Idents(obj)<-"seurat_clusters"
DimPlot(obj, reduction = "umap",group.by="CC.label",pt.size=0.1,
        cells=obj@meta.data$cellID[!is.na(obj@meta.data$CC.label)])
ggsave2(paste0("cyster.umap_subset",k,".png"),width=6, height=5,device="png")
DimPlot(obj, reduction = "umap",pt.size=0.1,label=T)+NoLegend()
ggsave2(paste0("umap.label_subset",k,".png"),width=6, height=5,device="png")
gene.list<-c("MADCAM1","TNFSF11","CXCL13","CR2","CXCL12","CCL19","CCL21","PDPN","CXCL9","CXCL10",
             "CH25H","IL7","LEPR","NR4A1","LEPR","CD34","TNFSF13B","INMT","CD274")
FeaturePlot(object = obj, features = gene.list, 
            cols = c("grey", "blue"), reduction = "umap",pt.size=0.001, order=T)
ggsave2(paste0("umap.goi_subset",k,".png"),width=10, height=10,device="png")
VlnPlot(obj, features = "CD274",pt.size = 0.001, combine=F,raster=F)
ggsave2(paste0("vln.cd274_subset",k,".png"),width=4, height=3,device="png")


#cluster DE
Idents(obj)<-"seurat_clusters"
# obj.markers <- FindAllMarkers(object = obj, test.use = "MAST")
obj.markers <- FindAllMarkers(object = obj, test.use = "MAST",max.cells.per.ident=1000)
write.csv(obj.markers,paste0("cluster.markers_subset",k,".csv"))
write.csv(obj.markers %>% group_by(cluster) %>% top_n(20, avg_log2FC),paste0("cluster.markers_subset",k,".top20.csv"))
top6 <- obj.markers %>% group_by(cluster) %>% top_n(6, avg_log2FC)
DoHeatmap(object = obj, features = top6$gene, label = TRUE)  #slim.col.label to TRUE prints cluster IDs instead of cells
ggsave2(paste0("cluster.heatmap_subset",k,".png"),width=10, height=10,device="png")
obj@meta.data$my.clusters <- Idents(obj)  # Store cluster identities in object@meta.data$my.clusters
top4 <- obj.markers[obj.markers$cluster %in% c(0:5),] %>% group_by(cluster) %>% top_n(4, avg_log2FC)
FeaturePlot(object = obj2, features = unique(top4$gene),pt.size=0.001,order=T,
            cols = c("grey", "blue"), reduction = "umap")
ggsave2(paste0("umap.clustermarkers0_",k,".png"),width=20,height=20,device="png")
rm(list=setdiff(ls(), c("obj","obj.markers","k")))
save.image(paste0("temp.analyzed1_subset",k,".RData"))

#rename clusters
obj<-subset(obj,idents=c(0:5))
new.cluster.ids <- c("MRC","TRC","INMT_TRC","CD34_SC","FDC","MRC_ADAMDEC1")
names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj, new.cluster.ids)
obj@meta.data$my.clusters2 <- Idents(obj) 
DimPlot(obj, reduction = "umap",pt.size=0.1)
ggsave2(paste0("umap.label2_",k,".png"),width=6, height=5,device="png")
FeatureScatter(obj,"CXCL13","CD274")
ggsave2(paste0("cd274.scatter_",k,".png"),width=3, height=3,device="png")
dim(obj)
save.image(paste0("temp.analyzed2_subset",k,".RData"))

#########
##Condition DE
obj<-SLO_all.combined
Idents(obj) <- "my.clusters2" #setting idents to condition metadata
fdc_mrc.DE <- FindMarkers(obj, ident.1 = "FDC_TRC", ident.2 = "MRC", min.pct=0,logfc.threshold = -Inf,test.use = "MAST",max.cells.per.ident=10000)
write.csv(fdc_mrc.DE ,"FDC_mrc_DE.csv")

fdc_epi.DE <- FindMarkers(obj, ident.1 = "FDC_TRC", ident.2 = "Epithelial", min.pct=0,logfc.threshold = -Inf,test.use = "MAST",max.cells.per.ident=10000)
write.csv(fdc_epi.DE ,"FDC_epi_DE.csv")

fdc_trc.DE <- FindMarkers(obj, ident.1 = "FDC", ident.2 = "TRC", min.pct=0,logfc.threshold = -Inf,test.use = "MAST",max.cells.per.ident=10000)
write.csv(fdc_trc.DE ,"FDC_trc_DE.csv")

obj@meta.data$my.clusters3<-"cells"
obj@meta.data$my.clusters3[obj@meta.data$seurat_clusters==4]<-"fdc"
Idents(obj) <- "my.clusters3" #setting idents to condition metadata
fdc_cells.DE <- FindMarkers(obj, ident.1 = "fdc", ident.2 = "cells", min.pct=0,logfc.threshold = -Inf,test.use = "MAST",max.cells.per.ident=10000)
write.csv(fdc_cells.DE ,"FDC_cells.csv")

test<-FetchData(obj,vars=c("CD274"))
obj@meta.data$my.clusters4<-"pdl1_lo"
obj@meta.data$my.clusters4[test$CD274==0]<-"pdl1_neg"
obj@meta.data$my.clusters4[test$CD274>0.5]<-"pdl1_hi"
Idents(obj) <- "my.clusters4" #setting idents to condition metadata
pdl1.DE <- FindMarkers(obj, ident.1 = "pdl1_hi", ident.2 = "pdl1_neg", min.pct=0,logfc.threshold = -Inf,test.use = "MAST",max.cells.per.ident=10000)
write.csv(pdl1.DE ,"pdl1_vs.csv")

DimPlot(obj, cols.highlight=c("#D95F02","darkgreen"),
  cells.highlight = list("PD-L1 high"=Cells(subset(obj, idents="pdl1_hi")),"PD-L1 neg"=Cells(subset(obj, idents="pdl1_neg"))),
  sizes.highlight=c(1,0.2))+NoLegend()+ 
  NoAxes()+theme(panel.border = element_blank())
ggsave2("pdl1_hi_umap.png",width=5, height=5,device="png")

brewer.pal(n = 2, name = "Dark2")

Idents(obj)<-"Species"
hum<-subset(obj,idents="Hs")
hum@meta.data$my.clusters3<-"cells"
hum@meta.data$my.clusters3[hum@meta.data$seurat_clusters==4]<-"fdc"
Idents(hum) <- "my.clusters3" #setting idents to condition metadata
hum_fdc_cells.DE <- FindMarkers(hum, ident.1 = "fdc", ident.2 = "cells", min.pct=0,logfc.threshold = -Inf,test.use = "MAST",max.cells.per.ident=10000)
write.csv(hum_fdc_cells.DE ,"hum_FDC_cells.csv")
mus<-subset(obj,idents="Mus")
mus@meta.data$my.clusters3<-"cells"
mus@meta.data$my.clusters3[mus@meta.data$seurat_clusters==4]<-"fdc"
Idents(mus) <- "my.clusters3" #setting idents to condition metadata
mus_fdc_cells.DE <- FindMarkers(mus, ident.1 = "fdc", ident.2 = "cells", min.pct=0,logfc.threshold = -Inf,test.use = "MAST",max.cells.per.ident=10000)
write.csv(mus_fdc_cells.DE ,"mus_FDC_cells.csv")

save.image("analyzed.RData")
