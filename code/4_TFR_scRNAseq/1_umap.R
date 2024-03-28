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

#Harmony integration
T_all <- T_all %>% Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(vars.to.regress=c("nCount_RNA", "percent.mt","HSP.score1","percent.Rpl","percent.Rps"),verbose = FALSE) %>% 
  RunPCA(pc.genes = T_all@var.genes, npcs = 20, verbose = FALSE)
T_all.combined<-RunHarmony(T_all,c("condition","batch","mouse_ID"),plot_convergence=T)
T_all.combined <- RunUMAP(T_all.combined, reduction = "harmony", dims = 1:20) #change based on elbow plot
T_all.combined <- FindNeighbors(T_all.combined, reduction = "harmony", dims = 1:20) #change based on elbow plot
T_all.combined <- FindClusters(T_all.combined, resolution = 0.3) #adjust resolution (bigger=more clusters), initially used 0.2
DimPlot(T_all.combined, reduction = "umap")
for (i in c("condition","Phase","mouse_ID","gender","batch","seurat_clusters")){
  Idents(T_all.combined)<-i
  DimPlot(T_all.combined, reduction = "umap",pt.size=0.1)
  ggsave2(paste0(i,".umap.png"),width=6, height=5,device="png")
}
for (i in c("nCount_RNA","nFeature_RNA","percent.mt","percent.Rps","percent.Rpl","CC.Difference","HSP.score1","S.Score","G2M.Score")){
  FeaturePlot(T_all.combined, features= i,split.by = "condition",pt.size=0.1, order=T)
  ggsave2(paste0(i,".umap.png"),width=10, height=5,device="png")
}
rm(list=setdiff(ls(), c("T_all.combined", "metadata","genes.meta")))
save.image("merged.RData")

T_all.markers <- FindAllMarkers(object = T_all.combined, test.use = "MAST")
write.csv(T_all.markers,"cluster.markers.csv")
write.csv(T_all.markers %>% group_by(cluster) %>% top_n(20, avg_log2FC),"cluster.markers.top20.csv")
top6 <- T_all.markers %>% group_by(cluster) %>% top_n(6, avg_log2FC)
top6 <- T_all.markers %>% group_by(cluster) %>% top_n(4, avg_log2FC)
DoHeatmap(object = T_all.combined, features = top6$gene, label = TRUE)  #slim.col.label to TRUE prints cluster IDs instead of cells
ggsave2("cluster.heatmap.png",width=10, height=10,device="png")
T_all.combined@meta.data$my.clusters <- Idents(T_all.combined)  # Store cluster identities in object@meta.data$my.clusters

#rename clusters
new.cluster.ids <- c("TFR","Sostdc1","TFH-Tcf1","TFH-Exhausted","TFH-Activated","TFH-CM","TFH-Effector","TFH-ISG","NA")
names(new.cluster.ids) <- levels(T_all.combined)
T_all.combined <- RenameIdents(T_all.combined, new.cluster.ids)
T_all.combined@meta.data$my.clusters2 <- Idents(T_all.combined)  
T_all.combined<-subset(T_all.combined,idents=c("TFR","Sostdc1","TFH-Tcf1","TFH-Exhausted","TFH-Activated","TFH-CM","TFH-Effector","TFH-ISG"))
T_all.combined@meta.data$my.clusters2  <- factor(T_all.combined@meta.data$my.clusters2, levels = c("TFR","Sostdc1","TFH-Tcf1","TFH-Exhausted","TFH-Activated","TFH-CM","TFH-Effector","TFH-ISG"))
DimPlot(T_all.combined, reduction = "umap",pt.size=0.1,label=T)+NoLegend()
ggsave2("umap.label2.png",width=5.5, height=5,device="png")
dim(T_all.combined)
save.image("merged2.RData")

###Subset 
k="TFR"
obj<-T_all.combined
obj<-subset(obj,idents=c("TFR"))

#Harmony integration
obj <- obj %>% Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(vars.to.regress=c("nCount_RNA", "percent.mt","HSP.score1","percent.Rpl","percent.Rps"),verbose = FALSE) %>% 
  RunPCA(pc.genes = obj@var.genes, npcs = 20, verbose = FALSE)
obj<-RunHarmony(obj,c("batch","mouse_ID"),plot_convergence=T)
ElbowPlot(obj)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:15) #change based on elbow plot
obj <- RunTSNE(obj, reduction = "harmony", dims = 1:15) #change based on elbow plot
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:15) #change based on elbow plot
obj <- FindClusters(obj, resolution = 0.2) #adjust resolution (bigger=more clusters), initially used 0.2
# DimPlot(obj, reduction = "tsne")
obj2<-obj
DimPlot(obj2, reduction = "umap")
FeaturePlot(object = obj2, features = "Pdcd1", reduction = "umap",pt.size=0.001, order=T)
FeaturePlot(object = obj2, features = c("CR2","CD274"), reduction = "umap",pt.size=0.001, order=T)
VlnPlot(obj2, features = "Pdcd1",pt.size = 0.001, combine=F,raster=F)
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

#cluster DE
Idents(obj)<-"seurat_clusters"
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
new.cluster.ids <- c("TFR_1","TFR_2","TFR_3")
names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj, new.cluster.ids)
obj@meta.data$my.clusters2 <- Idents(obj) 
DimPlot(obj, reduction = "umap",pt.size=0.1)
ggsave2(paste0("umap.label2_",k,".png"),width=6, height=5,device="png")
dim(obj)


##Condition DE
test<-FetchData(obj,vars=c("Pdcd1"))
hist(test$Pdcd1)
obj@meta.data$my.clusters4<-"pdl1_mid"
obj@meta.data$my.clusters4[test$Pdcd1<1]<-"pd1_neg"
obj@meta.data$my.clusters4[test$Pdcd1>3]<-"pd1_hi"
Idents(obj) <- "my.clusters4" #setting idents to condition metadata
tfr_pd1.DE <- FindMarkers(obj, ident.1 = "pd1_hi", ident.2 = "pd1_neg", min.pct=0,logfc.threshold = -Inf,test.use = "MAST")
write.csv(tfr_pd1.DE,"tfr_pd1.csv")


FeaturePlot(obj, "Pdcd1") & scale_color_gradientn(colors = plasma(n = 10, direction = -1))
ggsave2("tfr_pd1_umap.png",width=4, height=4,device="png")

Idents(obj)<-"condition"
m564<-subset(obj,idents="m564")
m564@meta.data$my.clusters3<-"cells"
m564@meta.data$my.clusters3[m564@meta.data$seurat_clusters==4]<-"fdc"
Idents(m564) <- "my.clusters4" #setting idents to condition metadata
m564_pd1.DE <- FindMarkers(m564, ident.1 = "pd1_hi", ident.2 = "pd1_neg", min.pct=0,logfc.threshold = -Inf,test.use = "MAST")
write.csv(m564_pd1.DE ,"m564_FDC_cells.csv")
B6<-subset(obj,idents="AID")
B6@meta.data$my.clusters3<-"cells"
B6@meta.data$my.clusters3[B6@meta.data$seurat_clusters==4]<-"fdc"
Idents(B6) <- "my.clusters4" #setting idents to condition metadata
B6_pd1.DE <- FindMarkers(B6, ident.1 = "pd1_hi", ident.2 = "pd1_neg", min.pct=0,logfc.threshold = -Inf,test.use = "MAST")
write.csv(B6_pd1.DE ,"B6_FDC_cells.csv")


save.image("analyzed.RData")
