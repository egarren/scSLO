
rm(list=ls())
library(Seurat)
library(cowplot)
library(destiny)
library(ggplot2)
library(ggbeeswarm)
library(ggthemes)
library(scater)
library(plyr)
library(dplyr)
library(Matrix)
library(knitr)
library(viridis)


#PCA
obj.sce <- runPCA(obj.sce, ncomponents = 50)
pca <- reducedDim(obj.sce, "PCA")
head(pca)
dim(pca)
destiny.PCA<-data.frame(PC1<-pca[, 1],PC2<-pca[, 2])
colnames(destiny.PCA) <- c("destiny.PC1", "destiny.PC2")
rownames(destiny.PCA)<-colnames(obj.sce)
destiny.PCA<-as.matrix(destiny.PCA)

rm(list=setdiff(ls(), c("destiny.PCA","obj.sce","pca")))
load("scanpy.pseudotime.RData")

#adding to seurat metadata
Idents(obj)<-"my.clusters2"
obj[["destiny.PCA"]] <- CreateDimReducObject(embeddings = destiny.PCA, key = "destinyPC_", assay = DefaultAssay(obj))
DimPlot(obj, reduction= "destiny.PCA")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.ticks=element_blank(),axis.line=element_blank(),axis.text=element_blank()) 
ggsave2("destiny.PC.png",width=5, height=4,device="png")
#rank by PC1
obj@meta.data$destinyPC1<-rank(destiny.PCA[,1]) 
VlnPlot(obj, features = "destinyPC1",pt.size=0.01)+ NoLegend()+labs(y="PC1 Rank")+
  theme(axis.title.x=element_blank(),plot.title=element_blank())
ggsave2("cluster.by.destinyPC1.png",width=6, height=4,device="png")

#Diffusion on PC
cellLabels <- obj.sce$seurat_clusters
rownames(pca) <- cellLabels
sigmas <- find_sigmas(pca, verbose = FALSE)  # find optimal sigma
dm.pc <- DiffusionMap(pca, sigma = optimal_sigma(sigmas))
destiny.diffmap<-data.frame(DC1<-eigenvectors(dm.pc)[, 1],DC2 <-eigenvectors(dm.pc)[, 2])
colnames(destiny.diffmap) <- c("destiny.DC1", "destiny.DC2")
rownames(destiny.diffmap)<-colnames(obj.sce)
destiny.diffmap$destiny.DC1<-destiny.diffmap$destiny.DC1*1000
destiny.diffmap$destiny.DC2<-destiny.diffmap$destiny.DC2*1000
destiny.diffmap<-as.matrix(destiny.diffmap)
obj[["destiny.diffmap"]] <- CreateDimReducObject(embeddings = destiny.diffmap, key = "destinyDC_", 
                                                            assay = DefaultAssay(obj))
#rank by DC1
obj@meta.data$destinyDC1<-rank(eigenvectors(dm.pc)[,1])  
VlnPlot(obj, features = "destinyDC1",pt.size=0.01)
ggsave2("cluster.by.destinyDC1.png",width=6, height=4,device="png")

#Diffusiontime
test<-obj@meta.data[obj@meta.data$my.clusters2=="Naive",]
which.min(test$destinyDC1)
dpt <- DPT(dm.pc, tips = 5258) #or dm.pc
obj@meta.data$destiny.pseudotime<-dpt$dpt
obj@meta.data$destiny.pseudo.rank <- rank(obj@meta.data$destiny.pseudotime) 

#comparison plots
Idents(obj) <- "my.clusters2"
DimPlot(obj, reduction= "destiny.diffmap")+ labs(x="DC_1",y="DC_2")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.ticks=element_blank(),axis.line=element_blank(),axis.text=element_blank()) #
ggsave2("destiny.diffmap.png",width=6, height=4,device="png")
FeaturePlot(obj, features= "destiny.pseudotime", cols= viridis(100, begin = 0), reduction = "destiny.diffmap")+ 
  labs(x="DC_1",y="DC_2",color="Pseudotime")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),plot.title = element_blank(),
        axis.ticks=element_blank(),axis.line=element_blank(),axis.text=element_blank())
ggsave2("destiny.diffmap.pseudo.png",width=5.5, height=4,device="png")
FeaturePlot(obj, features= "destiny.pseudotime", cols= viridis(100, begin = 0))+labs(color="Pseudotime")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),plot.title = element_blank(),
        axis.ticks=element_blank(),axis.line=element_blank(),axis.text=element_blank())
ggsave2("umap.destiny.pseudo.png",width=5.5, height=4,device="png")
VlnPlot(obj, features = "destinyDC_1",pt.size=0)+ NoLegend()+labs(y="DC_1")+
  theme(axis.title.x=element_blank(),plot.title=element_blank())
ggsave2("cluster.by.destinyDC1.png",width=4, height=4,device="png")
VlnPlot(obj, features = "destiny.pseudotime",pt.size=0)+ NoLegend()+labs(y="Pseudotime")+
  theme(axis.title.x=element_blank(),plot.title=element_blank())
ggsave2("cluster.by.destinypseudo.png",width=4, height=4,device="png")
VlnPlot(obj, features = "destiny.pseudo.rank",pt.size=0)+ NoLegend()+labs(y="Pseudotime Rank")+
  theme(axis.title.x=element_blank(),plot.title=element_blank())
ggsave2("cluster.by.destinypseudo.rank.png",width=4, height=4,device="png")
FeatureScatter(obj, feature1 = "destinyDC_1", feature2 = "PC_1")+labs(color="")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),plot.title = element_blank(),
        axis.ticks=element_blank(),axis.line=element_blank(),axis.text=element_blank())
ggsave2("destiny.DC.PC.png",width=6, height=4,device="png")

save.image("destiny.scanpy.pseudotime.RData")
