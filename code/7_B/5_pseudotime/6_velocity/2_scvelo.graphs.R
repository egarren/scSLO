rm(list=ls())
library(plyr)
library(Seurat)
library(dplyr)
library(Matrix)
library(knitr)
library(viridis)
library(ggplot2)
library(cowplot)
library(ggpubr)


scanpy.dpt <- read.table("sc.velo.txt", sep = "\t")
scanpy.dpt2 <- scanpy.dpt[2:nrow(scanpy.dpt),2:4, drop = FALSE]
colnames(scanpy.dpt2) <- c("velocity","velocity.confidence","velocity.pseudotime")
scanpy.dpt2<-as.data.frame(sapply(scanpy.dpt2,as.numeric))
rownames(scanpy.dpt2) <- scanpy.dpt[2:nrow(scanpy.dpt),1]
obj <- AddMetaData(obj, metadata = scanpy.dpt2)
obj@meta.data$scanpy.pseudo.rank <- rank(obj@meta.data$velocity.pseudotime) 

#plot
Idents(obj) <- "my.clusters2"
for(i in c("velocity","velocity.confidence","velocity.pseudotime")){
  FeaturePlot(obj, features= i, cols= viridis(100, begin = 0))+#labs(color="Pseudotime")+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),plot.title = element_blank(),
          axis.ticks=element_blank(),axis.line=element_blank(),axis.text=element_blank())
  ggsave2(paste0("./velocity3/umap.",i,".png"),width=5.5, height=4,device="png")
}
for(i in c("velocity","velocity.confidence","velocity.pseudotime","scanpy.pseudo.rank")){
  VlnPlot(obj, features = i,pt.size=0.0001)+ #NoLegend()+#labs(y="Pseudotime")+
    theme(axis.title.x=element_blank(),plot.title=element_blank())+labs(y=i)+
    stat_compare_means(label="p.signif",method = "t.test")
  ggsave2(paste0("./velocity3/vln.",i,".png"),width=6, height=4,device="png")
}

save.image("scv.graphed.RData")
