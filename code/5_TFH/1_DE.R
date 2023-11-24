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


###Subset 
k="TFH"
obj<-T_all.combined
obj<-subset(obj,idents=c("Sostdc1","TFH-CM","TFH-Activated","TFH-Tcf1","TFH-ISG","TFH-Effector","TFH-Exhausted"))


#Harmony integration
obj <- SCTransform(obj, method = "glmGamPoi", vars.to.regress = c("nCount_RNA", "percent.mt","HSP.score1","percent.Rpl","percent.Rps"), verbose = FALSE)
obj<-RunPCA(obj)
obj<-RunHarmony(obj,c("batch","mouse_ID"),plot_convergence=T)
ElbowPlot(obj)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:15) #change based on elbow plot
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:15) #change based on elbow plot
obj <- FindClusters(obj, resolution = 0.2) #adjust resolution (bigger=more clusters), initially used 0.2
DimPlot(obj, reduction = "umap")
DimPlot(obj, reduction = "umap",group.by="my.clusters3")
FeaturePlot(object = obj, features = "Pdcd1", reduction = "umap",pt.size=0.001, order=T)
VlnPlot(obj, features = "Pdcd1",pt.size = 0.001, combine=F,raster=F)
for (i in c("condition","mouse_ID","seurat_clusters","GSE","Sample","Chemistry","my.clusters2","my.clusters3")){
  Idents(obj)<-i
  DimPlot(obj, reduction = "umap",pt.size=0.1)
  ggsave2(paste0(i,".umap_subset",k,".png"),width=6, height=5,device="png")
}
for (i in c("nCount_RNA","nFeature_RNA","percent.mt","percent.Rps","percent.Rpl","CC.Difference","HSP.score1","S.Score","G2M.Score")){
  FeaturePlot(obj, features= i,pt.size=0.1, order=T)
  ggsave2(paste0(i,".umap_subset",k,".png"),width=5, height=5,device="png")
}
rm(list=setdiff(ls(), c("obj", "metadata","k")))
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
FeaturePlot(object = obj, features = unique(top4$gene),pt.size=0.001,order=T,
            cols = c("grey", "blue"), reduction = "umap")
ggsave2(paste0("umap.clustermarkers0_",k,".png"),width=20,height=20,device="png")
rm(list=setdiff(ls(), c("obj","obj.markers","k")))
save.image(paste0("temp.analyzed1_subset",k,".RData"))

#rename clusters
new.cluster.ids <- c("TFH_Tcf1","TFH_Activated","TFH_Exhausted","Sostdc1_1","TFH_Il31ra","TFH_CM","TFH_Effector","Sostdc1_2","TFH_ISG")
names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj, new.cluster.ids)
obj@meta.data$my.clusters3 <- Idents(obj) 
DimPlot(obj, reduction = "umap",pt.size=0.1,label=T)+
  NoLegend()+ NoAxes()+theme(panel.border = element_blank())
ggsave2(paste0("umap.label2_",k,".png"),width=6, height=5,device="png")
dim(obj)

obj2<-obj
DefaultAssay(obj2)<-"RNA"
obj2[["SCT"]]<-NULL
obj.sce <- as.SingleCellExperiment(obj2)
save(obj.sce,file="sce.filtered.RData")
show_col(hue_pal()(9))
hue_pal()(9)


obj2@meta.data$cellID2<-paste0(obj2@meta.data$mouse_ID,"_",gsub("^.*\\_","",obj2@meta.data$cellID))
rownames(obj2@meta.data)<-obj2@meta.data$cellID2
DefaultAssay(obj2)<-"RNA"
obj.sce <- as.SingleCellExperiment(obj2)
save(obj.sce,file="sce.monocle.RData")


##Condition DE
DefaultAssay(obj)<-"RNA"
test<-FetchData(obj,vars=c("Pdcd1"))
hist(test$Pdcd1)
obj@meta.data$my.clusters4<-"pdl1_mid"
obj@meta.data$my.clusters4[test$Pdcd1<1]<-"pd1_neg"
obj@meta.data$my.clusters4[test$Pdcd1>3]<-"pd1_hi"
Idents(obj) <- "my.clusters4" #setting idents to condition metadata
tfh_pd1.DE <- FindMarkers(obj, ident.1 = "pd1_hi", ident.2 = "pd1_neg", min.pct=0,logfc.threshold = -Inf,test.use = "MAST")
write.csv(tfh_pd1.DE,"tfh_pd1.csv")


FeaturePlot(obj, "Pdcd1") & scale_color_gradientn(colors = plasma(n = 10, direction = -1))
ggsave2("tfh_pd1_umap.png",width=4, height=4,device="png")

Idents(obj)<-"condition"
m564<-subset(obj,idents="m564")
Idents(m564) <- "my.clusters4" #setting idents to condition metadata
m564_pd1.DE <- FindMarkers(m564, ident.1 = "pd1_hi", ident.2 = "pd1_neg", min.pct=0,logfc.threshold = -Inf,test.use = "MAST")
write.csv(m564_pd1.DE ,"m564_pd1_DE.csv")
B6<-subset(obj,idents="AID")
Idents(B6) <- "my.clusters4" #setting idents to condition metadata
B6_pd1.DE <- FindMarkers(B6, ident.1 = "pd1_hi", ident.2 = "pd1_neg", min.pct=0,logfc.threshold = -Inf,test.use = "MAST")
write.csv(B6_pd1.DE ,"B6_pd1_DE.csv")


save.image("analyzed.RData")


###Gene correlations
goi<-c("Selplg","Vsir","Cxcr5","Foxp3","Ly6a","Pdcd1","Il7r","Cxcr3","Ccr6","Sell","Il2ra",
       "Itgb1","Cd40lg","Id3","Lag3","Sostdc1")
for(i in goi){
  for(k in c("B6","m564","obj")){
    obj2<-get(k)
    df<-as.data.frame(t(as.matrix(obj2[["RNA"]]@data)))
    mod_score<-as.numeric(df[[i]])
    mat<-cbind(mod_score,df)
    cor.scores<-data.frame(gene=rep(NA,ncol(mat)),pval=rep(NA,ncol(mat)),rho=rep(NA,ncol(mat)))
    for(j in 1:ncol(mat)){
      cor<-cor.test(mat[[j]],mat$mod_score,method="spearman")
      cor.scores$pval[j]<-cor$p.value
      cor.scores$rho[j]<-cor$estimate
      cor.scores$gene[j]<-colnames(mat)[j]
    }
    cor.scores[is.na(cor.scores)]<-0
    assign(paste0(i,".",k,".gene.correl"),cor.scores[cor.scores$gene!="mod_score"&cor.scores$gene!=i,])
  }
  
  #by cluster
  for(l in unique(obj@meta.data$my.clusters3)){
    Idents(obj)<-"my.clusters3"
    clust<-subset(obj,idents=l)
    for(k in c("AID","m564","ALL")){
      Idents(clust)<-"condition"
      if(k %in% c("AID","m564")){obj2<-subset(clust,idents=k)}else{obj2<-clust}
      df<-as.data.frame(t(as.matrix(obj2[["RNA"]]@data)))
      mod_score<-as.numeric(df[[i]])
      mat<-cbind(mod_score,df)
      cor.scores<-data.frame(gene=rep(NA,ncol(mat)),pval=rep(NA,ncol(mat)),rho=rep(NA,ncol(mat)))
      for(j in 1:ncol(mat)){
        cor<-cor.test(mat[[j]],mat$mod_score,method="spearman")
        cor.scores$pval[j]<-cor$p.value
        cor.scores$rho[j]<-cor$estimate
        cor.scores$gene[j]<-colnames(mat)[j]
      }
      cor.scores[is.na(cor.scores)]<-0
      assign(paste0(i,".clust",l,".",k,".gene.correl"),cor.scores[cor.scores$gene!="mod_score"&cor.scores$gene!=i,])
      colnames(cor.scores)<-c("gene",paste0(k,".pval"),paste0(k,".rho"))
    }
    
  }
  save.image(paste0(i,".temp.gene.correl.RData"))
}


save.image("analyzed2.RData")

dir.create("./gene.correl")
setwd("./gene.correl")
#correl compare
for(i in ls(pattern=".gene.correl")){
  x<-get(i)
  x[is.na(x)]<-0
  for(j in ls(pattern=".gene.correl")){
    if(i !=j){
      for(k in goi){ 
        if(grepl(k,i) & grepl(k,j)){
          y<-get(j)
          y[is.na(y)]<-0
          df<-left_join(x[,c("rho","gene")],y[,c("rho","gene")],by="gene",keep=F)
          colnames(df)<-c("x","gene","y")
          rownames(df)<-df$gene
          df$diff= abs(df$x-df$y)
          df$sum = df$x + df$y
          df<-df %>% mutate(rank_diff = rank(-diff), rank_sum = rank(sum), rank_sum_dec = rank(-sum)) 
          gene.list<-df$gene[df$rank_diff<10 | df$rank_sum <5 |df$rank_sum_dec<10]
          df$label<-"NS"
          df$label[df$x>0.3 & df$y>0.3]<-"up"
          df$label[df$x< -0.3 & df$y< -0.3]<-"down"
          df$label[df$diff>1]<-"diff"
          df<-df[order(-df$rank_diff),]
          col_up=viridis(3)[1]
          col_down=viridis(3)[2]
          p2 <- ggplot(df, aes(x,y)) + geom_point(aes(colour=label),size=0.5)+#geom_point(fill=NA,colour=alpha("black",0.5),pch=21,size=3) +
            theme_classic()+
            scale_color_manual(values = c("up" = col_up, "down" = col_down, "NS" = "grey80","diff"="black"))+
            labs(x=i,y=j)+
            theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position="none",
                  axis.ticks=element_blank(),axis.line=element_blank())+ #axis.text=element_blank(),
            geom_hline(yintercept=0,color=alpha("black",0.5))+geom_vline(xintercept=0,color=alpha("black",0.5))
          if(length(gene.list)>0){p2 <- LabelPoints(plot = p2, points = gene.list, repel = TRUE,xnudge=0,ynudge=0,size=2.5,segment.size=0.1)}
          p2
          ggsave2(paste0(i,"_",j,"_gene_compare.png"),width=2.5, height=2.5,device="png")
        }
      }
    }
  }
}

load("analyzed2.RData")
goi<-c("Selplg","Vsir","Cxcr5","Foxp3","Ly6a","Il7r","Cxcr3","Ccr6","Sell","Il2ra",
       "Itgb1","Cd40lg","Id3","Lag3","Sostdc1")
#Scatter graphs
Idents(obj)<-"my.clusters3"
for(h in goi){ 
  FeatureScatter(object = obj, feature1 = "Pdcd1", feature2 = h,pt.size=0.1)+
    xlab("Pdcd1")+ylab(h)+ NoLegend()+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.ticks=element_blank(),axis.line=element_blank())
  ggsave2(paste0(h,".pdcd1.correl.scatter.png"),width=2.5, height=2.7,device="png")
}
