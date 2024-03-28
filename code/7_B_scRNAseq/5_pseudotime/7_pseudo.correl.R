rm(list=ls())
library(org.Mm.eg.db)
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
# library(VISION)
library(msigdbr)
library(KEGGREST)
# library(msigdb)
# library(goseq)
# library(nicethings)
library(DOSE)
library(ggpubr)
library(cogena)
library(viridisLite)
library(viridis)
library(gridExtra)

temp.obj<-obj

load("../scanpy.pseudotime.RData")
temp.obj@meta.data$scanpy.pseudotime<-obj@meta.data$scanpy.pseudotime
temp.obj@meta.data$scanpy.pseudo.rank <- rank(obj@meta.data$scanpy.pseudotime)
load("../destiny.scanpy.pseudotime.RData")
temp.obj@meta.data$destiny.pseudotime<-obj@meta.data$destiny.pseudotime
temp.obj@meta.data$destiny.pseudo.rank <- rank(obj@meta.data$destiny.pseudotime)
load("../temp.slingshot.DE2.RData")
temp.obj@meta.data$slingshot.pseudotime<-obj@meta.data$slingshot.pseudotime
temp.obj@meta.data$slingshot.pseudo.rank <- rank(obj@meta.data$slingshot.pseudotime)
load("../temp2.monocle.RData")
temp.obj@meta.data$monocle.pseudotime<-cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
temp.obj@meta.data$monocle.pseudo.rank <- rank(temp.obj@meta.data$monocle.pseudotime)
load("../scv.graphed.RData")
temp.obj@meta.data$velocity.pseudotime<-obj@meta.data$velocity.pseudotime
temp.obj@meta.data$velocity.pseudo.rank <- rank(temp.obj@meta.data$monocle.pseudotime)
obj<-temp.obj
save.image("temp0.pseudo.correl.RData")
save(obj,file="temp.pseudo.RData")

#Plots
Idents(obj)<-"my.clusters2"
for(h in c("scanpy","destiny","slingshot","monocle","velocity")){ #,
  FeaturePlot(obj, features= paste0(h,".pseudotime"), cols= viridis(100, begin = 0))+labs(color="Pseudotime")+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),plot.title = element_blank(),
          axis.ticks=element_blank(),axis.line=element_blank(),axis.text=element_blank())
  ggsave2(paste0(h,".umap.pseudo.png"),width=5.5, height=4,device="png")
  if(h != "slingshot"){
    VlnPlot(obj, features = paste0(h,".pseudotime"),pt.size=0)+ NoLegend()+labs(y="Pseudotime")+
      theme(axis.title.x=element_blank(),plot.title=element_blank())
    ggsave2(paste0(h,".vln.pseudo.png"),width=4, height=4,device="png")
    VlnPlot(obj, features = paste0(h,".pseudotime"), split.by = "Species",
            group.by = "my.clusters2",cols=c("red","grey"),pt.size = 0,combine=F)
    ggsave2(paste0(h,".vln2.pseudo.png"),width=5, height=4,device="png")
  }
  VlnPlot(obj, features = paste0(h,".pseudo.rank"),pt.size=0)+ NoLegend()+labs(y="Pseudotime Rank")+
    theme(axis.title.x=element_blank(),plot.title=element_blank())
  ggsave2(paste0(h,".vln.pseudo.rank.png"),width=4, height=4,device="png")
  if(h != "velocity"){
    FeaturePlot(obj, features =paste0(h,".pseudotime"), split.by = "Species",
                cols=viridis(100, begin = 0),order=T,pt.size=0.001,combine=T)
    ggsave2(paste0(h,".umap2.pseudo.png"),width=6.5, height=3,device="png")
    p<-FeaturePlot(obj, features = paste0(h,".pseudotime"), split.by = "Species",cols=viridis(100, begin = 0),
                   order=T,pt.size=0.001,combine=F)
    for(j in 1:length(p)) {
      p[[j]] <- p[[j]] + NoLegend()+NoAxes()+
        theme(panel.border = element_rect(colour = "black", size=1),
              plot.title=element_blank(),axis.title.y.right=element_blank(),
              axis.line=element_blank())
    }
    cowplot::plot_grid(p[[1]],p[[2]],ncol=2)
    ggsave2(paste0(h,".umap2.pseudo.clean.png"),width=6, height=2.75,device="png")
    VlnPlot(obj, features = paste0(h,".pseudo.rank"), split.by = "Species",
            group.by = "my.clusters2",cols=c("red","grey"),pt.size = 0,combine=F)
    ggsave2(paste0(h,".vln2.pseudo.rank.png"),width=5, height=4,device="png")
  }
}

#Pseudotime correlations
df<-as.data.frame(t(as.matrix(obj[["RNA"]]@data)))
for(h in c("scanpy.pseudotime","slingshot.pseudotime","destiny.pseudotime","monocle.pseudotime","velocity.pseudotime")){ #
  pseudo_score<-as.numeric(obj@meta.data[[h]])
  mat<-cbind(pseudo_score,df)
  ggscatter(mat,x="Cd74",y="Foxp3",add="reg.line",conf.int=T,cor.coef=T,cor.method="spearman")
  correlations<-apply(matrix,1,function(x){cor(pseudo_score,x)})
  cors<-apply(mat,2,cor.test,pseudo_score,method="spearman") #non-parametric (pearson for parametric)
  cor.scores<-data.frame(gene=rep(NA,ncol(mat)),pval=rep(NA,ncol(mat)),rho=rep(NA,ncol(mat)))
  for(j in 1:ncol(mat)){
    cor<-cor.test(mat[[j]],mat$pseudo_score,method="spearman")
    cor.scores$pval[j]<-cor$p.value
    cor.scores$rho[j]<-cor$estimate
    cor.scores$gene[j]<-colnames(mat)[j]
  }
  cor.scores[is.na(cor.scores)]<-0
  assign(paste0(h,".pseudo.correl"),cor.scores[cor.scores$gene!="pseudo_score",])
  write.csv(file=paste0(h,".pseudo.correl.csv"),cor.scores[cor.scores$gene!="pseudo_score",])
  #by cluster
  for(l in unique(obj@meta.data$my.clusters2)){
    Idents(obj)<-"my.clusters2"
    obj2<-subset(obj,idents=l)
    df<-as.data.frame(t(as.matrix(obj2[["RNA"]]@data)))
    pseudo_score<-as.numeric(obj2@meta.data[[h]])
    if(sum(!is.na(pseudo_score))>20){
      mat<-cbind(pseudo_score,df)
      cor.scores<-data.frame(gene=rep(NA,ncol(mat)),pval=rep(NA,ncol(mat)),rho=rep(NA,ncol(mat)))
      for(j in 1:ncol(mat)){
        cor<-cor.test(mat[[j]],mat$pseudo_score,method="spearman")
        cor.scores$pval[j]<-cor$p.value
        cor.scores$rho[j]<-cor$estimate
        cor.scores$gene[j]<-colnames(mat)[j]
      }
      cor.scores[is.na(cor.scores)]<-0
      assign(paste0(h,".clust",l,".pseudo.correl"),cor.scores[cor.scores$gene!="pseudo_score",])
      write.csv(file=paste0(h,".clust",l,".pseudo.correl.csv"),cor.scores[cor.scores$gene!="pseudo_score",])
    }
  }
  #by condition
  for(l in unique(obj@meta.data$Tx)){
    Idents(obj)<-"Tx"
    obj2<-subset(obj,idents=l)
    df<-as.data.frame(t(as.matrix(obj2[["RNA"]]@data)))
    pseudo_score<-as.numeric(obj2@meta.data[[h]])
    if(sum(!is.na(pseudo_score))>20){
      mat<-cbind(pseudo_score,df)
      cor.scores<-data.frame(gene=rep(NA,ncol(mat)),pval=rep(NA,ncol(mat)),rho=rep(NA,ncol(mat)))
      for(j in 1:ncol(mat)){
        cor<-cor.test(mat[[j]],mat$pseudo_score,method="spearman")
        cor.scores$pval[j]<-cor$p.value
        cor.scores$rho[j]<-cor$estimate
        cor.scores$gene[j]<-colnames(mat)[j]
      }
      cor.scores[is.na(cor.scores)]<-0
      assign(paste0(h,".",l,".condition.pseudo.correl"),cor.scores[cor.scores$gene!="pseudo_score",])
      write.csv(file=paste0(h,".",l,".condition.pseudo.correl.csv"),cor.scores[cor.scores$gene!="pseudo_score",])
    }
  }
  save.image(paste0(h,".temp.pseudo.correl.RData"))
}
save.image("temp.pseudo.correl.RData")


#correl compare
for(i in ls(pattern=".pseudo.correl")){
  x<-get(i)
  x[is.na(x)]<-0
  for(j in ls(pattern=".pseudo.correl")){
    if(i !=j){
      for(k in c("scanpy","destiny","slingshot","monocle","velocity")){ 
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
          ggsave2(paste0(i,"_",j,"_pseudo_compare.png"),width=2.5, height=2.5,device="png")
        }
      }
    }
  }
}

#psuedotime comparisons
dfs<-paste0(c("scanpy.pseudotime","slingshot.pseudotime","destiny.pseudotime","monocle.pseudotime","velocity.pseudotime"),".pseudo.correl")
for(i in dfs){
  x<-get(i)
  x[is.na(x)]<-0
  for(j in dfs){
    if(i !=j){
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
      ggsave2(paste0(i,"_",j,"_pseudo_compare.png"),width=2.5, height=2.5,device="png")
    }
  }
}

#Scatter graphs
DefaultAssay(obj)<-"RNA"
Idents(obj)<-"my.clusters2"
for(h in c("scanpy","destiny","slingshot","monocle","velocity")){ 
  FeatureScatter(object = obj, feature1 = paste0(h,".pseudotime"), feature2 = "CD274")+
    scale_y_continuous(trans='log10')+
    xlab("Pseudotime")+ylab("Average expression")+ NoLegend()+ggtitle("CD274")+ 
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.ticks=element_blank(),axis.line=element_blank())
  ggsave2(paste0(h,".pseudo.correl.CD274.scatter.png"),width=5, height=5,device="png")
}

h="slingshot"
FeatureScatter(object = obj, feature1 = paste0(h,".pseudotime"), feature2 = "CD274")+
  scale_y_continuous(trans='log10')+
  geom_smooth(color="black")+
  xlab("Pseudotime")+ylab("Average expression")+ NoLegend()+ggtitle("CD274")+ 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.ticks=element_blank(),axis.line=element_blank(),title=element_blank(),
        plot.title=element_blank())
ggsave2(paste0(h,".pseudo.correl.CD274.scatter2.png"),width=3, height=3,device="png")

for(h in c("scanpy","destiny","slingshot","monocle","velocity")){ 
  p<-list()
  for(i in c("CD274","PDCD1","AICDA","S1PR2")){ #"Thy1","Cd40lg","Ikzf2","Foxp3","Il1r2"
    p[[i]]<-  FeatureScatter(object = obj, feature1 = paste0(h,".pseudotime"), feature2 = i)+
      scale_y_continuous(trans='log10')+ #xlab("Pseudotime")+ylab("Average expression")+ 
      NoLegend()+ggtitle(i)+ 
      theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
            axis.ticks=element_blank(),axis.line=element_blank(),axis.title=element_blank())
  }
  cowplot::plot_grid(plotlist=p,ncol=4)
  ggsave2(paste0(h,".pseudo.correl.scatter.png"),width=12, height=3,device="png")
}

df<-slingshot.pseudotime.pseudo.correl
df$rank<-rank(df$rho,ties.method="random")
df<-df[order(df$rank),]
head(df$gene,n=100)
tail(df$gene,n=100)
gene.list<-c("CD52","RPL32","FCMR","CD37","SELL","CD74","IGHM","MS4A1","JUND","FAU","CD274",
             "C14orf119","KDELR2","CD63","FABP5","VDAC1","ICAM3","PRDM1","JCHAIN","CD27")
df<- df %>% mutate(gene_label = ifelse(gene %in% gene.list , gene, NA))
df$trend<-"none"
df$trend[df$rho>0.1]<-"up"
df$trend[df$rho< -0.1]<-"down"
df$gene_label[df$trend !="none" & is.na(df$gene_label)]<-""
ggplot(df,aes(x=rank,y=rho,label=gene_label,color=trend))+geom_point()+
  scale_color_manual(values=c("red","black","blue"))+
  geom_text_repel(size = 3, segment.size = 0.2,seed = 1, show.legend=F, box.padding=0.3,max.overlaps=Inf) +
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position="none",
        axis.ticks=element_blank(),axis.line=element_blank(),title=element_blank())
ggsave2("pseudo.correl.ranks.png",width=3, height=3,device="png")

