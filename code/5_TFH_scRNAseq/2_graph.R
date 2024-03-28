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
library(ggalluvial)
library(tidyverse)
library(gridExtra)
library(corrplot)
library(msigdbr)
my.ttest <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}


obj2<-obj
top8 <- obj.markers %>% group_by(cluster) %>% top_n(8, avg_log2FC)
top6 <- obj.markers %>% group_by(cluster) %>% top_n(6, avg_log2FC)
top4 <- obj.markers %>% group_by(cluster) %>% top_n(4, avg_log2FC)
top2 <- obj.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)
k="TFH"

# ## tSNE/UMAP
# # Visualization
Idents(obj2)<-"my.clusters3"
DimPlot(obj2, reduction = "umap")
DimPlot(obj2, reduction = "umap",pt.size=0.001)+ #NoLegend()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
                   axis.text=element_blank(),axis.ticks=element_blank(),axis.line=element_blank())
ggsave2(paste0("umap_",k,".png"),width=4.5, height=3.5,device="png")

DimPlot(obj2, reduction = "umap",pt.size=0.1,label=T)+
  NoLegend()+ NoAxes()+theme(panel.border = element_blank())
ggsave2(paste0("umap.label_subset",k,".png"),width=3, height=3,device="png")
ggsave2(paste0("umap.label_subset",k,".png"),width=4, height=4,device="png")

FeaturePlot(obj2, "Pdcd1") + #NoLegend()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text=element_blank(),axis.ticks=element_blank(),axis.line=element_blank())& 
  scale_color_gradientn(colors = plasma(n = 10, direction = -1))
ggsave2(paste0(k,"_pd1_umap.png"),width=4, height=3.7,device="png")
ggsave2(paste0(k,"_pd1_umap.png"),width=5, height=4.5,device="png")

FeaturePlot(obj, "Pdcd1",order=T)+ 
  NoAxes()+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
                 plot.title = element_blank()) & 
  scale_color_gradientn(colors = plasma(n = 10, direction = -1))
ggsave2("pd1_umap.png",width=4.5, height=4,device="png")


#UMAP by cluster marker
gene.list<-c("Tcf7","Sostdc1","Lmo4","Ext1","Klf2","Slamf7","Ccl5","Ifit3","Il31ra","Sell") #fdc_0
p<-FeaturePlot(object = obj2, features = gene.list,pt.size=0.001,order=T,
               cols = c("grey", "blue"), reduction = "umap",combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
}
cowplot::plot_grid(plotlist=p,ncol=4)
ggsave2(paste0("umap.clustermarkers_",k,".png"),width=8.5,height=7,device="png")

#vln cluster marker
VlnPlot(obj2, features = "Pdcd1",group.by = "my.clusters3",pt.size = 0.01)
ggsave2(paste0("vlnplot.pd1_",k,".png"),width=3, height=4,device="png")
VlnPlot(obj2, features = "Pdcd1",group.by = "my.clusters3",pt.size = 0)
ggsave2(paste0("vlnplot.pd1_",k,".png"),width=5, height=3,device="png")
VlnPlot(obj2, features = "Pdcd1",group.by = "my.clusters3",pt.size = 0)+
  stat_summary(fun.data=mean_cl_normal, mult=1,geom="pointrange", color="black")+ theme(legend.position = "none")
ggsave2(paste0("vlnplot.pd1_mean_",k,".png"),width=3, height=4,device="png")
VlnPlot(obj2, features = "Pdcd1",group.by = "my.clusters3",split.by="condition",pt.size = 0.01)
ggsave2(paste0("vlnplotsplit.pd1_",k,".png"),width=4, height=3,device="png")
p<-VlnPlot(obj2, features = gene.list,group.by = "my.clusters3",pt.size = 0.01, combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]]+ NoLegend()+theme(axis.title=element_blank())
}
cowplot::plot_grid(plotlist = p,ncol=5)
ggsave2(paste0("vlnplot.clustermarkers_",k,".png"),width=10, height=8,device="png")
RidgePlot(obj2, features = "Pdcd1",group.by = "my.clusters3")+ theme(legend.position = "none")
ggsave2(paste0("ridgeplot.pd1_",k,".png"),width=4, height=4,device="png")


#condition expression correlation
obj2@meta.data$group_by<-paste0(obj2@meta.data$condition,"-",obj2@meta.data$my.clusters3)
Idents(obj2)<-"group_by"
df<-AverageExpression(obj2)$RNA
res<-cor(as.matrix(df))
res<-as.data.frame(res)
res<-res[rownames(res) %like% "AID",colnames(res) %like% "m564"]
colnames(res)<-gsub("m564-","",colnames(res))
rownames(res)<-gsub("AID-","",rownames(res))
colnames(res)<-gsub("-","_",colnames(res))
rownames(res)<-gsub("-","_",rownames(res))
order<-names(table(obj2@meta.data$my.clusters3))
res<-res[order,order]
png(paste0(k,"_corrplot.png"), width = 5, height = 5, res = 200,units="in")
corrplot(as.matrix(res),
         is.corr = F,
         outline = T,
         col = rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#FFFFFF"))(100)),
         tl.col = "black",
         cl.pos="n",
         na.label = " ")
dev.off()

#dot plot
Idents(obj2)<-"my.clusters3"
top6$gene<-make.unique(top6$gene,sep="--")
DotPlot(obj2, features = top6$gene)+ 
  coord_flip()+theme(legend.title = element_blank(),axis.title=element_blank())+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))#+ RotatedAxis()
ggsave2(paste0("dotplot_",k,".png"),width=5, height=9,device="png")


#by condition
FeaturePlot(object = obj2, features = "Pdcd1", cols = c("grey", "blue"), 
            split.by="condition",pt.size=0.001, order=T)
ggsave2("umap.pd1.condition.png",width=6, height=3,device="png")

# Cluster frequency comparison
ggplot(obj2@meta.data, aes(x = my.clusters3,fill = condition)) +  
  geom_bar(aes(y = (..count..)/sum(..count..)),position='dodge')+ 
  scale_y_continuous(labels = percent)+ 
  labs(y="Frequency", x = "Cluster")
ggsave2("clusterfreq.png",width=4, height=2.5,device="png")

#bycondition
cluster.Counts <- table(obj2@meta.data$my.clusters3,obj2@meta.data$condition)
cluster.prop <- as.data.frame(scale(cluster.Counts,scale=colSums(cluster.Counts),center=FALSE)*100) 
ggplot(cluster.prop, aes(fill=Var1,y=Freq, x=Var2,alluvium=Var1,stratum=Var1)) + 
  geom_lode()+geom_flow()+geom_stratum(alpha=0) +theme_classic()+
  theme(legend.title = element_blank(),axis.title=element_blank()) +
  coord_flip()
ggsave2("cluster.prop.flow.condition.png",width=6, height=1.5,device="png")

###Volcano
for(j in ls(pattern=".DE")){
  res<-get(j)
  thresh_p_val_adj <- 1e-10
  thresh_lfc <-0.3
  plt_df<- res %>% rownames_to_column(var = "gene") %>% 
    mutate(up_in = ifelse(p_val_adj >= thresh_p_val_adj, "NS", ifelse(avg_log2FC > thresh_lfc, "up", ifelse(avg_log2FC < -thresh_lfc, "down", "NS"))),
           rank_pval = rank(p_val_adj), rank_lfc_inc = rank(avg_log2FC), rank_lfc_dec = rank(-abs(avg_log2FC)),
           lp = -log10(p_val_adj)) %>% 
    arrange(-abs(avg_log2FC))
  if (sum(plt_df$up_in != "NS") >= 30) {
    plt_df <- plt_df %>% mutate(gene_label = ifelse(up_in != "NS" & (rank_pval < 10 | rank_lfc_inc < 10 | rank_lfc_dec < 10), gene, NA))
  } else {
    plt_df <- plt_df %>% mutate(gene_label = ifelse((up_in != "NS" | rank_pval < 10 | rank_lfc_inc < 10 | rank_lfc_dec < 10), gene, NA))
  }
  plt_df$gene_label[plt_df$gene=="Pdcd1"]<-"Pdcd1"
  plt_df$gene_label[plt_df$gene=="Sostdc1"]<-"Sostdc1"
  plt_df$lp <- pmin(plt_df$lp, 300)
  plt_df$avg_log2FC <- pmin(plt_df$avg_log2FC, 3)
  ggplot(plt_df, aes(x = avg_log2FC, y = lp, color = up_in, label = gene_label)) +
    geom_point(size = 0.5) +
    geom_text_repel(color="black",size = 3, segment.size = 0.2,seed = 1, show.legend=F, box.padding=0.5) +
    geom_hline(yintercept=10,linetype=2,size=0.2)+
    geom_vline(xintercept=0.3,linetype=2,size=0.2)+
    geom_vline(xintercept= -0.3,linetype=2,size=0.2)+
    scale_color_manual(values = c("up" = viridis_pal(option="plasma")(8)[1], "down" = viridis_pal(option="plasma")(8)[6], "NS" = "grey80")) +
    scale_x_continuous(limits = c(-max(plt_df$avg_log2FC), max(plt_df$avg_log2FC)), expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(limits = c(0, max(plt_df$lp)), expand = expansion(mult = c(0.05, 0.01))) +
    labs(x=NULL,y=NULL) + 
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),legend.position="none",#legend.title=element_text(size=8),
          plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
          axis.title = element_text(size = 10, color = "black"))
  ggsave2(paste0(j,"_volcano.png"),width=3, height=3,device="png")
}


#FC vs FC
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
p2 <- ggplot(df, aes(b6FC,m564FC)) + geom_point(aes(colour=label),size=0.5)+#geom_point(fill=NA,colour=alpha("black",0.5),pch=21,size=3) +
  theme_classic()+
  scale_color_manual(values = c("up" = "#0D0887FF", "down" = "#F48849FF", "NS" = "grey80","diff"="black"))+
  scale_x_continuous(limits = c(-2.5,3), expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(limits = c(-2.5,3), expand = expansion(mult = c(0.01, 0.01))) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position="none",
        axis.ticks=element_blank(),axis.line=element_blank())+ #axis.text=element_blank(),
  geom_hline(yintercept=0,color=alpha("black",0.5))+geom_vline(xintercept=0,color=alpha("black",0.5))
if(length(gene.list)>0){p2 <- LabelPoints(plot = p2, points = gene.list, repel = TRUE,xnudge=0,ynudge=0,size=2.5,segment.size=0.1)}
p2
ggsave2("FCvsFC.png",width=2.7, height=2.7,device="png")


