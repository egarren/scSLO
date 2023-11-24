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
my.ttest <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}



obj2<-SLO_all.combined
top8 <- SLO_all.markers[SLO_all.markers$cluster %in% c(0:12),] %>% group_by(cluster) %>% top_n(8, avg_log2FC)
top6 <- SLO_all.markers[SLO_all.markers$cluster %in% c(0:12),] %>% group_by(cluster) %>% top_n(6, avg_log2FC)
top4 <- SLO_all.markers[SLO_all.markers$cluster %in% c(0:12),] %>% group_by(cluster) %>% top_n(4, avg_log2FC)
top2 <- SLO_all.markers[SLO_all.markers$cluster %in% c(0:12),] %>% group_by(cluster) %>% top_n(2, avg_log2FC)
k="all"

obj2<-obj
top8 <- obj.markers[obj.markers$cluster %in% c(0:5),] %>% group_by(cluster) %>% top_n(8, avg_log2FC)
top6 <- obj.markers[obj.markers$cluster %in% c(0:5),] %>% group_by(cluster) %>% top_n(6, avg_log2FC)
top4 <- obj.markers[obj.markers$cluster %in% c(0:5),] %>% group_by(cluster) %>% top_n(4, avg_log2FC)
top2 <- obj.markers[obj.markers$cluster %in% c(0:5),] %>% group_by(cluster) %>% top_n(2, avg_log2FC)
k="fdc"

# ## tSNE/UMAP
# # Visualization
DimPlot(obj2, reduction = "umap")
DimPlot(obj2, reduction = "umap",pt.size=0.001)+ #NoLegend()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
                   axis.text=element_blank(),axis.ticks=element_blank(),axis.line=element_blank())
ggsave2(paste0("umap_",k,".png"),width=7, height=5,device="png")

DimPlot(obj2, reduction = "umap",pt.size=0.1,label=T)+
  NoLegend()+ NoAxes()+theme(panel.border = element_blank())
ggsave2(paste0("umap.label_subset",k,".png"),width=5, height=5,device="png")

DimPlot(obj2, reduction = "umap",pt.size=0.1,split.by="Organ2",shuffle=T)+ 
  NoLegend() + NoAxes()+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
ggsave2(paste0("umap.byorgan.",k,".png"),width=10,height=2,device="png")

DimPlot(obj2, reduction = "umap",pt.size=0.1,split.by="Organ2",shuffle=T,raster=F)+ 
  NoLegend() + NoAxes()+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
ggsave2(paste0("umap.byorgan.rasteroff.",k,".png"),width=10,height=2,device="png")

DimPlot(obj2, reduction = "umap",pt.size=0.1,split.by="Species",shuffle=T)+ 
  NoLegend() + NoAxes()+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
ggsave2(paste0("umap.byspecies.",k,".png"),width=4,height=2.25,device="png")


## Cluster Analysis
Idents(obj2) <- "my.clusters2"
VlnPlot(obj2, features = c("CD274"), group.by = "my.clusters2", pt.size=0.0001)
ggsave2(paste0("vlnplot.cd274.cluster.",k,".png"),width=6, height=4,device="png")
# DefaultAssay(obj2) <- "integrated" 
obj2<-BuildClusterTree(obj2)
png(paste0("cluster.tree_",k,".png"),width=4,height=6,units="in",res=300)
PlotClusterTree(obj2,font=1)
dev.off()
# # Cluster heatmap
# DefaultAssay(obj2) <- "integrated"
DoHeatmap(object = obj2, features = top8$gene, label = F)  #slim.col.label to TRUE prints cluster IDs instead of cells, ,size=5,angle=45,hjust=0.5
ggsave2(paste0("cluster.heatmap_",k,".png"),width=6, height=10,device="png")
DoHeatmap(object = obj2, features = top4$gene, label = F)  #slim.col.label to TRUE prints cluster IDs instead of cells, ,size=5,angle=45,hjust=0.5
ggsave2(paste0("cluster.heatmap2_",k,".png"),width=6, height=6,device="png")
# DefaultAssay(obj2) <- "RNA"

#UMAP by cluster marker
FeaturePlot(object = obj2, features = "CXCL16")

gene.list<-c("MADCAM1","CXCL13","CR2","CCL19","CCL21","PDPN",
             "CD34","TNFSF13B","INMT","KLRK1",
             "CD3E","CD19","PTPRC","PECAM1","ITGAX","ITGAM","CD74","CD68",
             "MARCO","PTX3","SIGLEC1","PROX1","VWF","MYH11",
             "AGT","EPCAM","CD274")
gene.list<-c("MADCAM1","TNFSF11","CXCL13","CR2","CXCL12","CCL19","CCL21","PDPN","CXCL9","CXCL10",
             "CH25H","IL7","LEPR","NR4A1","LEPR","CD34","TNFSF13B","INMT","NKG2D",
             "CD4","CD8A","CD3E","CD19","PTPRC","PECAM1","ACTA2","ITGAX","ITGAM","CD74","CD68",
             "CXCR5","SDC1","MARCO","PTX3","SIGLEC1","PROX1","DCN","VWF","TNFSF13B","VCAM1","MYH11",
             "APOE","AGT","FBN1","PTX3","ATF3","CD5L","EPCAM","CD274")
gene.list<-c("CR2","GPC6","HTR1F","PRICKLE2","PAPPA","NRXN3","C3","IRF7","CD274")
gene.list<-c("CR2","CXCL13","MFGE8",
             "GPC6","HTR1F","PRICKLE2","PAPPA","NRXN3","SLC4A8","CSN2","PTGDS","CXCL14","C3","COL3A1",
             "TAGLN","IRF7","CXCL10","ISG15","CCL7","CD274")
gene.list<-c("CR2","CXCL13","C3","PAPPA","CRYM","TNFSF11","CXCL14","VWC2","FCGR2A","NRXN3","CXCL12","IL6")
gene.list<-c("CR2","CXCL13","C3","PAPPA","CRYM","RIMS1","FCGR2A","CD34","PDZD11","EPSTI1","CXCL14","IL6")
gene.list<-c("CR2","CXCL13","C3","GREM1","TNFSF11","RIMS1","VWC2","FCGR2A","CXCL14","IL6","CXCL12","ECM1")
gene.list<-c("TENM2","INMT","CCL19","CD34","CR2","ADAMDEC1") #new2 (stroma)
gene.list<-c("C3","PAPPA","CRYM","RIMS1","FCGR2A","CD34","PDZD11","EPSTI1") #fdc_0
FeaturePlot(object = obj2, features = "CD274", reduction = "umap",pt.size=0.001, order=T)
p<-FeaturePlot(object = obj2, features = gene.list,pt.size=0.001,order=T,
               cols = c("grey", "blue"), reduction = "umap",combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
}
cowplot::plot_grid(plotlist=p,ncol=3)
ggsave2(paste0("umap.clustermarkers_",k,".png"),width=5.5,height=4,device="png")
cowplot::plot_grid(plotlist=p,ncol=4)
ggsave2(paste0("umap.clustermarkers_",k,".png"),width=8.5,height=4.5,device="png")
cowplot::plot_grid(plotlist=p,ncol=9)
ggsave2(paste0("umap.clustermarkers_",k,".png"),width=20,height=6,device="png")
#vln cluster marker
p<-VlnPlot(obj2, features = gene.list,group.by = "my.clusters2",pt.size = 0.01, combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]]+ NoLegend()+theme(axis.title=element_blank())
}
cowplot::plot_grid(plotlist = p,ncol=5)
ggsave2(paste0("vlnplot.clustermarkers_",k,".png"),width=9, height=9,device="png")
# ggsave2("vlnplot.clustermarkers.png",width=9, height=4.5,device="png")
#vln cluster marker
p<-RidgePlot(obj2, features = gene.list,group.by = "my.clusters2", combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]]+ NoLegend()+theme(axis.title=element_blank())
}
cowplot::plot_grid(plotlist = p,ncol=5)
ggsave2(paste0("ridgeplot.clustermarkers_",k,".png"),width=9, height=6,device="png")


#Vln plot comparison of genes
gene.list<-c("VCAM1","CD81","CD82","TNFRSF9","MFGE8","ICAM1","DLL1","CD274","CD40","PDCD1LG2")
data<-FetchData(obj2,vars=gene.list,slot="data")
long_data<-melt(data)
long_data$value[long_data$value==0]<-NA
long_data$variable<-reorder(long_data$variable,-long_data$value,mean)
ggplot(long_data,aes(x=variable,y=value))+geom_violin()+geom_jitter(size=0.01,color=adjustcolor("black",alpha.f=0.2))

#pheatamp
obj2@meta.data$group_by<-paste0(obj2@meta.data$Species,"-",obj2@meta.data$my.clusters2)
Idents(obj2)<-"group_by"
df<-AverageExpression(obj2)$RNA
df<-df[gene.list,]
df<-as.data.frame(df)
meta<-unique(obj2@meta.data[,c("Species","my.clusters2","group_by")])
meta$group_by<-gsub("_","-",meta$group_by)
rownames(meta)<-meta$group_by
meta$group_by<-NULL
meta$Species[meta$Species=="Hs"]<-"Human"
meta$Species[meta$Species=="Mus"]<-"Mouse"
colnames(meta)<-c("Species","Cluster")
p <- pheatmap(df,
              scale="row",
              cluster_rows = TRUE,cluster_cols = TRUE,show_colnames = F, treeheight_row=3,treeheight_col=3, #clustering_method = "ward.D2",
              # color = rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#FFFFFF"))(100)),
              color = colorRampPalette(c("#45628f", "#F7F7F7", "#de425b"))(50),
              fontsize = 8,fontsize_col = 8,fontsize_row=8,
              annotation_col=meta,annotation_names_col=F,annotation_names_row=F,
              annotation_colors=list(Species=c("Human"=viridis(3)[1],"Mouse"=viridis(3)[2]),
                                               Cluster=c("MRC"=hue_pal()(6)[1],"TRC"=hue_pal()(6)[2],
                                                    "INMT_TRC"=hue_pal()(6)[3],"CD34_SC"=hue_pal()(6)[4],
                                                    "FDC"=hue_pal()(6)[5],"MRC_ADAMDEC1"=hue_pal()(6)[6])),
              # annotation_row=annot.row,
              silent = T)
png("heatmap.png", width = 5, height = 2.5, res = 200,units="in")
grid.arrange(p$gtable)
dev.off()

#species expression correlation
obj2@meta.data$group_by<-paste0(obj2@meta.data$Species,"-",obj2@meta.data$my.clusters2)
Idents(obj2)<-"group_by"
df<-AverageExpression(obj2)$RNA
res<-cor(as.matrix(df))
res<-as.data.frame(res)
res<-res[rownames(res) %like% "Hs",colnames(res) %like% "Mus"]
colnames(res)<-gsub("Mus-","",colnames(res))
rownames(res)<-gsub("Hs-","",rownames(res))
colnames(res)<-gsub("-","_",colnames(res))
rownames(res)<-gsub("-","_",rownames(res))
order<-c("MRC","TRC","INMT_TRC","CD34_SC","FDC","MRC_ADAMDEC1")
res<-res[order,order]
png("corrplot.png", width = 4, height = 4, res = 200,units="in")
corrplot(as.matrix(res),
         is.corr = F,
         outline = T,
         col = rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#FFFFFF"))(100)),
         tl.col = "black",
         na.label = " ")
dev.off()
png("corrplot_legend.png", width = 6, height = 6, res = 200,units="in")
corrplot(as.matrix(res),
         is.corr = F,
         outline = T,
         col = rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#FFFFFF"))(100)),
         tl.col = "black",
         na.label = " ")
dev.off()

#scatterplot
obj2@meta.data$group_by<-paste0(obj2@meta.data$Species,"-",obj2@meta.data$my.clusters2)
Idents(obj2)<-"group_by"
df<-AverageExpression(obj2)$RNA
df<-as.data.frame(df)
colnames(df)<-gsub("-","_",colnames(df))
# df$diff<-abs(df$Mus_FDC-df$Hs_FDC)
df$diff2<-"NS"
df$diff2[df$Mus_FDC-df$Hs_FDC > 5]<-"UP"
df$diff2[df$Mus_FDC-df$Hs_FDC < -5]<-"DOWN"
df$gene<-rownames(df)
cor.test(df$Mus_FDC,df$Hs_FDC)
# df<- df %>% mutate(gene_label = ifelse(diff>5 , gene, NA))
gene.list<-c("CLU","CR2","CXCL13","FTH1","MFGE8","APOE","TMSB10","MADCAM1","VCAM1")
df<- df %>% mutate(gene_label = ifelse(gene %in% gene.list , gene, NA))
options(scipen = 999)
ggplot(df,aes(Mus_FDC,Hs_FDC,color=diff2,label=gene_label))+geom_point()+theme_classic()+ 
  geom_text_repel(color="black",size = 3 ) +
  scale_color_manual(values=c(viridis(3)[1],"#999999", viridis(3)[2]))+
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10')+
  xlab("Mouse") + ylab("Human")+ 
  theme(legend.position = "none")
ggsave2("fdc_compare.png",width=4, height=3.5,device="png")


#dot plot
top6$gene<-make.unique(top6$gene,sep="--")
DotPlot(obj2, features = top6$gene)+ 
  coord_flip()+theme(legend.title = element_blank(),axis.title=element_blank())#+ RotatedAxis()
ggsave2(paste0("dotplot_",k,".png"),width=12, height=6,device="png")


#by condition
FeaturePlot(object = obj2, features = "CD274", cols = c("grey", "blue"), 
            split.by="Organ2",pt.size=0.001, order=T)
ggsave2("umap.cd274.organ2.png",width=20, height=3,device="png")

# Cluster frequency comparison
ggplot(obj2@meta.data, aes(x = my.clusters2,fill = Organ2)) +  
  geom_bar(aes(y = (..count..)/sum(..count..)),position='dodge')+ 
  scale_y_continuous(labels = percent)+ 
  labs(y="Frequency", x = "Cluster")
ggsave2("clusterfreq.png",width=4, height=2.5,device="png")
cluster.Counts <- table(obj2@meta.data$my.clusters2,obj2@meta.data$Organ2)
cluster.prop <- as.data.frame(scale(cluster.Counts,scale=colSums(cluster.Counts),center=FALSE)*100) 
ggplot(cluster.prop, aes(fill=Var1,y=Freq, x=Var2,alluvium=Var1,stratum=Var1)) + 
  geom_lode()+geom_flow()+geom_stratum(alpha=0) +theme_classic()+
  # scale_x_discrete(labels= c("WT","564Igi"))+
  theme(legend.title = element_blank(),axis.title=element_blank()) 
  #+ theme(legend.position = "none")
ggsave2("cluster.prop.flow.png",width=4, height=4,device="png")
write.csv(cluster.prop,file="rawvalue.csv")
#byspecies
cluster.Counts <- table(obj2@meta.data$my.clusters2,obj2@meta.data$Species)
cluster.prop <- as.data.frame(scale(cluster.Counts,scale=colSums(cluster.Counts),center=FALSE)*100) 
ggplot(cluster.prop, aes(fill=Var1,y=Freq, x=Var2,alluvium=Var1,stratum=Var1)) + 
  geom_lode()+geom_flow()+geom_stratum(alpha=0) +theme_classic()+
  # scale_x_discrete(labels= c("WT","564Igi"))+
  theme(legend.title = element_blank(),axis.title=element_blank()) +
  coord_flip()
#+ theme(legend.position = "none")
ggsave2("cluster.prop.flow.species.png",width=6, height=1.5,device="png")


###Volcano
for(j in ls(pattern=".DE")){
  res<-get(j)
  # res<-SLO_all.markers[SLO_all.markers$cluster==5,]
  thresh_p_val_adj <- 1e-10
  thresh_lfc <-0.3
  plt_df<- res %>% rownames_to_column(var = "gene") %>% 
    # filter(!(grepl("Rps", gene) | grepl("Rpl", gene)| grepl("mt.", gene)| grepl("H2.", gene))) %>% 
    mutate(up_in = ifelse(p_val_adj >= thresh_p_val_adj, "NS", ifelse(avg_log2FC > thresh_lfc, "FDC", ifelse(avg_log2FC < -thresh_lfc, "TRC", "NS"))),
           rank_pval = rank(p_val_adj), rank_lfc_inc = rank(avg_log2FC), rank_lfc_dec = rank(-abs(avg_log2FC)),
           lp = -log10(p_val_adj)) %>% 
    arrange(-abs(avg_log2FC))
  if (sum(plt_df$up_in != "NS") >= 30) {
    plt_df <- plt_df %>% mutate(gene_label = ifelse(up_in != "NS" & (rank_pval < 10 | rank_lfc_inc < 10 | rank_lfc_dec < 10), gene, NA))
  } else {
    plt_df <- plt_df %>% mutate(gene_label = ifelse((up_in != "NS" | rank_pval < 10 | rank_lfc_inc < 10 | rank_lfc_dec < 10), gene, NA))
  }
  plt_df$gene_label[plt_df$gene=="CD274"]<-"CD274"
  plt_df$lp <- pmin(plt_df$lp, 300)
  plt_df$avg_log2FC <- pmin(plt_df$avg_log2FC, 3)
  ggplot(plt_df, aes(x = avg_log2FC, y = lp, color = up_in, label = gene_label)) +
    geom_point(size = 0.5) +
    geom_text_repel(color="black",size = 3, segment.size = 0.2,seed = 1, show.legend=F, box.padding=0.5) +
    geom_hline(yintercept=10,linetype=2,size=0.2)+
    geom_vline(xintercept=0.3,linetype=2,size=0.2)+
    geom_vline(xintercept= -0.3,linetype=2,size=0.2)+
    scale_color_manual(values = c("FDC" = hue_pal()(6)[5], "TRC" = hue_pal()(6)[2], "NS" = "grey80")) +
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


#custom volcano
j="pdl1.DE"
res<-get(j)
thresh_p_val_adj <- 1e-20
thresh_lfc <-0.3
plt_df<- res %>% rownames_to_column(var = "gene") %>% 
  
  mutate(up_in = ifelse(p_val_adj >= thresh_p_val_adj, "NS", ifelse(avg_log2FC > thresh_lfc, "PD-L1 high", ifelse(avg_log2FC < -thresh_lfc, "PD-L1 neg", "NS"))),
         rank_pval = rank(p_val_adj), rank_lfc_inc = rank(avg_log2FC), rank_lfc_dec = rank(-abs(avg_log2FC)),
         lp = -log10(p_val_adj)) %>% 
  arrange(-abs(avg_log2FC))
table(plt_df$up_in)
plt_df$gene_label[plt_df$gene=="CD274"]<-"CD274"
plt_df$lp <- pmin(plt_df$lp, 300)
plt_df$avg_log2FC <- pmin(plt_df$avg_log2FC, 3)
ggplot(plt_df, aes(x = avg_log2FC, y = lp, color = up_in, label = gene_label)) +
  geom_point(size = 0.5) +
  geom_text_repel(color="black",size = 3, segment.size = 0.2,seed = 1, show.legend=F, box.padding=0.5) +
  geom_hline(yintercept=20,linetype=2,size=0.2)+
  geom_vline(xintercept=0.3,linetype=2,size=0.2)+
  geom_vline(xintercept= -0.3,linetype=2,size=0.2)+
  scale_color_manual(values = c("PD-L1 high" = "darkgreen", "PD-L1 neg" = "#D95F02", "NS" = "grey80")) +
  scale_x_continuous(limits = c(-max(plt_df$avg_log2FC), max(plt_df$avg_log2FC)), expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(limits = c(0, max(plt_df$lp)), expand = expansion(mult = c(0.05, 0.01))) +
  labs(x=NULL,y=NULL) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),legend.position="none",#legend.title=element_text(size=8),
        plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 10, color = "black"))
ggsave2(paste0(j,"_volcano.png"),width=3, height=2.5,device="png")


#custom volcano
j="fdc_cells.DE"
res<-get(j)
thresh_p_val_adj <- 1e-10
thresh_lfc <-0.3
plt_df$gene_label[plt_df$gene=="CD274"]<-"CD274"
gene.list<-c("VCAM1","CD82","TNFRSF9","MFGE8","DLL1","CD274","CD40","CR2",
             "CXCL13","CCL19","CCL21","INMT","CD82","FCER2","CTNNA3","BRINP3",
             "EGLN3","PI16",
             "IL6","MALT1","IL12A","PRNP","RAC2","CD80","CRYIB","TNFSF11","SOCS1","CD86","CEBPB","TCF7")
plt_df<- plt_df %>% mutate(gene_label = ifelse(gene %in% gene.list , gene, NA))
gene.list2<-unique(m_df$gene_symbol[m_df$gs_id=="M16253"]) #M18810
plt_df <- plt_df %>% mutate(label = ifelse(gene %in% gene.list2, "sig", "notsig"))
plt_df<-plt_df[order(plt_df$label),]
plt_df$lp <- pmin(plt_df$lp, 300)
plt_df$avg_log2FC <- pmin(plt_df$avg_log2FC, 3)
ggplot(plt_df, aes(x = avg_log2FC, y = lp, color = label, label = gene_label,size=label)) +
  geom_point() +
  geom_text_repel(color="black",size = 2, segment.size = 0.2,seed = 1, show.legend=T, box.padding=0.5) +
  geom_hline(yintercept=10,linetype=2,size=0.2)+
  geom_vline(xintercept=0.3,linetype=2,size=0.2)+
  geom_vline(xintercept= -0.3,linetype=2,size=0.2)+
  scale_size_manual(values=c(0.5,1.5))+
  scale_color_manual(values = c("sig" = hue_pal()(6)[5], "notsig" = "grey80")) +
  scale_x_continuous(limits = c(-max(plt_df$avg_log2FC), max(plt_df$avg_log2FC)), expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(limits = c(0, max(plt_df$lp)), expand = expansion(mult = c(0.05, 0.01))) +
  labs(x=NULL,y=NULL) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text=element_text(family="Arial"),
        panel.background = element_blank(),legend.position="none",#legend.title=element_text(size=8),
        plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 10, color = "black"))
ggsave2(paste0(j,"_volcano.png"),width=3, height=2.5,device="png")



#FC vs FC
hum_fdc_cells.DE$gene<-rownames(hum_fdc_cells.DE)
mus_fdc_cells.DE$gene<-rownames(mus_fdc_cells.DE)
df<-left_join(hum_fdc_cells.DE[,c("avg_log2FC","gene")],mus_fdc_cells.DE[,c("avg_log2FC","gene")],by="gene",keep=F)
colnames(df)<-c("humFC","gene","musFC")
rownames(df)<-df$gene
DE<-hum_fdc_cells.DE
df <- df %>% mutate(label = ifelse(gene %in% gene.list2, "sig", "notsig"))
df<-df[order(df$label),]
p2 <- ggplot(df, aes(musFC,humFC)) + geom_point(aes(colour=label,size=label))+#geom_point(fill=NA,colour=alpha("black",0.5),pch=21,size=3) +
  theme_classic()+
  scale_size_manual(values=c(0.5,1.5))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position="none",
        axis.ticks=element_blank(),axis.line=element_blank())+ #axis.text=element_blank(),
  geom_hline(yintercept=0,color=alpha("black",0.5))+geom_vline(xintercept=0,color=alpha("black",0.5))+#geom_abline(intercept = 0, slope = 1)+
  scale_color_manual(values = c("sig" = hue_pal()(6)[5], "notsig" = "grey70"))
if(length(gene.list)>0){p2 <- LabelPoints(plot = p2, points = gene.list, repel = TRUE,xnudge=0,ynudge=0,size=2.5,segment.size=0.1)}
p2
ggsave2("FCvsFC.png",width=3, height=3,device="png")
