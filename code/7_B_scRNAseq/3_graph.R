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


obj2<-obj
top8 <- obj.markers %>% group_by(cluster) %>% top_n(8, avg_log2FC)
top6 <- obj.markers %>% group_by(cluster) %>% top_n(6, avg_log2FC)
top4 <- obj.markers %>% group_by(cluster) %>% top_n(4, avg_log2FC)
top2 <- obj.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)
k="ALL"
Idents(obj2)<-"my.clusters2"

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

DimPlot(obj2, reduction = "umap",pt.size=0.1,split.by="GSE",shuffle=T)+ 
  NoLegend() + NoAxes()+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
ggsave2(paste0("umap.byGSE.",k,".png"),width=6,height=2,device="png")

DimPlot(obj2, reduction = "umap",pt.size=0.1,split.by="Species",shuffle=T)+ 
  NoLegend() + NoAxes()+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
ggsave2(paste0("umap.byspecies.",k,".png"),width=4,height=2.25,device="png")

FeaturePlot(obj, "CD274",order=T)+ 
  NoAxes()+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
                 plot.title = element_blank()) & 
  scale_color_gradientn(colors = rocket(n = 10, direction = -1))
ggsave2("pdl1_umap.png",width=4.5, height=4,device="png")

FeaturePlot(obj, "PDCD1",order=T)+ 
  NoAxes()+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
                 plot.title = element_blank()) & 
  scale_color_gradientn(colors = rocket(n = 10, direction = -1))
ggsave2("pd1_umap.png",width=4.5, height=4,device="png")

## Cluster Analysis
Idents(obj2) <- "my.clusters2"
VlnPlot(obj2, features = c("PDCD1"), group.by = "my.clusters2", pt.size=0.0001)
VlnPlot(obj2, features = c("CD274"), group.by = "my.clusters2", pt.size=0.0001)
ggsave2(paste0("vlnplot.cd274.cluster.",k,".png"),width=6, height=4,device="png")
DoHeatmap(object = obj2, features = top8$gene, label = F)  #slim.col.label to TRUE prints cluster IDs instead of cells, ,size=5,angle=45,hjust=0.5
ggsave2(paste0("cluster.heatmap_",k,".png"),width=6, height=10,device="png")
DoHeatmap(object = obj2, features = top4$gene, label = F)  #slim.col.label to TRUE prints cluster IDs instead of cells, ,size=5,angle=45,hjust=0.5
ggsave2(paste0("cluster.heatmap2_",k,".png"),width=6, height=6,device="png")

#UMAP by cluster marker
FeaturePlot(object = obj2, features = "CXCL16")

gene.list<-c("SELL","SERPINA9","VIM","MKI67","JCHAIN","MYC",
             "IGHG1","AICDA","TRAF1","S1PR2","CAMK1D","TOMM6",
             "MZB1","PLAC8")

FeaturePlot(object = obj2, features = "CD274", reduction = "umap",pt.size=0.001, order=T)
p<-FeaturePlot(object = obj2, features = gene.list,pt.size=0.001,order=T,
               cols = c("grey", "blue"), reduction = "umap",combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
}
cowplot::plot_grid(plotlist=p,ncol=5)
ggsave2(paste0("umap.clustermarkers_",k,".png"),width=11,height=7,device="png")


#Vln plot comparison of genes
gene.list<-c("VCAM1","CD81","CD82","TNFRSF9","MFGE8","ICAM1","DLL1","CD274","CD40","PDCD1LG2")
data<-FetchData(obj2,vars=gene.list,slot="data")
long_data<-melt(data)
long_data$value[long_data$value==0]<-NA
long_data$variable<-reorder(long_data$variable,-long_data$value,mean)
ggplot(long_data,aes(x=variable,y=value))+geom_violin()+geom_jitter(size=0.01,color=adjustcolor("black",alpha.f=0.2))

#pheatamp
Idents(obj2)<-"my.clusters2"
df<-AverageExpression(obj2,assays="RNA",features=gene.list,slot="data")$RNA
df<-as.data.frame(df)
colnames(df)<-gsub("-","_",colnames(df))
meta<-data.frame("Cluster"=unique(obj2@meta.data[,c("my.clusters2")]))
rownames(meta)<-meta$Cluster
p <- pheatmap(df,
              scale="row",
              cluster_rows = TRUE,cluster_cols = TRUE,show_colnames = F, treeheight_row=3,treeheight_col=3, #clustering_method = "ward.D2",
              color = colorRampPalette(c("#45628f", "#F7F7F7", "#de425b"))(50),
              fontsize = 8,fontsize_col = 8,fontsize_row=8,
              annotation_col=meta,annotation_names_col=F,annotation_names_row=F,
              annotation_colors=list(Cluster=c("Naive"=hue_pal()(14)[1],"GC_1"=hue_pal()(14)[2],
                                                "MBC_1"=hue_pal()(14)[3],"DZ"=hue_pal()(14)[4],
                                                "ASC_1"=hue_pal()(14)[5],"Activated"=hue_pal()(14)[6],
                                               "PC"=hue_pal()(14)[7],"LZ"=hue_pal()(14)[8],
                                               "MBC_2"=hue_pal()(14)[9],"GC_2"=hue_pal()(14)[10],
                                               "ASC_2"=hue_pal()(14)[11],"MZ_1"=hue_pal()(14)[12],
                                               "ASC_3"=hue_pal()(14)[13],"MZ_2"=hue_pal()(14)[14])),
              silent = T)
png("heatmap.png", width = 5, height = 2.5, res = 200,units="in")
grid.arrange(p$gtable)
dev.off()

#species expression correlation
obj2@meta.data$group_by<-paste0(obj2@meta.data$Species,"-",obj2@meta.data$my.clusters2)
Idents(obj2)<-"group_by"
df<-AverageExpression(obj2,assays="RNA",slot="data")$RNA
res<-cor(as.matrix(df))
res<-as.data.frame(res)
res<-res[rownames(res) %like% "Human",colnames(res) %like% "Mouse"]
colnames(res)<-gsub("Mouse-","",colnames(res))
rownames(res)<-gsub("Human-","",rownames(res))
colnames(res)<-gsub("-","_",colnames(res))
rownames(res)<-gsub("-","_",rownames(res))
order<-as.data.frame(table(obj2@meta.data$my.clusters2))$Var1
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


#dot plot
Idents(obj2)<-"my.clusters2"
top6$gene<-make.unique(top6$gene,sep="--")
DotPlot(obj2, features = top6$gene)+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  coord_flip()+theme(legend.title = element_blank(),axis.title=element_blank())#+ RotatedAxis()
ggsave2(paste0("dotplot_",k,".png"),width=8, height=12,device="png")

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
  theme(legend.title = element_blank(),axis.title=element_blank()) 
  ggsave2("cluster.prop.flow.png",width=4, height=4,device="png")
write.csv(cluster.prop,file="rawvalue.csv")
#byspecies
cluster.Counts <- table(obj2@meta.data$my.clusters2,obj2@meta.data$Species)
cluster.prop <- as.data.frame(scale(cluster.Counts,scale=colSums(cluster.Counts),center=FALSE)*100) 
ggplot(cluster.prop, aes(fill=Var1,y=Freq, x=Var2,alluvium=Var1,stratum=Var1)) + 
  geom_lode()+geom_flow()+geom_stratum(alpha=0) +theme_classic()+
  theme(legend.title = element_blank(),axis.title=element_blank()) +
  coord_flip()
ggsave2("cluster.prop.flow.species.png",width=6, height=1.5,device="png")


###Volcano
for(j in ls(pattern=".DE")){
  res<-get(j)
  thresh_p_val_adj <- 1e-10
  thresh_lfc <-1
  plt_df<- res %>% rownames_to_column(var = "gene") %>% 
    mutate(up_in = ifelse(p_val_adj >= thresh_p_val_adj, "NS", ifelse(avg_log2FC > thresh_lfc, "up", ifelse(avg_log2FC < -thresh_lfc, "down", "NS"))),
           rank_pval = rank(p_val_adj), rank_lfc_inc = rank(avg_log2FC), rank_lfc_dec = rank(-avg_log2FC),
           lp = -log10(p_val_adj)) %>% 
    arrange(-abs(avg_log2FC))
  if (sum(plt_df$up_in != "NS") >= 30) {
    plt_df <- plt_df %>% mutate(gene_label = ifelse(up_in != "NS" & (rank_pval < 20 | rank_lfc_inc < 20 | rank_lfc_dec < 20), gene, NA))
  } else {
    plt_df <- plt_df %>% mutate(gene_label = ifelse((up_in != "NS" | rank_pval < 20 | rank_lfc_inc < 20 | rank_lfc_dec < 20), gene, NA))
  }
  unique(plt_df$gene_label)
  plt_df$gene[plt_df$rank_lfc_dec<20]
  plt_df[plt_df$gene=="CD274",]
  col_up<-"red"
  col_down<-"red"
  if(grepl("pdl1",j)){
    plt_df$gene_label[plt_df$gene=="CD274"]<-"CD274"
    col_up<-"darkgreen"
    col_down<-"#D95F02"}
  if(grepl("pd1",j)){
    plt_df$gene_label[plt_df$gene=="PDCD1"]<-"PDCD1"
    col_up<-"#0D0887FF"
    col_down<-"#F48849FF"}
  if(j=="asc_gc.DE"){
    gene.list<-c("CD274","SLPI","JCHAIN","IGHG3","LTB","FCMR","LCK","BIK","SPIB","SERPINA9",
                 "SUGT","GCHFR","TCL1A","MARCKSL1","EDEM1","GPHN","ICOSLG","CAMK1D","ADGRL2",
                 "IGLC2","RAPH1","GSN","RELN","RFLNB","IGHG4","MAML2","VSIR","SPON1",
                 "STAG3","EPHX1","VNN2","TMEM108")
    plt_df<- plt_df %>% mutate(gene_label = ifelse(gene %in% gene.list , gene, NA))
    col_up<-hue_pal()(14)[11]
    col_down<-hue_pal()(14)[2]}
  if(j=="asc2_dz.DE"){
    gene.list<-c("CPEB3","IRF4","STMN1","HMGB1","AICDA","CD274","PTTG1","DEK",
                 "TP53INP1","CPEB2","DPED1","PKD2L2","SDC1","SSPN","PTGES","CDK1",
                 "TPX2","CDC20","H4C3","AURKB","ADGRL2","DGKG","CD274","ERN1",
                 "AFF1","PKD2L2","CPED1","DNM3","TRAM2","HMGN2","PARP8","HMGN1")
    plt_df<- plt_df %>% mutate(gene_label = ifelse(gene %in% gene.list , gene, NA))
    col_up<-hue_pal()(14)[11]
    col_down<-hue_pal()(14)[4]}
  if(j=="lz_dz.DE"){
    gene.list<-c("CDK1","H4C3","PCLAF","TOP2A","CCNA2","TYMS","CLSPN",
                 "BACH2","SOX5","AICDA","MKI67","EFNB1","MARCHF1","PDZD2",
                 "AFF3","SETBP1","GNG7","NTNG2","FGFR2","CDH13","CD274")
    plt_df<- plt_df %>% mutate(gene_label = ifelse(gene %in% gene.list , gene, NA))
    col_up<-hue_pal()(14)[8]
    col_down<-hue_pal()(14)[4]}
  if(j=="pd1.DE"){
    gene.list<-c("PDCD1","WIPF3","CTLA4","APOE","LCK","EGR2")
    plt_df<- plt_df %>% mutate(gene_label = ifelse(gene %in% gene.list , gene, NA))}
  if(j=="pdl1.DE"){
    gene.list<-c("IGHG3","CD274","TMSB4X","MALAT1","HMGB1","SLPI","JCHAIN","PON3",
                 "EDEM2","CPEB3","SDC1","TP53INP1","IRF4","PRG2","PKD2L2","CRELD2",
                 "SERF2","PDIA4","ATP5MG","HLA.A","PPIA","TCL1A","IGH3","IGHD",
                 "IGLC3","EIF5A","PCNP","HMGA1","CD27","H4C3","PEBP1","MZTB",
                 "SET","UBL5","HMGN2","EIF4A2","SUMO2","LDHB","PRR13","IFI16","LDHA")
    plt_df<- plt_df %>% mutate(gene_label = ifelse(gene %in% gene.list , gene, NA))}
  maxFC<-3.5
  plt_df$lp <- pmin(plt_df$lp, 300)
  plt_df$avg_log2FC <- pmin(plt_df$avg_log2FC, maxFC)
  plt_df$avg_log2FC <- pmax(plt_df$avg_log2FC, -maxFC)
  ggplot(plt_df, aes(x = avg_log2FC, y = lp, color = up_in, label = gene_label)) +
    geom_point(size = 0.5) +
    geom_text_repel(color="black",size = 3, segment.size = 0.2,seed = 1, show.legend=F, box.padding=0.1) +
    geom_hline(yintercept=10,linetype=2,size=0.2)+
    geom_vline(xintercept=1,linetype=2,size=0.2)+
    geom_vline(xintercept= -1,linetype=2,size=0.2)+
    scale_color_manual(values = c("up"=col_up, "down" = col_down, "NS" = "grey80")) +
    scale_x_continuous(limits = c(-maxFC, maxFC), expand = expansion(mult = c(0.01, 0.01))) +
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
for(k in c("pd1","pdl1")){
  for(i in ls(pattern=paste0(k,".DE"))){
    x<-get(i)
    for(j in ls(pattern=paste0(k,".DE"))){
      if(i !=j){
        y<-get(j)
        x$gene<-rownames(x)
        y$gene<-rownames(y)
        df<-left_join(x[,c("avg_log2FC","gene")],y[,c("avg_log2FC","gene")],by="gene",keep=F)
        colnames(df)<-c("x","gene","y")
        rownames(df)<-df$gene
        df$diff= abs(df$x-df$y)
        df$sum = df$x + df$y
        df<-df %>% mutate(rank_diff = rank(-diff), rank_sum = rank(sum), rank_sum_dec = rank(-sum)) 
        gene.list<-df$gene[df$rank_diff<10 | df$rank_sum <5 |df$rank_sum_dec<10]
        df$label<-"NS"
        df$label[df$x>2 & df$y>2]<-"up"
        df$label[df$x< -2 & df$y< -2]<-"down"
        df$label[df$diff>2]<-"diff"
        df<-df[order(-df$rank_diff),]
        if(k=="pdl1"){
          col_up<-"darkgreen"
          col_down<-"#D95F02"}
        if(k=="pd1"){
          col_up<-"#0D0887FF"
          col_down<-"#F48849FF"}
        p2 <- ggplot(df, aes(x,y)) + geom_point(aes(colour=label),size=0.5)+#geom_point(fill=NA,colour=alpha("black",0.5),pch=21,size=3) +
          theme_classic()+
          scale_color_manual(values = c("up" = col_up, "down" = col_down, "NS" = "grey80","diff"="black"))+
          labs(x=i,y=j)+
          theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position="none",
                axis.ticks=element_blank(),axis.line=element_blank())+ #axis.text=element_blank(),
          geom_hline(yintercept=0,color=alpha("black",0.5))+geom_vline(xintercept=0,color=alpha("black",0.5))
        if(length(gene.list)>0){p2 <- LabelPoints(plot = p2, points = gene.list, repel = TRUE,xnudge=0,ynudge=0,size=2.5,segment.size=0.1)}
        p2
        ggsave2(paste0(i,"_",j,"_FCvsFC.png"),width=2.5, height=2.5,device="png")
      }
    }
  }
}
