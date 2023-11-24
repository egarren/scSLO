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



obj2<-T_all.combined
top8 <- T_all.markers[T_all.markers$cluster %in% c(0:12),] %>% group_by(cluster) %>% top_n(8, avg_log2FC)
top6 <- T_all.markers[T_all.markers$cluster %in% c(0:12),] %>% group_by(cluster) %>% top_n(6, avg_log2FC)
top4 <- T_all.markers[T_all.markers$cluster %in% c(0:12),] %>% group_by(cluster) %>% top_n(4, avg_log2FC)
top2 <- T_all.markers[T_all.markers$cluster %in% c(0:12),] %>% group_by(cluster) %>% top_n(2, avg_log2FC)
k="all"

obj2<-obj
top8 <- obj.markers[obj.markers$cluster %in% c(0:5),] %>% group_by(cluster) %>% top_n(8, avg_log2FC)
top6 <- obj.markers[obj.markers$cluster %in% c(0:5),] %>% group_by(cluster) %>% top_n(6, avg_log2FC)
top4 <- obj.markers[obj.markers$cluster %in% c(0:5),] %>% group_by(cluster) %>% top_n(4, avg_log2FC)
top2 <- obj.markers[obj.markers$cluster %in% c(0:5),] %>% group_by(cluster) %>% top_n(2, avg_log2FC)
k="TFR"

# ## tSNE/UMAP
# # Visualization
Idents(obj2)<-"my.clusters2"
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

#UMAP by cluster marker
gene.list<-c("Foxp3","Sostdc1","Lmo4","Ext1","Klf2","Slamf7","Ccl5","Ifit3") #fdc_0
gene.list<-c("Foxp3","Sostdc1","Tcf7","Ext1","Klf2","Slamf7","Ccl5","Ifit3") #fdc_0
gene.list<-c("Tpi1","Selplg","C1qbp") #fdc_0
p<-FeaturePlot(object = obj2, features = gene.list,pt.size=0.001,order=T,
               cols = c("grey", "blue"), reduction = "umap",combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
}
cowplot::plot_grid(plotlist=p,ncol=3)
ggsave2(paste0("umap.clustermarkers_",k,".png"),width=5.5,height=2.1,device="png")
cowplot::plot_grid(plotlist=p,ncol=4)
ggsave2(paste0("umap.clustermarkers_",k,".png"),width=8,height=4.5,device="png")

#vln cluster marker
VlnPlot(obj2, features = "Pdcd1",group.by = "my.clusters2",pt.size = 0.01)
ggsave2(paste0("vlnplot.pd1_",k,".png"),width=3, height=4,device="png")
VlnPlot(obj2, features = "Pdcd1",group.by = "my.clusters2",pt.size = 0)
ggsave2(paste0("vlnplot.pd1_",k,".png"),width=5, height=3,device="png")
VlnPlot(obj2, features = "Pdcd1",group.by = "my.clusters2",pt.size = 0)+
  stat_summary(fun.data=mean_cl_normal, mult=1,geom="pointrange", color="black")
ggsave2(paste0("vlnplot.pd1_mean_",k,".png"),width=3, height=4,device="png")
VlnPlot(obj2, features = "Pdcd1",group.by = "my.clusters2",split.by="condition",pt.size = 0.01)
ggsave2(paste0("vlnplotsplit.pd1_",k,".png"),width=4, height=3,device="png")
p<-VlnPlot(obj2, features = gene.list,group.by = "my.clusters2",pt.size = 0.01, combine=F)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]]+ NoLegend()+theme(axis.title=element_blank())
}
cowplot::plot_grid(plotlist = p,ncol=5)
ggsave2(paste0("vlnplot.clustermarkers_",k,".png"),width=9, height=9,device="png")
#vln cluster marker
RidgePlot(obj2, features = "Pdcd1",group.by = "my.clusters2")
ggsave2(paste0("ridgeplot.pd1_",k,".png"),width=4, height=4,device="png")
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
obj2@meta.data$group_by<-paste0(obj2@meta.data$condition,"-",obj2@meta.data$my.clusters2)
Idents(obj2)<-"group_by"
df<-AverageExpression(obj2)$RNA
df<-df[gene.list,]
df<-as.data.frame(df)
meta<-unique(obj2@meta.data[,c("condition","my.clusters2","group_by")])
meta$group_by<-gsub("_","-",meta$group_by)
rownames(meta)<-meta$group_by
meta$group_by<-NULL
meta$condition[meta$condition=="Hs"]<-"Human"
meta$condition[meta$condition=="Mus"]<-"Mouse"
colnames(meta)<-c("condition","Cluster")
p <- pheatmap(df,
              scale="row",
              cluster_rows = TRUE,cluster_cols = TRUE,show_colnames = F, treeheight_row=3,treeheight_col=3, #clustering_method = "ward.D2",
              color = colorRampPalette(c("#45628f", "#F7F7F7", "#de425b"))(50),
              fontsize = 8,fontsize_col = 8,fontsize_row=8,
              annotation_col=meta,annotation_names_col=F,annotation_names_row=F,
              annotation_colors=list(condition=c("Human"=viridis(3)[1],"Mouse"=viridis(3)[2]),
                                               Cluster=c("MRC"=hue_pal()(6)[1],"TRC"=hue_pal()(6)[2],
                                                    "INMT_TRC"=hue_pal()(6)[3],"CD34_SC"=hue_pal()(6)[4],
                                                    "FDC"=hue_pal()(6)[5],"MRC_ADAMDEC1"=hue_pal()(6)[6])),
              silent = T)
png("heatmap.png", width = 5, height = 2.5, res = 200,units="in")
grid.arrange(p$gtable)
dev.off()

#condition expression correlation
obj2@meta.data$group_by<-paste0(obj2@meta.data$condition,"-",obj2@meta.data$my.clusters2)
Idents(obj2)<-"group_by"
df<-AverageExpression(obj2)$RNA
res<-cor(as.matrix(df))
res<-as.data.frame(res)
res<-res[rownames(res) %like% "AID",colnames(res) %like% "m564"]
colnames(res)<-gsub("m564-","",colnames(res))
rownames(res)<-gsub("AID-","",rownames(res))
order<-c("TFR","Sostdc1","TFH-Tcf1","TFH-Exhausted","TFH-Activated","TFH-CM","TFH-Effector","TFH-ISG")
order<-c("TFR-1","TFR-2","TFR-3")
res<-res[order,order]
png(paste0(k,"_corrplot.png"), width = 2, height = 2, res = 200,units="in")
corrplot(as.matrix(res),
         is.corr = F,
         outline = T,
         col = rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#FFFFFF"))(100)),
         tl.col = "black",
         cl.pos="n",
         na.label = " ")
dev.off()
png(paste0(k,"_corrplot_legend.png"), width = 5, height = 5, res = 200,units="in")
corrplot(as.matrix(res),
         is.corr = F,
         outline = T,
         col = rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#FFFFFF"))(100)),
         tl.col = "black",
         na.label = " ")
dev.off()

#scatterplot
obj2@meta.data$group_by<-paste0(obj2@meta.data$condition,"-",obj2@meta.data$my.clusters2)
Idents(obj2)<-"group_by"
df<-AverageExpression(obj2)$RNA
df<-as.data.frame(df)
colnames(df)<-gsub("-","_",colnames(df))
df$diff2<-"NS"
df$diff2[df$Mus_FDC-df$Hs_FDC > 5]<-"UP"
df$diff2[df$Mus_FDC-df$Hs_FDC < -5]<-"DOWN"
df$gene<-rownames(df)
cor.test(df$Mus_FDC,df$Hs_FDC)
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
FeaturePlot(object = obj2, features = "Pdcd1", cols = c("grey", "blue"), 
            split.by="condition",pt.size=0.001, order=T)
ggsave2("umap.cd274.organ2.png",width=20, height=3,device="png")

# Cluster frequency comparison
ggplot(obj2@meta.data, aes(x = my.clusters2,fill = condition)) +  
  geom_bar(aes(y = (..count..)/sum(..count..)),position='dodge')+ 
  scale_y_continuous(labels = percent)+ 
  labs(y="Frequency", x = "Cluster")
ggsave2("clusterfreq.png",width=4, height=2.5,device="png")

#bycondition
cluster.Counts <- table(obj2@meta.data$my.clusters2,obj2@meta.data$condition)
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
    mutate(up_in = ifelse(p_val_adj >= thresh_p_val_adj, "NS", ifelse(avg_log2FC > thresh_lfc, "FDC", ifelse(avg_log2FC < -thresh_lfc, "TRC", "NS"))),
           rank_pval = rank(p_val_adj), rank_lfc_inc = rank(avg_log2FC), rank_lfc_dec = rank(-abs(avg_log2FC)),
           lp = -log10(p_val_adj)) %>% 
    arrange(-abs(avg_log2FC))
  if (sum(plt_df$up_in != "NS") >= 30) {
    plt_df <- plt_df %>% mutate(gene_label = ifelse(up_in != "NS" & (rank_pval < 10 | rank_lfc_inc < 10 | rank_lfc_dec < 10), gene, NA))
  } else {
    plt_df <- plt_df %>% mutate(gene_label = ifelse((up_in != "NS" | rank_pval < 10 | rank_lfc_inc < 10 | rank_lfc_dec < 10), gene, NA))
  }
  plt_df$gene_label[plt_df$gene=="Pdcd1"]<-"Pdcd1"
  plt_df$lp <- pmin(plt_df$lp, 300)
  plt_df$avg_log2FC <- pmin(plt_df$avg_log2FC, 3)
  ggplot(plt_df, aes(x = avg_log2FC, y = lp, color = up_in, label = gene_label)) +
    geom_point(size = 0.5) +
    geom_text_repel(color="black",size = 3, segment.size = 0.2,seed = 1, show.legend=F, box.padding=0.5) +
    geom_hline(yintercept=10,linetype=2,size=0.2)+
    geom_vline(xintercept=0.3,linetype=2,size=0.2)+
    geom_vline(xintercept= -0.3,linetype=2,size=0.2)+
    scale_color_manual(values = c("FDC" = viridis_pal(option="plasma")(8)[1], "TRC" = viridis_pal(option="plasma")(8)[6], "NS" = "grey80")) +
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
j="tfr_pd1.DE"
res<-get(j)
thresh_p_val_adj <- 1e-10
thresh_lfc <-0.3
plt_df<- res %>% rownames_to_column(var = "gene") %>% 
   mutate(up_in = ifelse(p_val_adj >= thresh_p_val_adj, "NS", ifelse(avg_log2FC > thresh_lfc, "PD-L1 high", ifelse(avg_log2FC < -thresh_lfc, "PD-L1 neg", "NS"))),
         rank_pval = rank(p_val_adj), rank_lfc_inc = rank(avg_log2FC), rank_lfc_dec = rank(-abs(avg_log2FC)),
         lp = -log10(p_val_adj)) %>% 
  arrange(-abs(avg_log2FC))
table(plt_df$up_in)
head(plt_df[plt_df$up_in=="PD-L1 neg",],n=50)
gene.list<-c("Pdcd1","Hif1a","Lamp1","Ctsl","Gna13","Il4","Dnase1l3","Cd40lg","Bcl6","Tox2","Art2a",
             "Il7r","Cxcr3","Il2ra","Itgae","Pou2af1","S100a10","Serpina3g","Selplg","Anxa6","Cish","Capg")
plt_df<- plt_df %>% mutate(gene_label = ifelse(gene %in% gene.list , gene, NA))
plt_df$lp <- pmin(plt_df$lp, 150)
plt_df$avg_log2FC <- pmin(plt_df$avg_log2FC, 2)
ggplot(plt_df, aes(x = avg_log2FC, y = lp, color = up_in, label = gene_label)) +
  geom_point(size = 0.5) +
  geom_text_repel(color="black",size = 3, segment.size = 0.2,seed = 1, show.legend=F, box.padding=0.5) +
  geom_hline(yintercept=10,linetype=2,size=0.2)+
  geom_vline(xintercept=0.3,linetype=2,size=0.2)+
  geom_vline(xintercept= -0.3,linetype=2,size=0.2)+
  scale_color_manual(values = c("PD-L1 high" = "#0D0887FF", "PD-L1 neg" = "#F48849FF", "NS" = "grey80")) +
  scale_x_continuous(limits = c(-max(plt_df$avg_log2FC), max(plt_df$avg_log2FC)), expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(limits = c(0, max(plt_df$lp)), expand = expansion(mult = c(0.05, 0.01))) +
  labs(x=NULL,y=NULL) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),legend.position="none",#legend.title=element_text(size=8),
        plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 10, color = "black"))
ggsave2(paste0(j,"_volcano2.png"),width=3, height=2.5,device="png")


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
ggsave2("FCvsFC.png",width=2.5, height=2.5,device="png")


#custom volcano GSEA
genes.meta<-getBM(attributes=c("ensembl_gene_id", "mgi_symbol", "start_position", "end_position", "chromosome_name", 
                               "percentage_gene_gc_content", "external_gene_name", "gene_biotype","go_id","name_1006"),filters=
                    "mgi_symbol",values=list(rownames(obj@assays[["RNA"]]@meta.features)),
                  mart=useMart("ensembl", dataset = "mmusculus_gene_ensembl",host="ensembl.org"),useCache=F) #useast.

m_df = msigdbr(species = "Mus musculus")#, category = "C7") #H = hallmarks, C2=curated, C5=GO, C7=immune
pathways<-m_df %>% split(x = .$gene_symbol, f = .$gs_name)
msig.df<-m_df %>% dplyr::select(gs_name,entrez_gene)
m_df[m_df$gs_id=="M15399",]

for(i in c("M15399","M15816","M16940","M24299",
           "M2810","M14645","M40427","M13926","M16933",
           "GO:0072676","GO:0006090","GO:1903037","GO:0072678","GO:0006096","GO:0042110",
           "GO:0007159","GO:0040017","GO:0050852","GO:2000404","GO:0061615")){
  if(grepl("GO:",i)){gene.list2<-genes.meta$mgi_symbol[genes.meta$go_id==i]}else{
    gene.list2<-unique(m_df$gene_symbol[m_df$gs_id==i]) #M18810
  }
  if(length(gene.list2)>5){
    j="tfr_pd1.DE"
    res<-get(j)
    thresh_p_val_adj <- 1e-10
    thresh_lfc <-0.3
    plt_df<- res %>% rownames_to_column(var = "gene") %>% 
      mutate(up_in = ifelse(p_val_adj >= thresh_p_val_adj, "NS", ifelse(avg_log2FC > thresh_lfc, "FDC", ifelse(avg_log2FC < -thresh_lfc, "TRC", "NS"))),
             rank_pval = rank(p_val_adj), rank_lfc_inc = rank(avg_log2FC), rank_lfc_dec = rank(-abs(avg_log2FC)),
             lp = -log10(p_val_adj)) %>% 
      arrange(-abs(avg_log2FC))
    plt_df <- plt_df %>% mutate(gene_label = ifelse(up_in != "NS" & gene %in% gene.list2 &
                                                      (rank_pval < 100 | rank_lfc_inc < 100 | rank_lfc_dec < 100), gene, NA))
    plt_df <- plt_df %>% mutate(label = ifelse(gene %in% gene.list2, "sig", "notsig"))
    plt_df<-plt_df[order(plt_df$label),]
    plt_df$lp <- pmin(plt_df$lp, 150)
    plt_df$avg_log2FC <- pmin(plt_df$avg_log2FC, 3)
    ggplot(plt_df, aes(x = avg_log2FC, y = lp, color = label, label = gene_label,size=label)) +
      geom_point() +
      geom_text_repel(color="black",size = 2, segment.size = 0.2,seed = 1, show.legend=T, box.padding=0.5) +
      geom_hline(yintercept=10,linetype=2,size=0.2)+
      geom_vline(xintercept=0.3,linetype=2,size=0.2)+
      geom_vline(xintercept= -0.3,linetype=2,size=0.2)+
      scale_size_manual(values=c(0.5,1))+
      scale_color_manual(values = c("sig" = "black", "notsig" = "grey80")) +
      scale_x_continuous(limits = c(-max(plt_df$avg_log2FC), max(plt_df$avg_log2FC)), expand = expansion(mult = c(0.01, 0.01))) +
      scale_y_continuous(limits = c(0, max(plt_df$lp)), expand = expansion(mult = c(0.05, 0.01))) +
      labs(x=NULL,y=NULL) + 
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            text=element_text(size=5),
            panel.background = element_blank(),legend.position="none",#legend.title=element_text(size=8),
            plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
            axis.title = element_text(size = 10, color = "black"))
    ggsave2(paste0(j,"_volcano_gsea",i,".png"),width=1.7, height=1.5,device="png")
    
    gene.list<-df$gene[df$gene %in% gene.list2 & (df$rank_diff<100 | df$rank_sum <100 |df$rank_sum_dec<100)]
    df <- df %>% mutate(label = ifelse(gene %in% gene.list2, "sig", "notsig")) #genes to color
    df<-df[order(df$label),]
    df<- df %>% mutate(gene_label = ifelse(gene %in% gene.list , gene, NA))
    ggplot(df, aes(b6FC,m564FC,colour=label,size=label,label=gene_label)) + 
      geom_point()+
      theme_classic()+
      scale_size_manual(values=c(0.5,1))+
      geom_text_repel(color="black",size = 2, segment.size = 0.2,seed = 1, show.legend=T, box.padding=0.5) +
      theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position="none",
            axis.ticks=element_blank(),axis.line=element_blank())+ #axis.text=element_blank(),
      geom_hline(yintercept=0,color=alpha("black",0.5))+geom_vline(xintercept=0,color=alpha("black",0.5))+#geom_abline(intercept = 0, slope = 1)+
      scale_color_manual(values = c("sig" ="black", "notsig" = "grey70"))
    ggsave2(paste0(i,"_FCvsFC.png"),width=3, height=3,device="png")
    
    obj3 <- AddModuleScore(object = obj2,features = list(gene.list2),name = 'temp')
    FeaturePlot(obj3, "temp1") +
      NoLegend()+ NoAxes()+theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
                                 plot.title=element_blank()) & #NoLegend()+
      scale_color_gradientn(colors = viridis(n = 10, direction = -1))
    ggsave2(paste0(i,"_umap.png"),width=3.5, height=3.5,device="png")
    # ggsave2(paste0(k,"_pd1_umap.png"),width=5, height=4.5,device="png")
  }
}

FeaturePlot(obj3, "temp1") & scale_color_gradientn(colors = viridis(n = 10, direction = -1))
ggsave2("enrichment_legend.png",width=3.5, height=3.5,device="png")

