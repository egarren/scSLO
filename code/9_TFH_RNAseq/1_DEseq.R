rm(list=ls())
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(DEGreport)
library(tximport)
library(ggplot2)
library(ggrepel)
library(AnnotationHub)
library(ensembldb)
library(purrr)
library(dplyr)
library(stringr)
library(pheatmap)
library(EnhancedVolcano)
library(cowplot)
library(vsn)
library(Seurat)
save_pheatmap_png <- function(x, filename, width=1000, height=1400, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
library(DESeq2)


#transcript annotation
ah = AnnotationHub()
ms_ens <- query(ah, c("Mus musculus", "EnsDb"))
ms_ens <- ms_ens[["AH89211"]]
txdb <- transcripts(ms_ens, return.type = "data.frame") %>%dplyr::select(tx_id, gene_id)
txdb <- txdb[grep("ENSMUST", txdb$tx_id),]
genedb <- genes(ms_ens, return.type = "data.frame")  %>%dplyr::select(gene_id, symbol)
tx2gene <- inner_join(txdb, genedb)

## List all directories containing data  
samples <- list.files(path = "./counts", full.names = T, pattern="_out")
files <- file.path(samples, "quant.sf")
names(files) <- str_replace(samples, "./counts/LIB061786_", "") %>% str_replace("_out", "")

# Run tximport
txi <- tximport(files, type="salmon", tx2gene=tx2gene[,c("tx_id", "gene_id")], countsFromAbundance="lengthScaledTPM",ignoreTxVersion = TRUE)
attributes(txi)
data <- txi$counts %>% round() %>% data.frame()
write.csv(data,"count_mtx.csv")

#process matrices
meta<-data.frame(sample=str_split_i(colnames(txi$counts), "_", 2))
meta<-cbind(meta,str_split_fixed(meta$sample,"-",2))
colnames(meta)<-c("sample","mouse","cell.type")
meta$mouse<-as.numeric(meta$mouse)
meta2<-read.csv("meta.csv", header=T) #define path
meta2$genotype<-factor(meta2$genotype)
meta<-left_join(meta,meta2,by="mouse")
meta$genotype_celltype<-paste0(meta$genotype,"_",meta$cell.type)
rownames(meta)<-colnames(txi$counts)
write.csv(meta,"meta2.csv")

#explore multiple genes
goi<-c("Sostdc1", "Pdcd1", "Zfp703", "Bcl6", "Tbc1d4", "Maf", "Hif1a", "Tox2", "Tnfsf8", "Nav2", "Eea1",
       "Ppp1r14b", "Art2a", "Cxx5", "Ctsl", "Gna13", "Apoe","Aldoa", "Cebpa", "Egln3", "Gfra4", "Dap", "Gpm6b",
       "Spp1", "Fam20a","Mki67",
       "Il4", "Cd40lg", "Lrp5", "Anxa6", "Il7r", "Rasgrp2", "Ccr7", "Sell", "Slamf7", "Nkg7", "Eomes", 
       "Arhgef18",  "Spp1", "Selplg", "Cd74","Tnfrsf9", "Itgb7")
data2<-data[rownames(data)%in% tx2gene$gene_id[tx2gene$symbol %in% goi] ,]
data2<-as.data.frame(t(data2))
dict<-unique(tx2gene[tx2gene$symbol %in% goi,c("gene_id","symbol")]) 
colnames(data2)<-dplyr::recode(
  colnames(data2),!!!setNames(as.character(dict$symbol),dict$gene_id)
)
meta3<-meta
rownames(meta3)<-gsub("-",".",rownames(meta3))
mat<-t(data2)
keep<-rownames(meta3)[meta3$cell.type=="TFH"&!(meta3$mouse %in% c("1805","1806"))]
pheatmap(mat[,keep], 
             color = colorRampPalette(c("#45628f", "#F7F7F7", "#de425b"))(50),
             cluster_rows = T,
             show_rownames = T,
             annotation = meta3[,c("genotype")],
             border_color = NA, 
             fontsize = 10, 
             scale = "row", 
             fontsize_row = 10, 
             height = 20)

data3<-data2
colnames(data3)<-make.unique(colnames(data3))
data3$sample<-rownames(data3)
meta4<-meta
meta4$sample<-gsub("-",".",rownames(meta4))
data3<-left_join(data3,meta4,by="sample")
data3<-data3[data3$cell.type=="TFH",]
write.csv(data3[order(data3$genotype),],"gois.csv")

#explore data
ggplot(data) +
  geom_histogram(aes(x = TRA00276855_1794.TFH_S1), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")
mean_counts <- apply(data, 1, mean)        #The second argument '1' of 'apply' function indicates the function being applied to rows. Use '2' if applied to columns 
variance_counts <- apply(data, 1, var)
df <- data.frame(mean_counts, variance_counts)
ggplot(df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) + 
  scale_y_log10(limits = c(1,1e9)) +
  scale_x_log10(limits = c(1,1e9)) +
  geom_abline(intercept = 0, slope = 1, color="red")

#DEseq object
### Check that sample names match in both files
all(colnames(txi$counts) %in% rownames(meta))
all(colnames(txi$counts) == rownames(meta))
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ genotype_celltype)
dds
colData(dds)[sapply(colData(dds), is.character)] <- lapply(colData(dds)[sapply(colData(dds), is.character)], as.factor)
save.image("temp0.DEseq.RData")

#normalization
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="normalized_counts.txt", sep="\t", quote=F, col.names=NA)

#QC using transformation
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
ntd <- normTransform(dds)
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

#QC using PCA
df<-plotPCA(rld, intgroup="sample",returnData=T)
df<-left_join(df,meta,by="sample")
ggplot(df,aes(x=PC1,y=PC2,color=genotype,shape=cell.type,label=sample))+
  geom_text_repel(max.overlaps=1e10)+
  geom_point()+coord_fixed()
rld_mat <- assay(rld)    
rld_cor <- cor(rld_mat)    
pheatmap(rld_cor, annotation = meta)

#DGE
dds <- DESeq(dds)
plotDispEsts(dds) ## Plot dispersion estimates
sizeFactors(dds)
colSums(counts(dds)) #raw counts per sample
colSums(counts(dds, normalized=T)) #normalized counts per sample
design(dds) <- formula(~cell.type)
dds <- DESeq(dds)
res_cell.type <- results(dds,contrast=c("cell.type", "TFH", "TFR"))
save.image("temp.DEseq.RData")

#subset TFH and TFR
for(i in c("TFH","TFR")){
  keep<-rownames(meta[meta$cell.type==i,])
  dds2<-dds[,colnames(dds) %in% keep]
  colData(dds2)$genotype<-as.character(colData(dds2)$genotype)
  colData(dds2)[sapply(colData(dds2), is.character)] <- lapply(colData(dds2)[sapply(colData(dds2), is.character)], as.factor)
  assign(paste0("dds_",i),dds2)
  design(dds2) <- formula(~ genotype)
  dds2 <- DESeq(dds2)
  res2 <- results(dds2,contrast=c("genotype", "564homo","564homo.Cd21-PDL1"))
  assign(paste0("res_genotype_",i),res2)
}

#subset genotypes
for(i in c("564homo","564homo.Cd21-PDL1")){
  keep<-rownames(meta[meta$genotype==i,])
  dds2<-dds[,colnames(dds) %in% keep]
  colData(dds2)$cell.type<-as.character(colData(dds2)$cell.type)
  colData(dds2)[sapply(colData(dds2), is.character)] <- lapply(colData(dds2)[sapply(colData(dds2), is.character)], as.factor)
  assign(paste0("dds_",i),dds2)
  design(dds2) <- formula(~ cell.type)
  dds2 <- DESeq(dds2)
  res2 <- results(dds2,contrast=c("cell.type", "TFH", "TFR"))
  assign(paste0("res_celltype_",i),res2)
}
save.image("temp2.DEseq.RData")

#Volcano
for(i in ls(pattern="res_")){
  res_df<-get(i)
  res_df<- res_df[order(res_df$pvalue),]
  volcano<-data.frame(gene=rownames(res_df),log2FC=res_df$log2FoldChange,p_val_adj=res_df$padj)
  volcano<-left_join(volcano,unique(tx2gene[,c("gene_id","symbol")]),by=c("gene"="gene_id"),keep=F)
  volcano<-volcano[!is.na(volcano$p_val_adj),c("log2FC","p_val_adj","symbol")]
  volcano$gene<-volcano$symbol
  volcano$avg_log2FC = -volcano$log2FC

  thresh_p_val_adj <- 0.05
  thresh_lfc <-2
  plt_df<- volcano %>% 
    mutate(up_in = ifelse(p_val_adj >= thresh_p_val_adj, "NS", ifelse(avg_log2FC > thresh_lfc, "UP", ifelse(avg_log2FC < -thresh_lfc, "DOWN", "NS"))),
           rank_pval = rank(p_val_adj), rank_lfc_inc = rank(avg_log2FC), rank_lfc_dec = rank(-abs(avg_log2FC)),
           lp = -log10(p_val_adj)) %>% 
    arrange(-abs(avg_log2FC))
  table(plt_df$up_in)
  if (sum(plt_df$up_in != "NS") >= 30) {
    plt_df <- plt_df %>% mutate(gene_label = ifelse(up_in != "NS" & (rank_pval < 10 | rank_lfc_inc < 10 | rank_lfc_dec < 10), gene, NA))
  } else {
    plt_df <- plt_df %>% mutate(gene_label = ifelse((up_in != "NS" | rank_pval < 10 | rank_lfc_inc < 10 | rank_lfc_dec < 10), gene, NA))
  }
  plt_df$gene_label[plt_df$gene=="Sostdc1"]<-"Sostdc1"
  plt_df$lp <- pmin(plt_df$lp, 5)
  plt_df$avg_log2FC <- pmin(plt_df$avg_log2FC, 5)
  ggplot(plt_df, aes(x = avg_log2FC, y = lp, color = up_in, label = gene_label)) +
    geom_point(size = 0.5) +
    geom_text_repel(color="black",size = 3, segment.size = 0.2,seed = 1, show.legend=F, box.padding=0.5) +
    geom_hline(yintercept=1.3,linetype=2,size=0.2)+
    geom_vline(xintercept=2,linetype=2,size=0.2)+
    geom_vline(xintercept= -2,linetype=2,size=0.2)+
     scale_color_manual(values = c("UP" = "darkgreen", "DOWN" = "#D95F02", "NS" = "grey80")) +
    scale_x_continuous(limits = c(-max(plt_df$avg_log2FC), max(plt_df$avg_log2FC)), expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(limits = c(0, 5), expand = expansion(mult = c(0.05, 0.01))) +
     labs(x=NULL,y=NULL) + 
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),legend.position="none",#legend.title=element_text(size=8),
          plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
          axis.title = element_text(size = 10, color = "black"))
  ggsave2(paste0("output/",i,"_volcano.png"),width=3, height=2.5,device="png")
  write.csv(volcano,paste0("output/",i,".csv") )
}

###Visualization
tibble_meta <- meta %>% rownames_to_column(var="samplename") %>% as_tibble()
normalized_counts <- counts(dds, normalized=T) %>% data.frame() %>%rownames_to_column(var="gene") 
mm10annot <- tx2gene %>% dplyr::select(gene_id, symbol) %>% dplyr::distinct()
normalized_counts <- merge(normalized_counts, mm10annot, by.x="gene", by.y="gene_id")
normalized_counts <- normalized_counts %>%as_tibble()

# plot counts
for(i in c("Sostdc1","Pdcd1","Cd40lg")){
  for(j in c("cell.type","genotype","genotype_celltype")){
    d <- plotCounts(dds, gene=mm10annot[mm10annot$symbol == i, "gene_id"], intgroup=j, returnData=TRUE)
    ggplot(d, aes_string(x = j, y = "count", color = j)) + 
      geom_point(position=position_jitter(w = 0.1,h = 0)) +
      geom_text_repel(aes(label = rownames(d))) + 
      theme_bw() +
      ggtitle(i) +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave2(paste0("output/",i,"_",j,"_counts.png"),width=10, height=7,device="png")
  }
}

#heatmap
meta3<-meta
for(i in ls(pattern="res_")){
  res_df<-get(i)
  res_df<- res_df[order(res_df$pvalue),]
  volcano<-data.frame(gene=rownames(res_df),log2FC=res_df$log2FoldChange,p_val_adj=res_df$padj)
  volcano<-left_join(volcano,unique(tx2gene[,c("gene_id","symbol")]),by=c("gene"="gene_id"),keep=F)
  volcano<-volcano[!is.na(volcano$p_val_adj),]
  select_genes<- head(volcano,n=20)
  mat <- normalized_counts %>% dplyr::filter(gene %in% select_genes$gene)  
  mat<-as.data.frame(mat)
  rownames(mat)<-mat$symbol
  rownames(meta3)<-gsub("-",".",rownames(meta3))
  select_sample<-rownames(meta3)
  if(i=="res_celltype_564homo.Cd21-PDL1"){select_sample<-rownames(meta3)[meta3$genotype=="564homo.Cd21-PDL1"]}
  if(i=="res_celltype_564homo"){select_sample<-rownames(meta3)[meta3$genotype=="564homo"]}
  if(i=="res_genotype_TFH"){select_sample<-rownames(meta3)[meta3$cell.type=="TFH"]}
  if(i=="res_genotype_TFR"){select_sample<-rownames(meta3)[meta3$cell.type=="TFR"]}
  p1<-pheatmap(mat[,select_sample], 
               color = colorRampPalette(c("#45628f", "#F7F7F7", "#de425b"))(50),
               cluster_rows = T, 
               show_rownames = T,
               annotation = meta[,c("genotype","cell.type"),drop=F], 
               border_color = NA, 
               fontsize = 10, 
               scale = "row", 
               fontsize_row = 10, 
               height = 20)
  save_pheatmap_png(p1, paste0("output/",i,".heatmap.png"))
}

#distance matrix
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Genotype, vsd$Cell.type, sep="-")
rownames(sampleDistMatrix) <- paste(vsd$Genotype)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
p1<-pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
save_pheatmap_png(p1, paste0("output/sampledist.heatmap.png"))

#PCA
pcaData <- plotPCA(vsd, intgroup=c("genotype_celltype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
meta$name<-rownames(meta)
pcaData<-left_join(pcaData,meta[,c("genotype","cell.type","name")],by="name")
ggplot(pcaData, aes(PC1, PC2,  color=cell.type,shape=genotype,label=name)) +
  geom_point(size=3) +#+coord_fixed()+#geom_text()+
  scale_color_manual(values = c("TFH" = "#E97132", "TFR" = "#7030A0")) +
  scale_shape_manual(values = c(1,16)) +
  xlab( paste0("PC1 (",percentVar[1],"%)") ) +
  ylab( paste0("PC2 (",percentVar[2],"%)") ) +
  guides(color=guide_legend(title=""),shape=guide_legend(title=""))+
  theme(
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))
ggsave2("output/pca.png",width=4.5, height=2.5,device="png")

##heatmap
goi<-c("Sostdc1","Foxp3")
for(i in ls(pattern="res_")){
  res_df<-get(i)
  res_df<- res_df[order(res_df$pvalue),]
  volcano<-data.frame(gene=rownames(res_df),log2FC=res_df$log2FoldChange,p_val_adj=res_df$padj)
  volcano<-left_join(volcano,unique(tx2gene[,c("gene_id","symbol")]),by=c("gene"="gene_id"),keep=F)
  volcano<-volcano[!is.na(volcano$p_val_adj),]
  volcano<-volcano[volcano$p_val_adj<0.01,]
  goi<-union(goi,volcano$gene)
}
mat <- normalized_counts %>% dplyr::filter(gene %in% goi)  
mat<-as.data.frame(mat)
rownames(mat)<-mat$symbol
rownames(meta3)<-gsub("-",".",rownames(meta3))
select_sample<-rownames(meta3)
p1<-pheatmap(mat[,select_sample], 
             color = colorRampPalette(c("#45628f", "#F7F7F7", "#de425b"))(50),
             cluster_rows = T, 
             show_rownames = F,
             show_colnames = F,
             treeheight_row=0,
             treeheight_col=10,
             annotation = meta3[,c("genotype","cell.type"),drop=F], #"Cell.type",
             annotation_colors=list(cell.type=c("TFH" = "#E97132", "TFR" = "#7030A0"),
                                    genotype=c("564homo.Cd21-PDL1" = "#00B050", "564homo" = "#0070C0")),
             annotation_names_col=F,
             cutree_rows=3,
             border_color = NA, 
             fontsize = 10, 
             scale = "row", 
             fontsize_row = 10, 
             height = 20)
save_pheatmap_png(p1, "output/pub_heatmap.png")
mat2<-mat[,select_sample]
mat2<-mat2[p1$tree_row$order,p1$tree_col$order]
mat2$gene<-rownames(mat2)
cl<-cutree(p1$tree_row,3)
cl<-data.frame(cl)
cl$gene<-rownames(cl)
mat3<-left_join(mat2,cl,by="gene")
write.csv(mat3,"output/res_heatmap.csv")

#FC vs FC
x<-get("res_celltype_564homo")
y<-get("res_celltype_564homo.Cd21-PDL1")
x<-data.frame(gene=rownames(x),log2FC= -x$log2FoldChange,p_val_adj=x$padj)
y<-data.frame(gene=rownames(y),log2FC= -y$log2FoldChange,p_val_adj=y$padj)
df<-left_join(x[,c("log2FC","gene")],y[,c("log2FC","gene")],by="gene",keep=F)
df<-left_join(df,unique(tx2gene[,c("gene_id","symbol")]),by=c("gene"="gene_id"),keep=F)
rownames(df)<-make.names(df$symbol,unique=T)
df$diff<-abs(df$log2FC.x - df$log2FC.y)
sig<-df$symbol[df$diff>5]
df <- df %>% mutate(label = ifelse(symbol %in% sig, "sig", "notsig"))
df$label[df$log2FC.x>1 & df$log2FC.y>1]<-"UP"
df$label[df$log2FC.x< -1 & df$log2FC.y< -1]<-"DOWN"
df$label[df$log2FC.x>1 & df$log2FC.y< -1]<-"X"
df$label[df$log2FC.x< -1 & df$log2FC.y>1]<-"Y"
write.csv(df,"output/FC_FC.csv")
plt_df<-df
border=10
plt_df$log2FC.x<-pmin(plt_df$log2FC.x,border)
plt_df$log2FC.x<-pmax(plt_df$log2FC.x,-border)
plt_df$log2FC.y<-pmin(plt_df$log2FC.y,border)
plt_df$log2FC.y<-pmax(plt_df$log2FC.y,-border)
gene.list<-c("Foxp3","Cd40lg","Ccr8","Gzmk","Tnfrsf18",
             "Izkf4","Slamf6","Eomes","Capg",
             "Tnfsf11","Itgb1","Tnfsf8","Aldh1a2","Cilp2","Frmpd4","Areg","Ccn2","Hgf",
             "Tacr1","Pet117","Pabpc2","Mei4","Apcdd1",
             "Tnfrsf11a","Arhgap36","Smoc2","Klf15",
             "Zfp935","Elmod1","Frmpd4","Spata19","Dpysl5",
             "Aldh1a2","Fgf16")
plt_df <- plt_df %>% mutate(gene_label = ifelse((symbol %in% gene.list), symbol, NA))
ggplot(plt_df, aes(log2FC.x,log2FC.y,label=gene_label)) + geom_point(aes(colour=label),size=0.5)+#geom_point(fill=NA,colour=alpha("black",0.5),pch=21,size=3) +
  theme_classic()+coord_fixed()+
  labs(x="",y="")+
  scale_x_continuous(limits = c(-border,border), expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(limits = c(-border,border), expand = expansion(mult = c(0.01, 0.01))) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position="none",
        axis.ticks=element_blank(),axis.line=element_blank())+ #axis.text=element_blank(),
  geom_hline(yintercept=0,color=alpha("black",0.5))+geom_vline(xintercept=0,color=alpha("black",0.5))+
  geom_hline(yintercept=1,linetype=2,size=0.2)+
  geom_hline(yintercept= -1,linetype=2,size=0.2)+
  geom_vline(xintercept=1,linetype=2,size=0.2)+
  geom_vline(xintercept= -1,linetype=2,size=0.2)+
  geom_abline(intercept = 0, slope = 1,color=alpha("black",0.5))+
  geom_label_repel(color="black",size = 3, segment.size = 0.2,seed = 1, show.legend=F, box.padding=0.5)+
  scale_color_manual(values = c("sig" = "red", "UP" = "#E97132", "DOWN" = "#7030A0","notsig" = "grey70",
                                "Y" = "#00B050", "X" = "#0070C0"))
ggsave2("output/FCvsFC.png",width=4, height=4,device="png")




