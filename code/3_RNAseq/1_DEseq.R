## Setup
### Bioconductor and CRAN libraries used
## see https://github.com/hbctraining/DGE_workshop_salmon_online/tree/master/lessons
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
## Obtain a vector of all filenames including the path
files <- file.path(samples, "quant.sf")
## Since all quant files have the same name it is useful to have names for each element
names(files) <- str_replace(samples, "./counts/LIB056658_", "") %>% str_replace("_out", "")

# Run tximport
txi <- tximport(files, type="salmon", tx2gene=tx2gene[,c("tx_id", "gene_id")], countsFromAbundance="lengthScaledTPM",ignoreTxVersion = TRUE)
attributes(txi)
# Look at the counts
txi$counts %>% View()
# Write the counts to an object
data <- txi$counts %>% round() %>% data.frame()

#process matrices
meta<-read.csv("../meta4.csv", header=T) #define path
rownames(meta)<-colnames(txi$counts)
meta$Genotype<-factor(meta$Genotype)
meta$Cell.type<-factor(meta$Cell.type)
meta$PDL1<-factor(meta$PDL1)

#explore data
ggplot(data) +
  geom_histogram(aes(x = TRA00242316_S1), stat = "bin", bins = 200) +
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
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ Genotype)
dds

# #subset
remove<-c("TRA00242316_S1")
dds<-dds[,!colnames(dds) %in% remove]
meta<-meta[!rownames(meta)%in% remove,]

#filter
dds_all<-dds
keep1<-rowSums(counts(dds) >= 5) >= 9
keep4<-rownames(dds)%in% tx2gene$gene_id[tx2gene$symbol %in% c("Cd274")] #"Cr2",
dds<-dds[keep1 | keep4,]

#QC using transformation
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
ntd <- normTransform(dds)
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

#QC using PCA
rld <- rlog(dds, blind=TRUE)
plotPCA(rld, intgroup="Genotype")
rld_mat <- assay(rld)    
rld_cor <- cor(rld_mat)    
pheatmap(rld_cor, annotation = meta)

#DGE
dds <- DESeq(dds)
plotDispEsts(dds) ## Plot dispersion estimates
sizeFactors(dds)
colSums(counts(dds)) #raw counts per sample
colSums(counts(dds, normalized=T)) #normalized counts per sample
save.image("temp.DEseq.RData")

#multifactor analysis
design(dds) <- formula(~ Genotype)
dds <- DESeq(dds)
res_wt <- results(dds,contrast=c("Genotype","Cxcl13.PDL1_SRBC", "WT_SRBC"))
res_564 <- results(dds,contrast=c("Genotype", "564homo.Cxcl13.PDL1","564homo"))
res_564_wt <- results(dds,contrast=c("Genotype", "564homo","WT_SRBC"))
design(dds) <- formula(~ Cell.type)
dds <- DESeq(dds)
res_fdc <- results(dds,contrast=c("Cell.type", "FDC", "FRC"))
design(dds) <- formula(~ PDL1)
dds <- DESeq(dds)
res_pdl1 <- results(dds,contrast=c("PDL1", "no","yes"))

#Volcano
for(i in c("res_564","res_wt","res_fdc","res_pdl1","res_564_wt")){
  res_df<-get(i)
  res_df<- res_df[order(res_df$pvalue),]
  volcano<-data.frame(gene=rownames(res_df),log2FC=res_df$log2FoldChange,p_val_adj=res_df$padj)
  volcano<-left_join(volcano,unique(tx2gene[,c("gene_id","symbol")]),by=c("gene"="gene_id"),keep=F)
  volcano<-volcano[!is.na(volcano$p_val_adj),c("log2FC","p_val_adj","symbol")]
  volcano$gene<-volcano$symbol
  volcano$avg_log2FC = -volcano$log2FC

  thresh_p_val_adj <- 0.05
  thresh_lfc <-5
  plt_df<- volcano %>% #rownames_to_column(var = "gene") %>% 
    mutate(up_in = ifelse(p_val_adj >= thresh_p_val_adj, "NS", ifelse(avg_log2FC > thresh_lfc, "PD-L1 high", ifelse(avg_log2FC < -thresh_lfc, "PD-L1 neg", "NS"))),
           rank_pval = rank(p_val_adj), rank_lfc_inc = rank(avg_log2FC), rank_lfc_dec = rank(-abs(avg_log2FC)),
           lp = -log10(p_val_adj)) %>% 
    arrange(-abs(avg_log2FC))
  table(plt_df$up_in)
  if (sum(plt_df$up_in != "NS") >= 30) {
    plt_df <- plt_df %>% mutate(gene_label = ifelse(up_in != "NS" & (rank_pval < 10 | rank_lfc_inc < 10 | rank_lfc_dec < 10), gene, NA))
  } else {
    plt_df <- plt_df %>% mutate(gene_label = ifelse((up_in != "NS" | rank_pval < 10 | rank_lfc_inc < 10 | rank_lfc_dec < 10), gene, NA))
  }
  plt_df$gene_label[plt_df$gene=="Cd274"]<-"Cd274"
  plt_df$lp <- pmin(plt_df$lp, 5)
  plt_df$avg_log2FC <- pmin(plt_df$avg_log2FC, 15)
  ggplot(plt_df, aes(x = avg_log2FC, y = lp, color = up_in, label = gene_label)) +
    geom_point(size = 0.5) +
    geom_text_repel(color="black",size = 3, segment.size = 0.2,seed = 1, show.legend=F, box.padding=0.5) +
    geom_hline(yintercept=1.3,linetype=2,size=0.2)+
    geom_vline(xintercept=5,linetype=2,size=0.2)+
    geom_vline(xintercept= -5,linetype=2,size=0.2)+
    scale_color_manual(values = c("PD-L1 high" = "darkgreen", "PD-L1 neg" = "#D95F02", "NS" = "grey80")) +
    scale_x_continuous(limits = c(-max(plt_df$avg_log2FC), max(plt_df$avg_log2FC)), expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(limits = c(0, 10), expand = expansion(mult = c(0.05, 0.01))) +
     labs(x=NULL,y=NULL) + 
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),legend.position="none",#legend.title=element_text(size=8),
          plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
          axis.title = element_text(size = 10, color = "black"))
  ggsave2(paste0(i,"_volcano.png"),width=3, height=2.5,device="png")

  assign(paste0("DE.",i),volcano)
  write.csv(volcano,paste0(i,".csv") )
}



save.image("temp2.DEseq.RData")


###Visualization
tibble_meta <- meta %>% rownames_to_column(var="samplename") %>% as_tibble()
normalized_counts <- counts(dds, normalized=T) %>% data.frame() %>%rownames_to_column(var="gene") 
mm10annot <- tx2gene %>% dplyr::select(gene_id, symbol) %>% dplyr::distinct()
normalized_counts <- merge(normalized_counts, mm10annot, by.x="gene", by.y="gene_id")
normalized_counts <- normalized_counts %>%as_tibble()

# plot counts
for(i in c("Cd274","Cxcl13","Cr1l")){#"Cr2",
  for(j in c("Genotype","Cell.type","PDL1")){
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
for(i in c("res_564","res_wt","res_fdc","res_pdl1","res_564_wt","top_counts")){
  if(i == "top_counts"){
    select_genes<- normalized_counts[order(rowMeans(normalized_counts[,rownames(meta)]), decreasing=TRUE)[1:40],]
  }else{
    res_df<-get(i)
    res_df<- res_df[order(res_df$pvalue),]
    volcano<-data.frame(gene=rownames(res_df),log2FC=res_df$log2FoldChange,p_val_adj=res_df$padj)
    volcano<-left_join(volcano,unique(tx2gene[,c("gene_id","symbol")]),by=c("gene"="gene_id"),keep=F)
    volcano<-volcano[!is.na(volcano$p_val_adj),]
    select_genes<- head(volcano,n=20)
    }
  mat <- normalized_counts %>% dplyr::filter(gene %in% select_genes$gene)  
  mat<-as.data.frame(mat)
  rownames(mat)<-mat$symbol
  select_sample<-rownames(meta)
  if(i=="res_564"){select_sample<-rownames(meta)[meta$Genotype %in% c("564homo","564homo.Cxcl13.PDL1")]}
  if(i=="res_wt"){select_sample<-rownames(meta)[meta$Genotype %in% c("WT_SRBC","Cxcl13.PDL1_SRBC")]}
  p1<-pheatmap(mat[,select_sample], 
               color = colorRampPalette(c("#45628f", "#F7F7F7", "#de425b"))(50),
               cluster_rows = T, 
               show_rownames = T,
               annotation = meta[,c("Genotype"),drop=F], #"Cell.type",
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
pcaData <- plotPCA(vsd, intgroup=c("Genotype"), returnData=TRUE)#, "Cell.type","PDL1"
percentVar <- round(100 * attr(pcaData, "percentVar"))
meta$type<-"564homo"
meta$type[grepl("SRBC",meta$Genotype)]<-"B6"
meta$name<-rownames(meta)
pcaData<-left_join(pcaData,meta[,c("type","PDL1","name")],by="name")
ggplot(pcaData, aes(PC1, PC2,  color=PDL1,shape=type,label=name)) + #,shape=Cell.type,
  geom_point(size=3) +#geom_text()+
  scale_color_manual(values = c("yes" = "darkgreen", "no" = "#D95F02")) +
  scale_shape_manual(values = c(1,16)) +
  xlab( paste0("PC1 (",percentVar[1],"%)") ) +
  ylab( paste0("PC2 (",percentVar[2],"%)") ) +
  theme(#line = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))
ggsave2("pca.png",width=4, height=2.5,device="png")




# #likelihood ratio
dds_lrt <- DESeq(dds, test="LRT", reduced=~1)
res_LRT <- results(dds_lrt)
res_LRT_tb <- res_LRT %>%data.frame() %>%rownames_to_column(var="gene") %>% as_tibble()
sigLRT_genes <- res_LRT_tb %>% dplyr::filter(padj < 0.05)
nrow(sigLRT_genes)
clustering_sig_genes <- sigLRT_genes %>%arrange(padj) %>%head(n=1000)
cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]
clusters <- degPatterns(cluster_rlog, metadata = meta, time = "Genotype", col=NULL)
group1 <- clusters$df %>% dplyr::filter(cluster == 1)

