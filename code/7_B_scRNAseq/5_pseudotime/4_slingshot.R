
rm(list=ls())
library(gam)
library(slingshot)
library(Seurat)
library(mclust)
library(RColorBrewer)
library(clusterExperiment)
library(corrplot)
library(ggplot2)
library(ggbeeswarm)
library(ggthemes)
library(scales)
library(cowplot)
library(viridis)
library(plyr)
library(dplyr)
library(BUSpaRse)
library(tidyverse)
library(tidymodels)
library(Matrix)
library(SummarizedExperiment)
library(gam)
library(scater)


#' Assign a color to each cell based on some value
#' @param cell_vars Vector indicating the value of a variable associated with cells.
#' @param pal_fun Palette function that returns a vector of hex colors, whose
#' argument is the length of such a vector.
#' @param ... Extra arguments for pal_fun.
#' @return A vector of hex colors with one entry for each cell.
cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}
cell_colors_clust <- cell_pal(obj$my.clusters2, hue_pal())
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)

# Slingshot on Seurat clusters
obj.sce$slingPseudotime_1 <- NULL  # remove old slingshot pseudotime data
colData(obj.sce)$my.clusters2 <- as.character(obj@meta.data$my.clusters2)  # go from factor to character
table(obj@meta.data$my.clusters2)
obj.sce <- slingshot(obj.sce, clusterLabels = 'my.clusters2', reducedDim = 'UMAP',start.clus="Naive",end.clus=c("MBC_1","ASC_1","PC","MBC_2","ASC_2","ASC_3"))
#UMAP
png("sling.umap.pseudo.png",width=12, height=12,units="in",res=300)
plot(reducedDims(obj.sce)$UMAP, col = colors[cut(obj.sce$slingPseudotime_1,breaks=100)], pch=16, asp = 1)
lines(SlingshotDataSet(obj.sce), lwd = 2, type = 'lineages', col = 'black')
dev.off()
png("sling.umap.pseudo.cluster.png",width=12, height=12,units="in",res=300)
plot(reducedDims(obj.sce)$UMAP,  col = cell_colors_clust,pch=16, asp = 1)
lines(SlingshotDataSet(obj.sce), lwd = 2, type = 'lineages', col = 'black')
dev.off()
png("sling.umap.pseudo.curve.png",width=12, height=12,units="in",res=300)
plot(reducedDims(obj.sce)$UMAP, col = colors[cut(obj.sce$slingPseudotime_1,breaks=100)], pch=16, asp = 1)
lines(SlingshotDataSet(obj.sce), lwd = 2)
dev.off()
png("sling.umap.pseudo.curve.cluster.png",width=12, height=12,units="in",res=300)
plot(reducedDims(obj.sce)$UMAP,  col = cell_colors_clust,pch=16, asp = 1)
lines(SlingshotDataSet(obj.sce), lwd = 2)
dev.off()

# Plot Slingshot pseudotime vs cell stage. 
df<-colData(obj.sce)
ggplot(as.data.frame(colData(obj.sce)[,c("my.clusters2","slingPseudotime_1")]), aes(x = slingPseudotime_1, y = my.clusters2, colour = my.clusters2)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime")
ggsave2("sling.cluster.pseudo.state.png",width=4, height=4,device="png")
save(obj.sce, file = "slingshot.path.RData")

#plot
Idents(obj) <- "my.clusters2"
obj@meta.data$slingshot.pseudotime<-obj.sce$slingPseudotime_1
obj@meta.data$slingshot.pseudo.rank <- rank(obj@meta.data$slingshot.pseudotime) 
FeaturePlot(obj, features= "slingshot.pseudotime", cols= viridis(100, begin = 0))+labs(color="Pseudotime")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),plot.title = element_blank(),
        axis.ticks=element_blank(),axis.line=element_blank(),axis.text=element_blank())
ggsave2("umap.slingshot.pseudo.png",width=5.5, height=4,device="png")
VlnPlot(obj, features = "slingshot.pseudotime",pt.size=0)+ NoLegend()+labs(y="Pseudotime")+
  theme(axis.title.x=element_blank(),plot.title=element_blank())
ggsave2("cluster.by.slingshotpseudo.png",width=4, height=4,device="png")
VlnPlot(obj, features = "slingshot.pseudo.rank",pt.size=0)+ NoLegend()+labs(y="Pseudotime Rank")+
  theme(axis.title.x=element_blank(),plot.title=element_blank())
ggsave2("cluster.by.slingshotpseudo.rank.png",width=4, height=4,device="png")

#temporally expressed genes heatmap
t <- obj.sce$slingPseudotime_1
Y <- log1p(assays(obj.sce)$counts)
var100 <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:100]
Y <- Y[var100,]
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  suppressWarnings({
    tmp <- suppressWarnings(gam(z ~ lo(t), data=d))
  })
  p <- summary(tmp)[3][[1]][2,3]
  p
})
topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]
heatdata <- assays(obj.sce)$counts[topgenes, order(t, na.last = NA)]
heatclus <- obj.sce$my.clusters2[order(t, na.last = NA)]
png("sling.pseudo.heatmap.png",width=12, height=12,units="in",res=300)
heatmap(as.matrix(log1p(heatdata)), Colv = NA,ColSideColors = brewer.pal(9,"Set1")[heatclus])
dev.off()

#expression over time
png("sling.pseudo.CD274.overtime.png",width=4, height=4,units="in",res=300)
plotExpression(obj.sce, "CD274", x = "slingPseudotime_1", 
               colour_by = "my.clusters2", show_violin = FALSE,
               show_smooth = TRUE)
dev.off()

#lineage structure
rd<-reducedDims(obj.sce)$UMAP
cl<-obj.sce$my.clusters2
lin1 <- getLineages(rd, cl, start.clus = 'Naive')
png("slingshot.umap.cluster.lineage.png",width=8,height=9, units="in",res=300)
plot(rd,  asp = 1, pch = 16,col=cell_colors_clust) #col = brewer.pal(9,"Set1")[cl],
lines(SlingshotDataSet(lin1), lwd = 3, col = 'black')
dev.off()
lin2 <- getLineages(rd, cl, start.clus= 'Naive', end.clus = 'PC')
png("slingshot.umap.cluster.lineage2.png",width=8,height=9, units="in",res=300)
plot(rd, col = cell_colors_clust, asp = 1, pch = 16)
lines(SlingshotDataSet(lin2), lwd = 3, col = 'black', show.constraints = TRUE)
dev.off()
crv1 <- getCurves(lin1)
crv1
png("slingshot.umap.cluster.curves.png",width=8,height=9, units="in",res=300)
plot(rd, col = cell_colors_clust, asp = 1, pch = 16)
lines(SlingshotDataSet(crv1), lwd = 3, col = 'black')
dev.off()
save.image("temp2.slingshot.RData")

#DE
# Get top highly variable genes
DefaultAssay(obj)<-"RNA"
obj<- SCTransform(obj, vars.to.regress = c("nCount_RNA"))
head(obj[["RNA"]]@meta.features[obj[["RNA"]]@meta.features$vst.variable,])
dimnames(obj[["RNA"]]@meta.features)
test<-HVFInfo(obj[["RNA"]], method = 'vst')
top_hvg <- HVFInfo(obj[["RNA"]], method = 'vst') %>% 
  mutate(., bc = rownames(.)) %>% 
  arrange(desc(variance)) %>% 
  top_n(300, variance) %>% 
  pull(bc)
save.image("temp.slingshot.DE.RData")

# Prepare data for random forest
dat_use <- t(GetAssayData(obj, slot = "data")[top_hvg,])
dat_use_df <- cbind(slingPseudotime(obj.sce)[,2], dat_use) # Do curve 2, so 2nd columnn
colnames(dat_use_df)[1] <- "pseudotime"
dat_use_df <- as.data.frame(dat_use_df[!is.na(dat_use_df[,1]),])
dat_split <- initial_split(dat_use_df)
dat_train <- training(dat_split)
dat_val <- testing(dat_split)
model <- rand_forest(mtry = 200, trees = 1400, min_n = 15, mode = "regression") %>%
  set_engine("ranger", importance = "impurity", num.threads = 3) %>%
  fit(pseudotime ~ ., data = dat_train)
val_results <- dat_val %>% 
  mutate(estimate = predict(model, .[,-1]) %>% pull()) %>% 
  select(truth = pseudotime, estimate)
metrics(data = val_results, truth, estimate)
summary(dat_use_df$pseudotime)
var_imp <- sort(model$fit$variable.importance, decreasing = TRUE)
top_genes <- names(var_imp)[1:9]
png("sling.clusterbypseudo.png",width=7, height=8,units="in",res=300)
par(mfrow = c(3, 3))
pal <- viridis(100, end = 0.95)
for (i in seq_along(top_genes)) {
  colors <- pal[cut(dat_use[,top_genes[i]], breaks = 100)]
  plot(reducedDims(obj.sce)$UMAP, col = colors, 
       pch = 16, cex = 0.5, main = top_genes[i])
  lines(SlingshotDataSet(obj.sce), lwd = 2, col = 'black', type = 'lineages')
}
dev.off()
save.image("temp.slingshot.DE2.RData")


