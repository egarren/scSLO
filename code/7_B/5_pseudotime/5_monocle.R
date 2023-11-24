

rm(list=ls())
library(monocle3)
library(ggplot2)
library(dplyr)
library(Seurat)
library(cowplot)
library(viridis)
library(scales)
library(xlsx)

# Seurat to Monocle cds
seurat.object <- obj
# data <- as(as.matrix(GetAssayData(seurat.object, assay = "integrated", slot = "scale.data")), 'sparseMatrix')
data <- as(as.matrix(GetAssayData(seurat.object, assay = "RNA", slot = "data")), 'sparseMatrix')
pd <- data.frame(seurat.object@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
#Construct monocle cds
cds <- new_cell_data_set(expression_data = data, cell_metadata = pd, gene_metadata = fData)
reducedDim(cds,"UMAP")<-seurat.object@reductions[["umap"]]@cell.embeddings
reducedDim(cds,"PCA")<-seurat.object@reductions[["pca"]]@cell.embeddings
#add parition and clusters
# cds<-cluster_cells(cds)
recreate.partition <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partition) <- cds@colData@rownames
recreate.partition <- as.factor(recreate.partition)
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition
list_cluster <- seurat.object@meta.data$my.clusters2
names(list_cluster) <- seurat.object@assays[["RNA"]]@data@Dimnames[[2]]
cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
cds@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"
plot_cells(cds,reduction_method="UMAP",color_cells_by="cluster",label_cell_groups=F,cell_size=0.3)
rm(list=setdiff(ls(), c("cds","obj")))
save.image("temp1.monocle.RData")

get_earliest_principal_node <- function(cds, time_bin="Naive"){
  cell_ids <- which(colData(cds)[, "my.clusters2"] == time_bin)
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}

#map pseudotime
cds = learn_graph(cds)
cds = order_cells(cds, reduction_method = "UMAP",root_pr_nodes=get_earliest_principal_node(cds))
for(i in c("pseudotime","cluster","partition","my.clusters2")){
  plot_cells(cds, color_cells_by = i,label_cell_groups=F,cell_size=0.3) 
  ggsave2(paste0(i,".moncole.umap.png"),width=5.5,height=4,device="png")
}

#plot
Idents(obj) <- "my.clusters2"
obj@meta.data$monocle.pseudotime<-cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
obj@meta.data$monocle.pseudo.rank <- rank(obj@meta.data$monocle.pseudotime) 
FeaturePlot(obj, features= "monocle.pseudotime", cols= viridis(100, begin = 0))+labs(color="Pseudotime")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),plot.title = element_blank(),
        axis.ticks=element_blank(),axis.line=element_blank(),axis.text=element_blank())
ggsave2("umap.monocle.pseudo.png",width=5.5, height=4,device="png")
VlnPlot(obj, features = "monocle.pseudotime",pt.size=0)+ NoLegend()+labs(y="Pseudotime")+
  theme(axis.title.x=element_blank(),plot.title=element_blank())
ggsave2("cluster.by.monoclepseudo.png",width=4, height=4,device="png")
VlnPlot(obj, features = "monocle.pseudo.rank",pt.size=0)+ NoLegend()+labs(y="Pseudotime Rank")+
  theme(axis.title.x=element_blank(),plot.title=element_blank())
ggsave2("cluster.by.monoclepseudo.rank.png",width=4, height=4,device="png")
save.image("temp.monocle.graph.RData")

#monocle clustering
cds_2d <- preprocess_cds(cds, num_dim = 50)
cds_2d <- reduce_dimension(cds_2d)
cds_2d <- cluster_cells(cds_2d)
marker_test_res <- top_markers(cds_2d, group_cells_by="cluster", 
                               reference_cells=1000, cores=8)
top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(3, pseudo_R2)
top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
plot_genes_by_group(cds_2d,
                    top_specific_marker_ids,
                    group_cells_by="cluster",
                    ordering_type="cluster_row_col",
                    max.size=3)
ggsave2("monocle.cluster.genes.png",width=5.5,height=8,device="png")
cds_2d <- learn_graph(cds_2d)
cds_2d <- order_cells(cds_2d, root_pr_nodes=get_earliest_principal_node(cds))
for(i in c("pseudotime","cluster","partition","my.clusters2")){
  plot_cells(cds_2d, color_cells_by = i,label_cell_groups=F,cell_size=0.3) 
  ggsave2(paste0(i,".moncole.2d.umap.png"),width=5.5,height=4,device="png")
}

#3D plot
cds_3d <- reduce_dimension(cds, max_components = 3)
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)
cds_3d <- order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds))
for(i in c("pseudotime","partition","cluster","my.clusters2")){
  plot_3d <- plot_cells_3d(cds_3d, color_cells_by=i) #clusters, pseudotime, partition
  if(i=="my.clusters2"){plot_3d <- plot_cells_3d(cds_3d,color_cells_by=i,color_palette = hue_pal()(14))}
  htmlwidgets::saveWidget(plot_3d,file=paste0(i,".monocle.3d.umap.html"))
}
save.image("temp2.monocle.RData")

# Gene regression:
gene_fits <- fit_models(cds, model_formula_str = "~my.clusters2")
fit_coefs <- coefficient_table(gene_fits)
condition_terms <- fit_coefs %>% filter(term!="(Intercept)")
condition_terms <- condition_terms %>% mutate(q_value = p.adjust(p_value))
sig_genes <- condition_terms %>% filter (q_value < 0.05) %>% pull(gene_short_name)
evaluate_fits(gene_fits)
write.xlsx(as.data.frame(condition_terms),file="condition_terms.xlsx")
write.xlsx(sig_genes,file="sig_gene.xlsx")
save(condition_terms,file="condition_terms.RData")
save(sig_genes,file="sig_genes.RData")
save(gene_fits,file="gene_fits.RData")
save.image("temp3.monocle.RData")

#Co-regulated gene modules
pr_graph_test_res <- graph_test(cds, neighbor_graph="knn", cores=4)
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=1e-2)
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), cell_group=colData(cds)$my.clusters2)
                                # cell_group=partitions(cds)[colnames(cds)])
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
# colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))
pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6)
plot_cells(cds, 
           genes=gene_module_df %>% filter(module %in% c(8, 28, 33, 37)),
           group_cells_by="partition",
           color_cells_by="partition",
           show_trajectory_graph=FALSE)
save.image("temp6.monocle.RData")

# With graph autocorrelation:
pr_test_res <- graph_test(cds,  neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(pr_test_res, q_value < 0.05))
plot_cells(cds, genes=c("hlh-4", "gcy-8", "dac-1", "oig-8"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=c(0,10^seq(-6,-1))) #trajectory modules
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$cell.type)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")
plot_cells(cds,
           genes=gene_module_df %>% filter(module %in% c(27, 10, 7, 30)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

#gene dynamics along path
AFD_genes <- c("CD274")
AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes,]
plot_genes_in_pseudotime(AFD_lineage_cds,
                         color_cells_by="my.clusters2",
                         min_expr=0.5)

save.image("temp7.monocle.RData")






