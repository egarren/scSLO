rm(list=ls())
library(Seurat)
library(scater)
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

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

#load GSE192762 (mouse) data
load("merged.RData")
GSE192762_obj<-SLE.B.obj.combined
GSE192762_obj@meta.data<-GSE192762_obj@meta.data[,1:4]
GSE192762_obj@meta.data$Sample<-SLE.B.obj.combined@meta.data$mouse_ID

metadata<-read_excel("scrnaseq_metadata_all.xlsx") #define path
metadata[] <- lapply(metadata, as.character)

#load GSE136376 (mouse) data
samples<-metadata$Sample[metadata$GSE=="GSE136376"]
for (i in samples){
  test<-as.matrix(fread(paste0("/GSE136376/",i)),rownames=1)
  sum<-aggregate(test, list(row.names(test)), sum)
  rownames(sum)<-sum$Group.1
  sum2<-select(sum,-c(Group.1))
  mtx <- as(as.matrix(sum2), "sparseMatrix")   
  obj<-CreateSeuratObject(counts=mtx)
  obj@meta.data$Sample<-i
  assign(i,obj)
}
#merge
seurat.list<-lapply(samples,get)
GSE136376_obj<-merge(seurat.list[[1]], y=seurat.list[2:length(seurat.list)],add.cell.ids=samples,project="B")
save.image("obj2.RData")


#load GSE203132 (mouse) data
samples<-metadata$Sample[metadata$GSE=="GSE203132"]
for (i in samples){
  path<-paste0("/GSE203132/",i,"/outs/filtered_feature_bc_matrix")
  seurat.object<-CreateSeuratObject(counts = Read10X(data.dir = path,gene.column=1), min.cells = 3, min.features  = 200, project = i, assay = "RNA")
  seurat.object@meta.data$Sample<-i
  assign(i,seurat.object)
}
##Merge 
seurat.list<-lapply(samples,get)
GSE203132_obj<-merge(seurat.list[[1]], y=seurat.list[2:length(seurat.list)],add.cell.ids=samples,project="B")
#Ms entrez to Ms mgi
GSE203132_obj@assays[["RNA"]]@meta.features$original_ensembl<-rownames(GSE203132_obj@assays[["RNA"]]@meta.features)
genes.meta<-getBM(attributes=c("ensembl_gene_id", "mgi_symbol"),filters=
                    "ensembl_gene_id",values=list(rownames(GSE203132_obj@assays[["RNA"]]@meta.features)), mart=mouse,useCache=F) #useast.
m <- match(GSE203132_obj@assays[["RNA"]]@meta.features$original_ensembl, genes.meta$ensembl_gene_id)
GSE203132_obj@assays[["RNA"]]@meta.features<-cbind(GSE203132_obj@assays[["RNA"]]@meta.features,genes.meta[m,])
rownames(GSE203132_obj@assays[["RNA"]]@meta.features) <- make.names(GSE203132_obj@assays[["RNA"]]@meta.features$mgi_symbol,unique=T)
rownames(GSE203132_obj@assays[["RNA"]]@data) <- make.names(GSE203132_obj@assays[["RNA"]]@meta.features$mgi_symbol,unique=T)
rownames(GSE203132_obj@assays[["RNA"]]@counts) <- make.names(GSE203132_obj@assays[["RNA"]]@meta.features$mgi_symbol,unique=T)
GSE203132_obj <- subset(GSE203132_obj, features=rownames(GSE203132_obj[!(is.na(GSE203132_obj@assays[["RNA"]]@meta.features$mgi_symbol)),]))
GSE203132_obj <- subset(GSE203132_obj, features=rownames(GSE203132_obj[GSE203132_obj@assays[["RNA"]]@meta.features$mgi_symbol != "",]))
B_all<-GSE203132_obj
save.image("obj3.RData")

#merge mice
seurat.list<-lapply(c("GSE203132_obj","GSE136376_obj","GSE192762_obj"),get)
ms_all<-merge(seurat.list[[1]], y=seurat.list[2:length(seurat.list)],project="B")

#mgi to human gene name
ms_all@assays[["RNA"]]@meta.features$mgi_symbol<-rownames(ms_all@assays[["RNA"]]@meta.features)
hgnc_list<-getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = list(rownames(ms_all@assays[["RNA"]]@meta.features)) , 
                  mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host="dec2021.archive.ensembl.org"), attributesL = c("hgnc_symbol"), 
                  martL = useMart("ensembl", dataset = "hsapiens_gene_ensembl",host="dec2021.archive.ensembl.org"), uniqueRows=T)
m <- match(ms_all@assays[["RNA"]]@meta.features$mgi_symbol, hgnc_list$MGI.symbol)
ms_all@assays[["RNA"]]@meta.features<-cbind(ms_all@assays[["RNA"]]@meta.features,hgnc_list[m,])
ms_all@assays[["RNA"]]@meta.features$hgnc_symbol<-ms_all@assays[["RNA"]]@meta.features$HGNC.symbol
rownames(ms_all@assays[["RNA"]]@meta.features) <- make.names(ms_all@assays[["RNA"]]@meta.features$hgnc_symbol,unique=T)
rownames(ms_all@assays[["RNA"]]@data) <- make.names(ms_all@assays[["RNA"]]@meta.features$hgnc_symbol,unique=T)
rownames(ms_all@assays[["RNA"]]@counts) <- make.names(ms_all@assays[["RNA"]]@meta.features$hgnc_symbol,unique=T)
ms_all <- subset(ms_all, features=rownames(ms_all[!(is.na(ms_all@assays[["RNA"]]@meta.features$hgnc_symbol)),]))
ms_all <- subset(ms_all, features=rownames(ms_all[ms_all@assays[["RNA"]]@meta.features$hgnc_symbol != "",]))
save.image("obj4.RData")

#load GSE182221 (human) data
samples<-metadata$Sample[metadata$GSE=="GSE182221"]
for(i in samples){
  mtx<-ReadMtx(paste0("/GSE182221/",i,"/matrix.mtx.gz"),
               cells=paste0("/GSE182221/",i,"/barcodes.tsv.gz"),
               features=paste0("/GSE182221/",i,"/features.tsv.gz"),
               feature.column=1)
  seurat_obj<-CreateSeuratObject(counts=mtx)
  seurat_obj@meta.data$Sample<-i
  assign(i,seurat_obj)
}
#merge
seurat.list<-lapply(samples,get)
GSE182221_obj<-merge(seurat.list[[1]], y=seurat.list[2:length(seurat.list)],add.cell.ids=samples,project="B")
#Hu entrez to Hgnc
dim(GSE182221_obj)
GSE182221_obj@assays[["RNA"]]@meta.features$original_ensembl<-rownames(GSE182221_obj@assays[["RNA"]]@meta.features)
genes.meta<-getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"),filters=
                    "ensembl_gene_id",values=list(rownames(GSE182221_obj@assays[["RNA"]]@meta.features)),mart=human,useCache = F) #useast.
m <- match(GSE182221_obj@assays[["RNA"]]@meta.features$original_ensembl, genes.meta$ensembl_gene_id)
GSE182221_obj@assays[["RNA"]]@meta.features<-cbind(GSE182221_obj@assays[["RNA"]]@meta.features,genes.meta[m,])
rownames(GSE182221_obj@assays[["RNA"]]@meta.features) <- make.names(GSE182221_obj@assays[["RNA"]]@meta.features$hgnc_symbol,unique=T)
rownames(GSE182221_obj@assays[["RNA"]]@data) <- make.names(GSE182221_obj@assays[["RNA"]]@meta.features$hgnc_symbol,unique=T)
rownames(GSE182221_obj@assays[["RNA"]]@counts) <- make.names(GSE182221_obj@assays[["RNA"]]@meta.features$hgnc_symbol,unique=T)
GSE182221_obj <- subset(GSE182221_obj, features=rownames(GSE182221_obj[!(is.na(GSE182221_obj@assays[["RNA"]]@meta.features$hgnc_symbol)),]))
dim(GSE182221_obj)


#merge EVERYTHING
ALL_samples<-list(ms_all,GSE182221_obj)
B_all<-merge(ALL_samples[[1]], y=ALL_samples[2:length(ALL_samples)],project="SLO")
B_all@assays[["RNA"]]@meta.features$hgnc_symbol<-rownames(B_all@assays[["RNA"]]@meta.features)
save.image("obj5.RData")


# add sample metadata
metadata<-read_excel("scrnaseq_metadata_all.xlsx") #define path
metadata[] <- lapply(metadata, as.character)
B_all@meta.data<-B_all@meta.data %>% mutate(Sample=coalesce(Sample,orig.ident))
setdiff(metadata$Sample,unique(B_all@meta.data$Sample)) #check all samples present
B_all<-AddMetaData(B_all, metadata=rownames(B_all@meta.data),col.name = "cellID")
B_all@meta.data<-left_join(x = B_all@meta.data, y = metadata, by = "Sample")
rownames(B_all@meta.data)<-B_all@meta.data$cellID
rm(list=setdiff(ls(), c("B_all", "metadata")))
save.image("temp.merged.RData")





