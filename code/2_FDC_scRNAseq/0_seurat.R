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

wd<-"/scFDC/"
  
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

#load figshare data from
# https://figshare.com/articles/dataset/Human_tonsillar_stromal_cells_and_immune_cells/21325737
# https://figshare.com/articles/dataset/BRCs_and_immune_cells_from_murine_and_human_secondary_lymphoid_organs/21221291
# https://www.nature.com/articles/s41590-023-01503-3?utm_source=ni_etoc&utm_medium=email&utm_campaign=toc_41590_24_7&utm_content=20230708#data-availability
# https://www.nature.com/articles/s41590-023-01502-4#data-availability
setwd("/RNAseq_fastq/figshare")
df.list<-list()
for(i in list.files(pattern="rds")){
  name<-gsub("_sce.rds","",i)
  sce<-readRDS(i)
  # sce<-runPCA(sce)
  obj<- as.Seurat(sce, counts = NULL, data = "logcounts")
  obj<-RenameAssays(obj,originalexp='RNA')
  df.list[[i]]<-data.frame(dataset=unique(obj@meta.data$dataset),file=name)
  assign(name,obj)
}
df<-do.call("rbind",df.list)
write.csv(df,"figshare_meta.csv")
figshare_meta<-read_excel("figshare_meta2.xlsx")
datasets_keep<-figshare_meta$dataset[figshare_meta$file=="humBRC" & figshare_meta$Organ == "LN"]
humBRC<-subset(humBRC, subset= dataset %in% datasets_keep)

#merge figshare human
figshare_hu_samples<-unique(figshare_meta$file[figshare_meta$Species=="Hs"])
seurat.list<-lapply(figshare_hu_samples,get)
figshare_hu<-merge(seurat.list[[1]], y=seurat.list[2:length(seurat.list)],add.cell.ids=figshare_hu_samples,project="SLO")
figshare_hu@meta.data$Sample<-figshare_hu@meta.data$dataset
#Hu entrez to Hgnc
figshare_hu@assays[["RNA"]]@meta.features$original_ensembl_gene<-rownames(figshare_hu@assays[["RNA"]]@meta.features)
figshare_hu@assays[["RNA"]]@meta.features$original_ensembl<-gsub("\\..*","",rownames(figshare_hu@assays[["RNA"]]@meta.features))
genes.meta<-getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"),filters=
                    "ensembl_gene_id",values=list(figshare_hu@assays[["RNA"]]@meta.features$original_ensembl),mart=human,useCache = F) #useast.
m <- match(figshare_hu@assays[["RNA"]]@meta.features$original_ensembl, genes.meta$ensembl_gene_id)
figshare_hu@assays[["RNA"]]@meta.features<-cbind(figshare_hu@assays[["RNA"]]@meta.features,genes.meta[m,])
rownames(figshare_hu@assays[["RNA"]]@meta.features) <- make.names(figshare_hu@assays[["RNA"]]@meta.features$hgnc_symbol,unique=T)
rownames(figshare_hu@assays[["RNA"]]@data) <- make.names(figshare_hu@assays[["RNA"]]@meta.features$hgnc_symbol,unique=T)
rownames(figshare_hu@assays[["RNA"]]@counts) <- make.names(figshare_hu@assays[["RNA"]]@meta.features$hgnc_symbol,unique=T)
figshare_hu <- subset(figshare_hu, features=rownames(figshare_hu[!(is.na(figshare_hu@assays[["RNA"]]@meta.features$hgnc_symbol)),]))
figshare_hu <- subset(figshare_hu, features=rownames(figshare_hu[figshare_hu@assays[["RNA"]]@meta.features$hgnc_symbol != "",]))

#merge figshare mouse
figshare_ms_samples<-unique(figshare_meta$file[figshare_meta$Species=="Mus"])
seurat.list<-lapply(figshare_ms_samples,get)
figshare_ms<-merge(seurat.list[[1]], y=seurat.list[2:length(seurat.list)],add.cell.ids=figshare_ms_samples,project="SLO")
figshare_ms@meta.data$Sample<-figshare_ms@meta.data$dataset
#Ms entrez to Ms mgi
figshare_ms@assays[["RNA"]]@meta.features$original_ensembl_gene<-rownames(figshare_ms@assays[["RNA"]]@meta.features)
figshare_ms@assays[["RNA"]]@meta.features$original_ensembl<-gsub("\\..*","",rownames(figshare_ms@assays[["RNA"]]@meta.features))
genes.meta<-getBM(attributes=c("ensembl_gene_id", "mgi_symbol"),filters=
                    "ensembl_gene_id",values=list(figshare_ms@assays[["RNA"]]@meta.features$original_ensembl), mart=mouse,useCache=F) #useast.
m <- match(figshare_ms@assays[["RNA"]]@meta.features$original_ensembl, genes.meta$ensembl_gene_id)
figshare_ms@assays[["RNA"]]@meta.features<-cbind(figshare_ms@assays[["RNA"]]@meta.features,genes.meta[m,])
rownames(figshare_ms@assays[["RNA"]]@meta.features) <- make.names(figshare_ms@assays[["RNA"]]@meta.features$mgi_symbol,unique=T)
rownames(figshare_ms@assays[["RNA"]]@data) <- make.names(figshare_ms@assays[["RNA"]]@meta.features$mgi_symbol,unique=T)
rownames(figshare_ms@assays[["RNA"]]@counts) <- make.names(figshare_ms@assays[["RNA"]]@meta.features$mgi_symbol,unique=T)
figshare_ms <- subset(figshare_ms, features=rownames(figshare_ms[!(is.na(figshare_ms@assays[["RNA"]]@meta.features$mgi_symbol)),]))
figshare_ms <- subset(figshare_ms, features=rownames(figshare_ms[figshare_ms@assays[["RNA"]]@meta.features$mgi_symbol != "",]))
save.image(paste0(wd,"temp.load0.RData"))

#load STARsolo data
setwd("/RNAseq_fastq/STARsolo")
for(i in c("SRR14355140","SRR14355144","SRR14355148","SRR14355152")){
  mtx<-ReadMtx(paste0("./",i,"/Solo.out/Gene/filtered/matrix.mtx"),
               cells=paste0("./",i,"/Solo.out/Gene/filtered/barcodes.tsv"),
               features=paste0("./",i,"/Solo.out/Gene/filtered/features.tsv"),
               feature.column=1)
  seurat_obj<-CreateSeuratObject(counts=mtx)
  seurat_obj@meta.data$Sample<-i
  assign(i,seurat_obj)
}
save.image(paste0(wd,"temp.load1.RData"))

#load SCP data
setwd("/RNAseq_fastq/SCP/SCP2169/expression/644751ebb67df4748764a4a1")
mtx<-ReadMtx("matrix.mtx.gz",cells="barcodes.tsv.gz",features="features.tsv.gz")
obj<-CreateSeuratObject(counts=mtx)
meta_temp<-read.csv("../../metadata/HumanTonsil_metadata.csv")
meta_temp<-meta_temp[2:nrow(meta_temp),]
meta_temp$Sample<-"SCP2169"
rownames(meta_temp)<-meta_temp$NAME
obj<-AddMetaData(obj,meta_temp[,c("NAME","Sample")])
SCP2169<-obj

setwd("/RNAseq_fastq/SCP/SCP1423/expression/60c1451a771a5b0b9fc4f79e")
mtx<-ReadMtx("gene_sorted-MMRNA_PREDICT_3p_Paper_CD.mtx",cells="BarcodesCD_3p_paper.tsv",features="GenesCD_3p_paper.tsv")
obj<-CreateSeuratObject(counts=mtx)
meta_temp<-read.csv("../../metadata/AlexandriaUpload_AlexandriaMetadataConvention_PREDICT_3p_Paper_CD_final_cellOntologyLabels2.csv")
meta_temp<-meta_temp[2:nrow(meta_temp),]
meta_temp$Sample<-paste0("SCP1423_",meta_temp$donor_id)
rownames(meta_temp)<-meta_temp$NAME
write.csv(unique(meta_temp$Sample),file="SCP1423_meta.csv")
obj<-AddMetaData(obj,meta_temp[,c("NAME","Sample")])
SCP1423<-obj

setwd("/RNAseq_fastq/SCP/SCP1422/expression/60c1423c771a5b4ce3de841f")
mtx<-ReadMtx("gene_sorted-MMRNA_PREDICT_3p_Paper_FGIDFull.mtx",cells="PrefixedCustomBarcodes_FGID_3p_paper.tsv",features="CustomGenes_FGID_3p_paper.tsv")
obj<-CreateSeuratObject(counts=mtx)
meta_temp<-read.csv("../../metadata/AlexandriaUpload_AlexandriaMetadataConvention_PREDICT_3p_Paper_FGID_final_cellOntologyLabels2.csv")
meta_temp<-meta_temp[2:nrow(meta_temp),]
meta_temp$Sample<-paste0("SCP1422_",meta_temp$donor_id)
rownames(meta_temp)<-meta_temp$NAME
write.csv(unique(meta_temp$Sample),file="SCP1422_meta.csv")
obj<-AddMetaData(obj,meta_temp[,c("NAME","Sample")])
SCP1422<-obj

setwd("/RNAseq_fastq/SCP/SCP1186")
df<-read.table("./expression/LN_Figure6_RawData.txt")
m<- as(as.matrix(df), "sparseMatrix")  
obj<-CreateSeuratObject(counts=m)
meta_temp<-read.table("./metadata/LN_Figure6_Metadata.txt", sep = '\t',header=T)
meta_temp<-meta_temp[2:nrow(meta_temp),]
meta_temp$Sample<-paste0("SCP1186_",meta_temp$Array)
rownames(meta_temp)<-meta_temp$NAME
write.csv(unique(meta_temp$Sample),file="SCP1186_meta.csv")
unique(meta_temp[,c("Sample","Mouse.Phenotype")])
obj<-AddMetaData(obj,meta_temp[,c("NAME","Sample")])
SCP1186<-obj
save.image(paste0(wd,"temp.load2.RData"))

#load EBI matrix data
setwd("/RNAseq_fastq/EBI/E-MTAB-10206/Files/")

DC<-ReadMtx("matrix_DCs.mtx.gz",cells="barcodes_DCs.tsv.gz",features="features_DCs.tsv.gz")
DC_obj<-CreateSeuratObject(counts=DC)
DC_meta<-read.table("Human_DCs_meta.txt",header = T)
DC_meta$Sample<-paste0("E-MTAB-10206-",DC_meta$Donor)
rownames(DC_meta)<-DC_meta$Barcode
DC_obj<-AddMetaData(DC_obj,DC_meta[,c("Barcode","Sample")])

FRC<-ReadMtx("matrix.mtx.gz",cells="barcodes.tsv.gz",features="features.tsv.gz")
FRC_obj<-CreateSeuratObject(counts=FRC)
FRC_meta<-read.table("Human_FCRs_meta.txt",header = T)
FRC_meta$Sample<-paste0("E-MTAB-10206-",FRC_meta$Donor+3)
rownames(FRC_meta)<-FRC_meta$Barcode
FRC_obj<-AddMetaData(FRC_obj,FRC_meta[,c("Barcode","Sample")])


setwd("/RNAseq_fastq/EBI/E-MTAB-11644")
for(i in c("mandLN","mesLN","pLN")){
  mtx<-ReadMtx(paste0(i,"_matrix.mtx.gz"),cells=paste0(i,"_barcodes.tsv.gz"),features=paste0(i,"_features.tsv.gz"))
  seurat_obj<-CreateSeuratObject(counts=mtx)
  if(i=="mandLN"){seurat_obj@meta.data$Sample<-"E-MTAB-11644-1"}
  if(i=="mesLN"){seurat_obj@meta.data$Sample<-"E-MTAB-11644-2"}
  if(i=="pLN"){seurat_obj@meta.data$Sample<-"E-MTAB-11644-3"}
  assign(i,seurat_obj)
}


#Load GSM/GSE
setwd("/RNAseq_fastq/GSE213254_RAW")
for(i in c("GSM6576591","GSM6576592","GSM6576593")){
  mtx<-ReadMtx(paste0(i,"_matrix.tsv.gz"),cells=paste0(i,"_barcodes.tsv.gz"),features=paste0(i,"_features.tsv.gz"))
  seurat_obj<-CreateSeuratObject(counts=mtx)
  seurat_obj@meta.data$Sample<-i
  assign(i,seurat_obj)
}

#merge preloaded ms objects
seurat.list<-list(mandLN,mesLN,pLN,SCP1186,figshare_ms,GSM6576591,GSM6576592,GSM6576593)
pre_ms<-merge(seurat.list[[1]], y=seurat.list[2:length(seurat.list)],project="SLO")

setwd(wd)
save.image(paste0(wd,"temp.objects.RData"))

#load cellranger data
setwd(wd)
metadata<-read_excel("scrnaseq_metadata.xlsx") #define path
metadata[] <- lapply(metadata, as.character)
samples<-metadata$Sample
# check that all samples are present
for (i in samples){
  if(i!="tonsil1"){
  test<-read.table(paste0("/RNAseq_cellranger/FDC/",i,"/outs/filtered_feature_bc_matrix/features.tsv.gz"))}
}
for (i in samples){
  path<-paste0("/RNAseq_cellranger/FDC/",i,"/outs/filtered_feature_bc_matrix")
  if(i=="tonsil1"){path<-"/RNAseq_cellranger/FDC/tonsil1/outs/count/filtered_feature_bc_matrix"}
  seurat.object<-CreateSeuratObject(counts = Read10X(data.dir = path,gene.column=1), min.cells = 3, min.features  = 200, project = i, assay = "RNA")
  seurat_obj@meta.data$Sample<-i
  assign(i,seurat.object)
}
save.image("temp.objects2.RData")

##Merge SRA mouse count matrices
SRA_ms_samples<-metadata$Sample[metadata$Species=="Mus"]
seurat.list<-lapply(SRA_ms_samples,get)
SRA_ms<-merge(seurat.list[[1]], y=seurat.list[2:length(seurat.list)],add.cell.ids=SRA_ms_samples,project="SLO")
save.image("temp.objects3.RData")

#Ms entrez to Ms mgi
SRA_ms@assays[["RNA"]]@meta.features$original_ensembl<-rownames(SRA_ms@assays[["RNA"]]@meta.features)
genes.meta<-getBM(attributes=c("ensembl_gene_id", "mgi_symbol"),filters=
                    "ensembl_gene_id",values=list(rownames(SRA_ms@assays[["RNA"]]@meta.features)), mart=mouse,useCache=F) #useast.
m <- match(SRA_ms@assays[["RNA"]]@meta.features$original_ensembl, genes.meta$ensembl_gene_id)
SRA_ms@assays[["RNA"]]@meta.features<-cbind(SRA_ms@assays[["RNA"]]@meta.features,genes.meta[m,])
rownames(SRA_ms@assays[["RNA"]]@meta.features) <- make.names(SRA_ms@assays[["RNA"]]@meta.features$mgi_symbol,unique=T)
rownames(SRA_ms@assays[["RNA"]]@data) <- make.names(SRA_ms@assays[["RNA"]]@meta.features$mgi_symbol,unique=T)
rownames(SRA_ms@assays[["RNA"]]@counts) <- make.names(SRA_ms@assays[["RNA"]]@meta.features$mgi_symbol,unique=T)
SRA_ms <- subset(SRA_ms, features=rownames(SRA_ms[!(is.na(SRA_ms@assays[["RNA"]]@meta.features$mgi_symbol)),]))
SRA_ms <- subset(SRA_ms, features=rownames(SRA_ms[SRA_ms@assays[["RNA"]]@meta.features$mgi_symbol != "",]))

#merge SRA and preloaded mice
ms_all<-merge(SRA_ms, y=pre_ms,project="SLO")
save.image("temp.objects4.RData")

#mgi to human gene name
ms_all@assays[["RNA"]]@meta.features$mgi_symbol<-rownames(ms_all@assays[["RNA"]]@meta.features)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
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
save.image("temp.objects5.RData")

#merge SRA human
SRA_hu_samples<-metadata$Sample[metadata$Species=="Hs"]
seurat.list<-lapply(SRA_hu_samples,get)
SRA_hu<-merge(seurat.list[[1]], y=seurat.list[2:length(seurat.list)],add.cell.ids=SRA_hu_samples,project="SLO")

#Hu entrez to Hgnc
SRA_hu@assays[["RNA"]]@meta.features$original_ensembl<-rownames(SRA_hu@assays[["RNA"]]@meta.features)
genes.meta<-getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"),filters=
                    "ensembl_gene_id",values=list(rownames(SRA_hu@assays[["RNA"]]@meta.features)),mart=human,useCache = F) #useast.
m <- match(SRA_hu@assays[["RNA"]]@meta.features$original_ensembl, genes.meta$ensembl_gene_id)
SRA_hu@assays[["RNA"]]@meta.features<-cbind(SRA_hu@assays[["RNA"]]@meta.features,genes.meta[m,])
rownames(SRA_hu@assays[["RNA"]]@meta.features) <- make.names(SRA_hu@assays[["RNA"]]@meta.features$hgnc_symbol,unique=T)
rownames(SRA_hu@assays[["RNA"]]@data) <- make.names(SRA_hu@assays[["RNA"]]@meta.features$hgnc_symbol,unique=T)
rownames(SRA_hu@assays[["RNA"]]@counts) <- make.names(SRA_hu@assays[["RNA"]]@meta.features$hgnc_symbol,unique=T)
SRA_hu <- subset(SRA_hu, features=rownames(SRA_hu[!(is.na(SRA_hu@assays[["RNA"]]@meta.features$hgnc_symbol)),]))
SRA_hu <- subset(SRA_hu, features=rownames(SRA_hu[SRA_hu@assays[["RNA"]]@meta.features$hgnc_symbol != "",]))
save.image("temp.objects6.RData")

#merge EVERYTHING
SRA_hu@meta.data$Sample<-SRA_hu@meta.data$orig.ident
SCP2169@meta.data$Sample<-"SCP2169"
ALL_samples<-list(ms_all,SRA_hu,DC_obj,FRC_obj,SCP2169,SCP1423,SCP1422,SRR14355140,SRR14355144,SRR14355148,SRR14355152,figshare_hu)
SLO_all<-merge(ALL_samples[[1]], y=ALL_samples[2:length(ALL_samples)],project="SLO")
SLO_all@assays[["RNA"]]@meta.features$hgnc_symbol<-rownames(SLO_all@assays[["RNA"]]@meta.features)
save.image("temp.objects7.RData")

# add sample metadata
metadata<-read_excel("scrnaseq_metadata_all.xlsx") #define path
metadata[] <- lapply(metadata, as.character)
SLO_all@meta.data<-SLO_all@meta.data %>% mutate(Sample=coalesce(Sample,orig.ident))
setdiff(metadata$Sample,unique(SLO_all@meta.data$Sample)) #check all samples present
SLO_all<-AddMetaData(SLO_all, metadata=rownames(SLO_all@meta.data),col.name = "cellID")
SLO_all@meta.data<-left_join(x = SLO_all@meta.data, y = metadata, by = "Sample")
rownames(SLO_all@meta.data)<-SLO_all@meta.data$cellID
rm(list=setdiff(ls(), c("SLO_all", "metadata")))
save.image("temp.merged.RData")





