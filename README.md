# scSLO

This repository contains the code used in our paper: "PD-L1 controls germinal center dynamics"

## Installation guide
Install dependencies listed below

## Demo
Model data is included in the `data` directory.

## Instructions
1. Download `data` and `code` directories
2. Set `data` as the working directory
3. Download mouse peyer's patch scRNA-seq data from [GSE260776](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE260776)
4. Download human tonsil scRNA-seq data from [GSE262278](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262278)
5. Download mouse TFH and TFR bulk RNA-seq data from [GSE263423](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263423)
6. Download mouse TFH and TFR scRNA-seq data from [GSE157649](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157649)
7. Download mouse GCB scRNA-seq data from [GSE203132](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203132)
8. Download B6.Sle1yaa scRNA-seq data [GSE192762](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE192762)
9. Download remaining mouse and human scRNA-seq data from sources listed in Table S1 and Table S2.
10. Run each script in numerical order.  These scripts will generate the figures presented in our manuscript.\
\
NB: Expected total run time is 3-5 days

## System requirements and software
[cellranger](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-gex-count) v7.0.0\
[TrackMate](https://imagej.net/plugins/trackmate/) v7.11.1\
[salmon](https://github.com/COMBINE-lab/salmon) v1.10.1\
[STARsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md) v2.7.11a

<ins>Python (v3.12.2) packages</ins>\
scanpy_1.10.1\
anndata2ri_1.3.1\
velocyto_0.17.17\
scVelo_0.3.2

<ins>R (v4.2.3) packages</ins>\
TrackMateR_0.3.9\
plotly_4.10.4\
tximport_1.30.0\
DESeq2_1.42.1\
Seurat_5.0.0\
biomaRt_2.58.2\
harmony_0.1.0\
slingshot_2.10.0\
clusterProfiler_4.10.1\
msigdb_1.10.0\
pheatmap_1.0.12\
ggpubr_0.6.0\
ggplot2_3.5.0



