# scSLO

This repository contains the code used in our paper: "PD-L1 controls germinal center dynamics"

## Installation guide
Install dependencies listed below

## Demo
Model data is included in the `data` directory.

## Instructions
1. Download `data` and `code` directories
2. Set `data` as the working directory
3. Download mouse peyer's patch scRNA-seq data from [xxx](xxxx)
4. Download human tonsil scRNA-seq, scTCR-seq, and scBCR-seq data from [xxx](xxxx)
5. Download mouse TFH and TFR bulk RNA-seq data from [xxx](xxxx)
6. Download mouse TFH and TFR scRNA-seq data from [GSE157649](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157649)
7. Download mouse GCB scRNA-seq data from [GSE203132](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203132)
8. Download B6.Sle1yaa scRNA-seq data [GSE192762](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE192762)
9. Download remaining mouse and human scRNA-seq data from sources listed in Table S1 and Table S2.
10. Run each script in numerical order.  These scripts will generate the figures presented in our manuscript.\
\
NB: Expected total run time is 3-5 days

## System requirements and software
[cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/multi) v7.0.0

<ins>Python (v3.7.4) packages</ins>

<ins>R (v4.1.1) packages</ins>\
harmony_0.1.0\
Seurat_4.2.0        
clusterProfiler_4.6.0\
SPIA_2.50.0\
qgraph_1.9.2 \
ggfortify_0.4.15 \
ggpubr_0.4.0.999 \
ggplot2_3.3.6 \
pheatmap_1.0.12 


