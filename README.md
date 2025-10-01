# Analysis of CrebA mutant scRNA-seq and CrebA ChIP-seq 
This repository contains the analysis code for our [manuscript](https://www.biorxiv.org/content/10.1101/2025.06.12.659381v1.abstract). The code base for the analysis of the scRNA-seq can be found in [*analysis*](/analysis/) folder. Before running the analysis code, download the zip file from [here]() and unzip the content into a directory called `quantification`. These files are the output files from CellRanger. 

The code base for the analysis of ChIP-seq and scRNA-seq integration can be found in [*CrebA_ChIP_Analysis*](/CrebA_ChIP_Analysis/). 

## Processed data 
The raw count matrices and processed data (Seurat object in .rds format) for stage 10-12 (early) CrebA mutant embryos and stage 13-16 (late) CrebA mutant embryos can be found on [GEO (GSE306121)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306121). The Seurat version used is 4.3.0. 
|Stage| Seurat object |
| --------- | --------------- |
| stage 10-12 | [download from GEO](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE306121&format=file&file=GSE306121%5Fearly%5FCrebA%5Fseurat%5Fobject%2Erds) |
| stage 13-16 | [download from GEO](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM9193590&format=file&file=GSM9193590%5Flate%5FCrebA%5Fseurat%5Fobject%2Erds) |
