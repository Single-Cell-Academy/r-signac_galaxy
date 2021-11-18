# Introduction
This repository contains code for Galaxy wrappers for the r-signac package to analyze scATAC-seq data.

!!! This package is currently in open development so check back recently to find new functions integrated. All functions will be available via the Galaxy test and main toolsheds once completed. !!!

# Processing raw ATAC-seq data
To process your raw scATAC-seq fastq files, you can use CellRanger-ATAC, either from the command line or you can use our CellRanger-ATAC galaxy tool, which you can find here:

GenAP2 CellRanger-ATAC: https://github.com/Single-Cell-Academy/genap2_cellranger_ATAC

# External sites
Original Signac project: https://github.com/timoast/signac
Signac vignettess for R: https://satijalab.org/signac/articles/pbmc_vignette.html

# Additional information
Due to the related nature of Seurat and Signac, the structure and organization of the Galaxy and R code in this repo is attempting to follow similar principles established by the EBI single-cell tools. 
Some of the code for running Rscripts via the command line have been directly copied over from their repositories which can be found here:
r-seurat-scripts: https://github.com/ebi-gene-expression-group/r-seurat-scripts/blob/develop/seurat-read.R
Seurat galaxy wrappers: https://github.com/ebi-gene-expression-group/container-galaxy-sc-tertiary/tree/develop/tools/tertiary-analysis/seurat