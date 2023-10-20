devtools::install_github('satijalab/seurat-data')
library(SeuratData)
InstallData("pbmc3k")
data("pbmc3k")
pbmc3k = UpdateSeuratObject(object = pbmc3k)

source("R/main.R")
celltype.compositions <- calc_CTcomp(object = pbmc3k, annot.cols = "seurat_annotations")
no_NA_test <- calc_CTcomp(object = pbmc3k, annot.cols = "seurat_annotations", useNA = FALSE)
