# HiTME :dart: :facepunch:

New tool for classifying all cell types in the tumor microenvirontment

# Annotate cell types
You can annotate cells with scGate and ProjecTILs in a parallel for-loop.
The function takes as input the path to the directory containing the seurat objects saved as .rds files (preferably each sample saved as separate .rds file). The directory should not contain other .rds files.
```r
scGate_models_DB <- get_scGateDB(branch = "master", verbose = T, force_update = TRUE)
models.TME <- scGate_models_DB$human$TME_HiRes

path_ref <- "~/Dropbox/CSI/reference_atlases"
ref.maps.list <- list(file.path(path_ref, "CD8T_human_ref_v1.rds"),
                      file.path(path_ref, "CD4T_human_ref_v2.rds"),
                      file.path(path_ref, "DC_human_ref_v1.rds"),
                      file.path(path_ref, "MoMac_human_v1.rds"))

annotate_cells(dir = file.path(path_output, "test"),
               scGate.model = models.TME,
               ref.maps = ref.maps.list,
               mc.cores = 3)
```


# Calculate cell type compositions
calc_CTcomp calculates the cell type composition from a seurat object with a metadata column containing the cell type annotations (e.g. called "celltype").
```r
calc_CTcomp(object = panc8, annot.cols = "celltype")
```

To calculate cell subtype composition, it is adviced to create a separate metadata column for each subtype composition, e.g.
- One metadata column called "celltype"
- One metadata column called "CD8_subtypes" containing only CD8 subtype annotations (all other cells as "NA")
- One metadata column called "CD4_subtypes" ...
```r
calc_CTcomp(object = panc8, annot.cols = c("celltype", "CD8_subtypes", "CD4_subtypes"))
```

You can also calculate the cell type composition for each sample by indicating a metadata column containing the sample names:
```r
calc_CTcomp(object = panc8, annot.cols = "celltype", sample.col = "Sample")
```

## Example
```r
devtools::install_github('satijalab/seurat-data')
library(SeuratData)
options(timeout = max(300, getOption("timeout")))
InstallData("panc8")
data("panc8")
panc8 = UpdateSeuratObject(object = panc8)

###################################################################
# Calculate overall composition
celltype.compositions.overall <- calc_CTcomp(object = panc8, annot.cols = "celltype")

# Calculate sample-wise composition
celltype.compositions.sample_wise <- calc_CTcomp(object = panc8, annot.cols = "celltype", sample.col = "orig.ident")

###################################################################
# Plot overall composition
par(mar = c(8.1, 4.1, 4.1, 2.1)) 
barplot(unlist(celltype.compositions.overall[["celltype"]][["freq"]]),
        ylab = "Relative abundance (%)",
        las = 2)


# Plot sample-wise composition
df <- celltype.compositions.sample_wise[["celltype"]][["freq"]]
dflong <- reshape2::melt(t(df))

library(ggplot2)
ggplot(dflong, aes(x=Var1, y=value, color=Var1)) +
  geom_boxplot() + geom_point() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + labs(x = "", y = "Relative abundance (%)") + NoLegend()
```

|Overall composition|Composition by sample|
|:-:|:-:|
|![Rplot1](https://github.com/carmonalab/HiTME/assets/67605347/43c81cfc-e4a2-42eb-8b8f-bf70da52a271)|![Rplot2](https://github.com/carmonalab/HiTME/assets/67605347/d876846d-b7c9-4ff4-af5c-a4bfefcbc29f)|
