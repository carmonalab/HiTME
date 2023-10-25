# HiTME :dart: :facepunch:

New tool for classifying all cell types in the tumor microenvironment

# TO DO ❗
```r
# Installation for devs. From Rstudio Terminal, run the following code:
# R CMD build /directory_of_your_downloaded_project_/HiTME
# R CMD install HiTME_0.1.tar.gz

# Then load as usual
# library(HiTME)
```
- Provide better example dataset for public use.
- Update example plots
<br>

### Installation
To `annotate_cells` you need to install [ProjecTILs](https://github.com/carmonalab/ProjecTILs) from GitHub:
```r
remotes::install_github("carmonalab/ProjecTILs")
```
<br>

# Annotate sample cell types
You can `annotate_cells` with [scGate](https://github.com/carmonalab/scGate) and [ProjecTILs](https://github.com/carmonalab/ProjecTILs) in a parallel for-loop from samples stored on disk for large datasets.
The function takes as input the path to the directory containing the seurat objects saved as .rds files (preferably each sample saved as separate .rds file). The directory should not contain other .rds files.
For each `ProjecTIL` reference map (provided in named list), a new metadata column will be created called x_subtypes, where x = reference map name from list, e.g. CD8_subtypes, CD4_subtypes, ...

```r
# Create separate .rds files for each sample
obj.list <- SplitObject(obj, split.by = "Sample")

# Use helper functions to save/load your object list to/from disk
save_objs(obj.list, "./output/samples")
obj.list <- read_objs("./output/samples")
```

```r
# Example data
path_data <- file.path("~/Dropbox/CSI/Standardized_SingleCell_Datasets/ZhangY_2022_34653365/output/samples_subset")

# Define scGate model
scGate_models_DB <- get_scGateDB(branch = "master", verbose = T, force_update = TRUE)
models.TME <- scGate_models_DB$human$TME_HiRes

# Load ProjecTILs reference maps
path_ref <- "~/Dropbox/CSI/reference_atlases"
ref.maps <- list(CD8 = load.reference.map(file.path(path_ref, "CD8T_human_ref_v1.rds")),
                 CD4 = load.reference.map(file.path(path_ref, "CD4T_human_ref_v2.rds")),
                 DC = load.reference.map(file.path(path_ref, "DC_human_ref_v1.rds")),
                 MoMac = load.reference.map(file.path(path_ref, "MoMac_human_v1.rds")))

annotate_cells(dir = path_data,
               scGate.model = models.TME,
               ref.maps = ref.maps,
               ncores = 6)
```
<br>

# Calculate cell type compositions
`calc_CTcomp` calculates the cell type composition from one or multiple metadata column `annot.cols` containing the cell type annotations (e.g. called "scGate_multi", "CD8_subtypes", ...).
Accepted inputs are:
- `object` = a single Seurat object
- `object` = a list of Seurat objects
- `dir` = a directory containing Seurat objects saved as .rds (`object` parameter will be ignored)

To calculate cell subtype composition, it is adviced to create a separate metadata column for each subtype composition, e.g.
- One metadata column called "scGate_multi"
- One metadata column called "CD8_subtypes" containing only CD8 subtype annotations (all other cells as "NA")
- One metadata column called "CD4_subtypes" ...
```r
# Load data
path_data <- file.path("~/Dropbox/CSI/Standardized_SingleCell_Datasets/ZhangY_2022_34653365/output/samples_subset")
obj.list <- read_objs(path_data)

# For a single Seurat object (one sample)
celltype.compositions <- calc_CTcomp(obj.list[[1]])
# For a single Seurat object, multiple cell type columns
celltype.compositions <- calc_CTcomp(obj.list[[1]], annot.cols = c("scGate_multi", "CD8_subtypes"))

# For a single Seurat object (containing multiple samples)
obj <- merge(obj.list[[1]],obj.list[2:length(obj.list)])
celltype.compositions <- calc_CTcomp(obj, sample.col = "Sample")

# For a list of Seurat objects
celltype.compositions <- calc_CTcomp(obj.list)
celltype.compositions <- calc_CTcomp(obj.list, annot.cols = c("scGate_multi", "CD8_subtypes"))

# From a directory with Seurat objects (e.g. samples), without loading all into memory but instead looping over the single files
celltype.compositions <- calc_CTcomp(dir = path_data)
celltype.compositions <- calc_CTcomp(dir = path_data, annot.cols = c("scGate_multi", "CD8_subtypes"))
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
