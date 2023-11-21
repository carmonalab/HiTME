# HiTME :dart: :facepunch:

New tool for classifying all cell types in the tumor microenvironment

# TO DO ❗

``` r
# Installation for devs. From Rstudio Terminal, run the following code:
# R CMD build /directory_of_your_downloaded_project_/HiTME
# R CMD install HiTME_0.1.tar.gz

# Then load as usual
library(HiTME)
```

-   Provide better example dataset for public use.
-   Update example plots

<br>

### Installation

``` r
remotes::install_github("carmonalab/HiTME")
```

<br>

# Cell type annotation

**HiTME's `Run.HiTME` is a wrapper of [scGate](https://github.com/carmonalab/scGate) and [ProjecTILs](https://github.com/carmonalab/ProjecTILs) to classify cell types in single-cell RNA-seq experiments.**

The function takes as input `Seurat` objects (or list of them). These should be split by sample to avoid batch effects, or split internally in `Run.HitME` by indicating the parameter `split.by`.

This wrapper firstly runs [scGate](https://github.com/carmonalab/scGate) on TME (Tumor micronenvirontment) default models or alternatively on the models provided, resulting in a coarse cell type classification (CD4T, B cell, Dendritic cell...). Next, it runs [ProjecTILs](https://github.com/carmonalab/ProjecTILs) for a finer cell type classification (CD4+ TFH, Tex CD8+, cDC1...) based on the references provided on the cell types classified by [scGate](https://github.com/carmonalab/scGate) that are linked to a respective reference map.

``` r
library(scGate)
library(ProjecTILs)
library(HiTME)

# If multiple samples are within the same Seurat object, split by sample.
# obj.list <- SplitObject(obj, split.by = "Sample")

# Define scGate model if other than default is wanted
scGate_models_DB <- get_scGateDB(branch = "master")
models.TME <- scGate_models_DB$human$HiTME

# Load ProjecTILs reference maps
path_ref <- "~/reference_atlases"
ref.maps <- list(CD8 = load.reference.map(file.path(path_ref, "CD8T_human_ref_v1.rds")),
                 CD4 = load.reference.map(file.path(path_ref, "CD4T_human_ref_v2.rds")),
                 DC = load.reference.map(file.path(path_ref, "DC_human_ref_v1.rds")),
                 MoMac = load.reference.map(file.path(path_ref, "MoMac_human_v1.rds"))
                 )

# add scGate_link to ref.maps
# Include a slot in @misc with the cell name output by scGate
# By default scGate returns cell ontology ID

layer1.links <- list("CD8" = "CL:0000625",
                  "CD4" = "CL:0000624",
                  "DC" = "CL:0000451",
                  "MoMac" = "CL:0000576_CL:0000235"
                  )
                  
for(a in names(ref.maps)){
  ref.maps[[a]]@misc$layer1_link <- layer1.links[[a]]
}

# Run HiTME
annotated.obj <- Run.HiTME(object = obj,
                scGate.model = models.TME,
                ref.maps = ref.maps)

# Alternatively HiTME can be run splitting by sample
annotated.obj <- Run.HiTME(object = obj,
                          scGate.model = models.TME,
                          ref.maps = ref.maps,
                          split.by = "sample")
```

<br>

# Summarized cell annotation

`Run.HiTME` will return the Seurat object or list of them with new metadata indicating cell type annotation.

Annotated Seurat objects can be summarized into HiT objects using `get.HiTObject` function. For this function the grouping variable `group.by` resulting from `Run.HiTME` annotation or additional annotations need to be indicated. Compositional cell type distribution and aggregated transcriptomic profile are returned for each sample.

``` r
HiT_summary <- get.HiTObject(annotated.obj ,
                            group.by = list("layer1" = "scGate_multi",
                                            "layer2" = "functional.cluster"))
```

Alternatively, HiT summarizing object can be obtained directly using `Run.HiTME` with parameters `return.Seurat = FALSE`.

``` r
HiT_summary <- Run.HiTME(object = obj,
                scGate.model = models.TME,
                ref.maps = ref.maps,
                return.Seurat = FALSE)
```

## Hit Object content

The Hit object summarize the cell type annotation and contain the following slots:

1.  Seurat object metadata: `metadata`

2.  Cell type predictions for each cell in the data set: `predictions`

3.  Cell type composition for each layer of cell type prediction: `composition`. Including:

    3.1.  cell counts 
    
    3.2.  frequency 
    
    3.3.  CLR (Central log ratio)-transformed frequency
    

    ``` r
    # Run by

    get.celltype.composition(object = NULL,
                            group.by.composition = NULL,
                            split.by = NULL,
                            min.cells = 10,
                            useNA = FALSE,
                            clr_zero_impute_perc = 1
                            )
    ```

4.  Aggregated profile of predicted cell types: `aggregated_profile`. Including:

    4.1. Average and aggregated expression per cell type of all genes in the dataset and a subset of them.

    ``` r
    # Run by

    get.aggregated.profile(object,
                          group.by.aggregated = NULL,
                          gene.filter = NULL, # list of genes to subset
                          GO_accession = NULL, # GO accessions to fetch for subsetting
                          nHVG = 1000, # number of highly variable genes to subset
                          assay = "RNA",
                          layer = "data",
                          useNA = FALSE,
                        )

    ```

    4.2. Mean of signature scores per cell type, if additional signatures are provided, for example from [SignatuR](https://github.com/carmonalab/SignatuR).

    ``` r
    # Run by

    get.aggregated.signature(object,
                            group.by.aggregated = NULL,
                            name.additional.signatures = NULL,
                            fun = mean, # function for aggregating cell-wise signature scores
                            useNA = FALSE
                            )
    ```
