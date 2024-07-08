# HiTME :dart: :facepunch:

## High-resolution Tumor Micro-Environment cell type classification

HiTME is designed for precise cell type classification within the complex tumor microenvironment (TME), providing high accuracy and interpretability in cell type identification.

Find a vignette describing its main functions in [html](https://carmonalab.github.io/HiTME_CaseStudies/HiTME_demo.html) and its [code (repository)](https://github.com/carmonalab/HiTME_CaseStudies).

### Installation

``` r
remotes::install_github("carmonalab/HiTME")
```

<br>

# Cell type annotation

**HiTME is an R package that combines [scGate](https://github.com/carmonalab/scGate) and [ProjecTILs](https://github.com/carmonalab/ProjecTILs) to classify cell types in single-cell RNA-seq data at high resolution and with large flexibility (e.g. easy to include new cell types).**

The function takes as input `Seurat` objects (or list of them). These should be split by sample to avoid batch effects, or split internally in `Run.HitME` by indicating the parameter `split.by`.

This wrapper firstly runs [scGate](https://github.com/carmonalab/scGate) (easily customizable) marker-based classification, resulting in a coarse-grained cell type classification (CD4T, B cell, Dendritic cell...). Next, it runs for each broad cell type [ProjecTILs](https://github.com/carmonalab/ProjecTILs) for a finer cell type classification (CD4+ TFH, Tex CD8+, cDC1...) based on cell mapping onto expert-curated single-cell reference maps.

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
```

<br>

By default [scGate](https://github.com/carmonalab/scGate) (layer 1) will return the [cell ontology ID](https://www.ebi.ac.uk/ols4/ontologies/cl) for each predicted cell type. This ID will be then used to link each coarse cell type with its respective reference map for finer cell type classification using [ProjecTILs](https://github.com/carmonalab/ProjecTILs). Hence, we need to indicate each respective cell ontology ID(s) for each reference map.

If alternative cell type link are used between the coarse and finer cell type classification, this must be specified in `Run.HiTME` using `layer1_link` parameter.

``` r
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
```

``` r
# Run HiTME
annotated.obj <- Run.HiTME(object = obj,
                scGate.model = models.TME,
                ref.maps = ref.maps)

annotated.obj <- Run.HiTME(obj,
                            scGate.model = models.TME,
                            ref.maps = ref.maps,
                            # already split object
                            split.by = NULL,
                            # if splitting or providing list, whether to return a single merged object
                            remerge = FALSE,
                            # link between scGate and ProjecTILs
                            layer1_link = "CellOntology_ID",
                            # extra signatures to be computed per celltype
                            additional.signatures = additional.signatures, 
                            # paralelization parameters
                            ncores = 4,
                            progressbar = TRUE
                            )
```
