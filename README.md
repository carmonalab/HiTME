# HiTME :dart: :facepunch:

## High-resolution Tumor Micro-Environment cell type classification and compositional analysis for single-cell RNA-seq

HiTME is designed for precise cell type classification within the complex tumor microenvironment (TME), providing high accuracy and interpretability in cell type identification.

Find a vignette describing its main functions in [html](https://carmonalab.github.io/HiTME_CaseStudies/HiTME_demo.html) and its [code (repository)](https://github.com/carmonalab/HiTME_CaseStudies).

### Installation

``` r
remotes::install_github("carmonalab/HiTME")
```

### How to cite HiTME
Please note that the publication describing HiTME is currently in preparation. In the meantime, we kindly ask that you cite the two primary components of HiTME in your work:

- [scGate](https://github.com/carmonalab/scGate): Andreatta, Massimo, Ariel J. Berenstein, and Santiago J. Carmona. 2022. “scGate: Marker-Based Purification of Cell Types from Heterogeneous Single-Cell RNA-Seq Datasets.” Bioinformatics 38 (April): 2642–44. https://doi.org/10.1093/BIOINFORMATICS/BTAC141.
- [ProjecTILs](https://github.com/carmonalab/ProjecTILs): Andreatta, Massimo, Jesus Corria-Osorio, Sören Müller, Rafael Cubas, George Coukos, and Santiago J. Carmona. 2021. “Interpretation of t Cell States from Single-Cell Transcriptomics Data Using Reference Atlases.” Nature Communications 2021 12:1 12 (May): 1–19. https://doi.org/10.1038/s41467-021-23324-4.

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

<br>

# Summarized cell annotation

`Run.HiTME` will return the Seurat object or list of them with new metadata indicating cell type annotation.

Annotated Seurat objects can be summarized into **HiT objects** using `get.HiTObject` function. For this function the grouping variable `group.by` resulting from `Run.HiTME` annotation or additional annotations need to be indicated. Compositional cell type distribution and aggregated transcriptomic profile (pseudobulk) are returned for each sample.

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

1.  Seurat object metadata (dataframe): `metadata`

2.  Cell type predictions for each cell in the data set (list): `predictions`

3.  Cell type composition for each layer of cell type prediction: `composition`. Including:


    3.1. cell counts

    3.2. frequency

    3.3. CLR (Centred log ratio)-transformed counts (useful for downstream analyses such as PCA/[Logratio analysis](https://doi.org/10.1146/annurev-statistics-042720-124436) )

4.  Aggregated profile of predicted cell types: `aggregated_profile`. Including:

    4.1. Average and aggregated expression per cell type of all genes in the dataset and a subset of them.

    4.2. Mean of [UCell](https://github.com/carmonalab/UCell) scores per cell type, if additional signatures are provided, for example from [SignatuR](https://github.com/carmonalab/SignatuR).
