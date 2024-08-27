# HiTME :dart: :facepunch:

<p align="center">

<img src="docs/HiTME_logo.png" height="100"/>

</p>

## High-resolution Tumor Micro-Environment cell type classification

HiTME is designed for precise cell type classification within the complex tumor microenvironment (TME), providing high accuracy and interpretability in cell type identification.

Find a vignette describing its main functions in [html](https://carmonalab.github.io/HiTME_CaseStudies/HiTME_GetStarted.html) and its [code (repository)](https://github.com/carmonalab/HiTME_CaseStudies).

## Installation

``` r
# install.packages("remotes")
remotes::install_github("carmonalab/HiTME")
```

## Usage

``` r
library(ProjecTILs)
library(HiTME)

# Fetch reference maps
ref.maps <- ref.maps <- get.reference.maps(collection = "human")

# Cell type classification on a Seurat object
query <- Run.HiTME(query,
                   ref.maps = ref.maps[["human"]])

# Seurat object metadata has been updated with cell type classification at different granularity levels (layers)
```

### How to cite HiTME
Please note that the publication describing HiTME is currently in preparation. In the meantime, we kindly ask that you cite the two primary components of HiTME in your work:

- [scGate](https://github.com/carmonalab/scGate): Andreatta, Massimo, Ariel J. Berenstein, and Santiago J. Carmona. 2022. “scGate: Marker-Based Purification of Cell Types from Heterogeneous Single-Cell RNA-Seq Datasets.” Bioinformatics 38 (April): 2642–44. https://doi.org/10.1093/BIOINFORMATICS/BTAC141.
- [ProjecTILs](https://github.com/carmonalab/ProjecTILs): Andreatta, Massimo, Jesus Corria-Osorio, Sören Müller, Rafael Cubas, George Coukos, and Santiago J. Carmona. 2021. “Interpretation of t Cell States from Single-Cell Transcriptomics Data Using Reference Atlases.” Nature Communications 2021 12:1 12 (May): 1–19. https://doi.org/10.1038/s41467-021-23324-4.

<br>

# Cell type annotation by HiTME

**HiTME is an R package that combines [scGate](https://github.com/carmonalab/scGate) and [ProjecTILs](https://github.com/carmonalab/ProjecTILs) to classify cell types in single-cell RNA-seq data at high resolution and with large flexibility (e.g. easy to include new cell types).**

The function takes as input `Seurat` objects (or list of them). These should be split by sample to avoid batch effects, or split internally in `Run.HitME` by indicating the parameter `split.by`.

**HiTME** firstly runs [scGate](https://github.com/carmonalab/scGate) (easily customizable) marker-based classification, resulting in a coarse-grained cell type classification (CD4T, B cell, Dendritic cell...). Next, it runs for each broad cell type [ProjecTILs](https://github.com/carmonalab/ProjecTILs) for a finer cell type classification (CD4+ TFH, Tex CD8+, cDC1...) based on cell mapping onto expert-curated single-cell reference maps. Finally, cell subtype are further classified based on gene programs such as cell cycling, IFN or HSP-response scoring, using [UCell](https://github.com/carmonalab/UCell).



<p align="center">

<img src="docs/HiTME_logo.png" height="100"/>

</p>
