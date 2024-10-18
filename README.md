# HiTME :dart: :facepunch:

<p align="center">

<img src="docs/HiTME_logo.png" height="140"/>

</p>

## High-resolution Tumor Micro-Environment cell type classification

**HiTME** is a robust tool for precise **cell type classification** in single-cell RNA-seq data, offering high resolution and flexibility. It provides accurate and interpretable cell type identification, with the ability to easily incorporate new cell types.

#### Key features of HiTME:

1.  **User-friendly:** Seamlessly integrates with the Seurat pipeline.

2.  **Hierarchical granularity**: Three interconnected levels of cell type annotation.

3.  **Tailored for the Tumor Micro-Environment**: Pre-configured with relevant cell types and granularity levels.

4.  **Flexible**: Easily customizable to add new cell types or states at any level.

5.  **Fast**: Processes samples of 5-10K cells in \~1 minute.

<br>

## Installation

``` r
# install.packages("remotes")
remotes::install_github("carmonalab/HiTME")
```

<br>

## Usage

``` r
library(ProjecTILs)
library(HiTME)

# Fetch reference maps
ref.maps <- get.reference.maps(collection = "human")

# Cell type classification on a Seurat object
query <- Run.HiTME(query,
                   ref.maps = ref.maps[["human"]],
                   species = "human" # for mouse data set to "mouse"
                   )

# Seurat object metadata has been updated with cell type classification at
# different granularity levels
```

Find a vignette describing its main functions in [html](https://carmonalab.github.io/HiTME_CaseStudies/HiTME_GetStarted.html) and its [code (repository)](https://github.com/carmonalab/HiTME_CaseStudies).

<br>

## How to cite HiTME

Please note that the publication describing HiTME is currently in preparation. In the meantime, we kindly ask that you cite the two primary components of HiTME in your work:

-   [scGate](https://github.com/carmonalab/scGate): Andreatta, Massimo, Ariel J. Berenstein, and Santiago J. Carmona. 2022. "scGate: Marker-Based Purification of Cell Types from Heterogeneous Single-Cell RNA-Seq Datasets." Bioinformatics 38 (April): 2642--44. <https://doi.org/10.1093/BIOINFORMATICS/BTAC141>.
-   [ProjecTILs](https://github.com/carmonalab/ProjecTILs): Andreatta, Massimo, Jesus Corria-Osorio, Sören Müller, Rafael Cubas, George Coukos, and Santiago J. Carmona. 2021. "Interpretation of t Cell States from Single-Cell Transcriptomics Data Using Reference Atlases." Nature Communications 2021 12:1 12 (May): 1--19. <https://doi.org/10.1038/s41467-021-23324-4>.

<br>

<p align="center">

<img src="docs/HiTME_logo.png" height="100"/>

</p>
