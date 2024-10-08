% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{plot.geneGating}
\alias{plot.geneGating}
\title{Plot gene expression along scGate model}
\usage{
\method{plot}{geneGating}(object = NULL, scGate.model = NULL, group.by = NULL, split.by = NULL)
}
\arguments{
\item{object}{Seurat object or list of them, if list is provided one plot per Seurat object is returned}

\item{scGate.model}{scGate model to plot their gating genes}

\item{group.by}{Whether to group gene expression by a variable in metadata}

\item{split.by}{Whether to generate one plot for each grouping variable provided}
}
\value{
Violin plots to evaluate the gene expression of scGate models gating genes.
}
\description{
Plot gene expression along scGate model
}
