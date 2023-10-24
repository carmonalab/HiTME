

#' Annotate cell types in parallel
#'
#' @param dir The directory where your sample .rds files(seurat objects) are located
#' @param scGate.model The scGate model to use (use get_scGateDB() to get a list of available models)
#' @param ref.maps A named list of the ProjecTILs reference maps to use
#' @param ncores The number of cores to use
#' @param progressbar Whether to show a progressbar or not
#'
#' @return No return in R. The input files will be overwritten.
#' @export annotate_cells
#'
#' @examples
#' path_data <- file.path("~", "Dropbox", "CSI", "Standardized_SingleCell_Datasets", "ZhangY_2022_34653365", "output", "test")
#'
#' # Define scGate model
#' scGate_models_DB <- get_scGateDB(branch = "master", verbose = T, force_update = TRUE)
#' models.TME <- scGate_models_DB$human$TME_HiRes
#'
#' # Load ProjecTILs reference maps
#' path_ref <- "~/Dropbox/CSI/reference_atlases"
#' ref.maps <- list(CD8 = load.reference.map(file.path(path_ref, "CD8T_human_ref_v1.rds")),
#'                  CD4 = load.reference.map(file.path(path_ref, "CD4T_human_ref_v2.rds")),
#'                  DC = load.reference.map(file.path(path_ref, "DC_human_ref_v1.rds")),
#'                  MoMac = load.reference.map(file.path(path_ref, "MoMac_human_v1.rds")))
#'annotate_cells(dir = path_data,
#'               scGate.model = models.TME,
#'               ref.maps = ref.maps,
#'               ncores = 6)
annotate_cells <- function(dir, scGate.model, ref.maps,
                           ncores = 1, progressbar = TRUE){
  files <- list.files(dir)
  if (!all(endsWith(files, '.rds'))) {
    stop(paste("There are some files not ending on '.rds'"))
  }

  if (!is.list(ref.maps)) {
    stop("Please provide ref.maps as named list, containing the reference map(s), even if it is just one. The name will be used for the metadata column containing the cell type annotations")
  }

  param <- BiocParallel::MulticoreParam(workers = ncores, progressbar = progressbar)

  print("Running scGate")
  BiocParallel::bplapply(
    X = files,
    BPPARAM =  param,
    FUN = function(file) {
         path <- file.path(dir, file)
         x <- readRDS(path)
         x <- scGate(x, model=models.TME)
         saveRDS(x, path)
      }
  )
  print("Finished scGate")

  for (i in 1:length(ref.maps)) {
    ref.map.name <- names(ref.maps)[i]
    print(paste("Running ProjecTILs with map:  ", ref.map.name))
    BiocParallel::bplapply(
      X = files,
      BPPARAM =  param,
      FUN = function(file) {
        path <- file.path(dir, file)
        x <- readRDS(path)
        x <- ProjecTILs.classifier(x, ref.maps[[i]])
        x@meta.data[[paste0(ref.map.name, "_subtypes")]] <- x@meta.data[["functional.cluster"]]
        x@meta.data[["functional.cluster"]] <- NULL
        saveRDS(x, path)
      }
    )
    print(paste("Finished ProjecTILs with map:  ", ref.map.name))
  }
}



#' Calculate cell type composition
#'
#' @param object A seurat object
#' @param sample.col The Seurat object metadata column containing sample names
#' @param annot.cols The Seurat object metadata column(s) containing celltype annotations (provide as character vector, containing the metadata column name(s))
#' @param min.cells Set a minimum threshold for number of cells to calculate relative abundance (e.g. less than 10 cells -> no relative abundnace will be calculated)
#' @param useNA Whether to include not annotated cells or not (labelled as "NA" in the annot.cols). Can be defined separately for each annot.cols (provide single boolean or vector of booleans)
#' @param rename.Multi.to.NA Whether to rename cells labelled as "Multi" by scGate to "NA"
#'
#' @return Cell type compositions as a list of data.frames containing cell counts, relative abundance (freq) and clr-transformed freq (freq_clr), respectively.
#' @export calc_CTcomp
#'
#' @examples
#' devtools::install_github('satijalab/seurat-data')
#'                          library(SeuratData)
#'                          options(timeout = max(300, getOption("timeout")))
#'                          InstallData("panc8")
#'                          data("panc8")
#'                          panc8 = UpdateSeuratObject(object = panc8)
#' # Calculate overall composition
#' celltype.compositions.overall <- calc_CTcomp(object = panc8, annot.cols = "celltype")
#'
#' # Calculate sample-wise composition
#' celltype.compositions.sample_wise <- calc_CTcomp(object = panc8, annot.cols = "celltype", sample.col = "orig.ident")
#'
#' # Plot overall composition
#' par(mar = c(8.1, 4.1, 4.1, 2.1))
#' barplot(unlist(celltype.compositions.overall[["celltype"]][["freq"]]),
#'         ylab = "Relative abundance (%)",
#'         las = 2)
#' # Plot sample-wise composition
#' df <- celltype.compositions.sample_wise[["celltype"]][["freq"]]
#' dflong <- reshape2::melt(t(df))
#' library(ggplot2)
#' ggplot(dflong, aes(x=Var1, y=value, color=Var1)) +
#'   geom_boxplot() + geom_point() + theme_classic() +
#'   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + labs(x = "", y = "Relative abundance (%)") + NoLegend()
calc_CTcomp <- function(object, sample.col = NULL, annot.cols = "scGate_multi",
                        min.cells = 10, useNA = FALSE, rename.Multi.to.NA = FALSE){
  if (length(useNA) == 1) {
    useNA <- ifelse(useNA == TRUE, "always", "no")
  } else if (length(useNA) == length(annot.cols)) {
    x <- c()
    for (i in useNA) {
      x <- c(x, ifelse(useNA == TRUE, "always", "no"))
    }
    useNA <- x
  } else {stop("useNA has not the same length as annot.cols")}

  celltype.compositions <- list()

  for (i in 1:length(annot.cols)) {
    if (length(object@meta.data[[annot.cols[i]]]) == 0) {
      stop(paste("Annotation metadata column", annot.cols[i],"could not be found"))
    }

    if (rename.Multi.to.NA) {
      object@meta.data[[annot.cols[i]]][object@meta.data[[annot.cols[i]]] == "Multi"] <- NA
    }

    if (is.null(sample.col)) {
      comp_table <- table(object@meta.data[[annot.cols[i]]],
                          useNA = useNA[i])
      comp_table_freq <- prop.table(comp_table+1) * 100 # To get percentage
    } else {
      if (length(object@meta.data[[sample.col]]) == 0) {
        stop(paste("Sample column could not be found"))
      }
      comp_table <- table(object@meta.data[[annot.cols[i]]],
                          object@meta.data[[sample.col]],
                          useNA = useNA[i])
      comp_table_freq <- prop.table(comp_table+1, margin = 2) * 100 # To get percentage
    }

    if (sum(comp_table) < min.cells) {
      warning(paste("There are less than", min.cells,
                    "cells detected. This is too few to calculate a reasonable celltype composition, please check or adjust parameter min.cells."))
      next
    }

    ## clr-transform
    comp_table_clr <- Hotelling::clr(t(comp_table_freq))

    ## Append
    celltype.compositions[[annot.cols[i]]][["cell_counts"]] <- as.data.frame.matrix(t(comp_table))
    celltype.compositions[[annot.cols[i]]][["freq"]] <- as.data.frame.matrix(t(comp_table_freq))
    celltype.compositions[[annot.cols[i]]][["freq_clr"]] <- as.data.frame.matrix(comp_table_clr)
  }
  return(celltype.compositions)
}


#' Save object list to disk, in parallel
#'
#' @param obj.list A list of Seurat objects
#' @param dir File location
#' @param ncores Number of CPU cores to use for parallelization
#' @param progressbar Whether to show a progressbar or not
#'
#' @return A list of Seurat objects saved to disk as separate .rds files
#' @export save_objs
#'
#' @examples
#' save_objs(obj.list, "./output/samples")
save_objs <- function(obj.list, dir, ncores = 6, progressbar = T){
  param <- BiocParallel::MulticoreParam(workers = ncores, progressbar = progressbar)
  invisible(BiocParallel::bplapply(
              X = obj.list,
              BPPARAM =  param,
              FUN = function(x) {
                sample_name <- unique(x$Sample)
                file_name <- file.path(dir, sprintf("%s.rds", sample_name))
                saveRDS(x, file_name)
                }))
}


#' Read all .rds files in a directory and return list of Seurat objects, in parallel
#'
#' @param dir File location
#' @param ncores Number of CPU cores to use for parallelization
#' @param progressbar Whether to show a progressbar or not
#'
#' @return A list of Seurat objects read into R
#' @export read_objs
#'
#' @examples
#' obj.list <- read_objs("./output/samples")
read_objs <- function(dir, ncores = 6, progressbar = T){
  param <- BiocParallel::MulticoreParam(workers = ncores, progressbar = progressbar)
  file_names <- list.files(dir)
  obj.list <- BiocParallel::bplapply(
    X = file_names,
    BPPARAM =  param,
    FUN = function(x) {
      readRDS(file.path(dir, x))
      })
  names(obj.list) <- stringr::str_remove_all(file_names, '.rds')
  return(obj.list)
}
