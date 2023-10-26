

#' Annotate cell types in parallel
#'
#' @param dir The directory where your sample .rds files(seurat objects) are located
#' @param scGate.model The scGate model to use (use get_scGateDB() to get a list of available models)
#' @param ref.maps A named list of the ProjecTILs reference maps to use
#' @param ncores The number of cores to use
#' @param progressbar Whether to show a progressbar or not
#'
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom scGate scGate
#' @importFrom ProjecTILs ProjecTILs.classifier
#' @return No return in R. The input files will be overwritten.
#' @export annotate_cells
#'
#' @examples
#' library(ProjecTILs)
#'
#' path_data <- file.path("~", "Dropbox", "CSI", "Standardized_SingleCell_Datasets", "ZhangY_2022_34653365", "output", "test")
#'
#' # Define scGate model
#' scGate_models_DB <- scGate::get_scGateDB(branch = "master", verbose = T, force_update = TRUE)
#' models.TME <- scGate_models_DB$human$TME_HiRes
#'
#' # Load ProjecTILs reference maps
#' path_ref <- "~/Dropbox/CSI/reference_atlases"
#' ref.maps <- list(CD8 = load.reference.map(file.path(path_ref, "CD8T_human_ref_v1.rds")),
#'                  CD4 = load.reference.map(file.path(path_ref, "CD4T_human_ref_v2.rds")),
#'                  DC = load.reference.map(file.path(path_ref, "DC_human_ref_v1.rds")),
#'                  MoMac = load.reference.map(file.path(path_ref, "MoMac_human_v1.rds")))
#' annotate_cells(dir = path_data,
#'                scGate.model = models.TME,
#'                ref.maps = ref.maps,
#'                ncores = 6)
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
         x <- scGate::scGate(x, model=models.TME)
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
        x <- ProjecTILs::ProjecTILs.classifier(x, ref.maps[[i]])
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
#' @param object A seurat object or a list of seurat objects
#' @param dir Directory containing the sample .rds files
#' @param split.by The Seurat object metadata column containing sample names
#' @param annot.cols The Seurat object metadata column(s) containing celltype annotations (provide as character vector, containing the metadata column name(s))
#' @param min.cells Set a minimum threshold for number of cells to calculate relative abundance (e.g. less than 10 cells -> no relative abundnace will be calculated)
#' @param useNA Whether to include not annotated cells or not (labelled as "NA" in the annot.cols). Can be defined separately for each annot.cols (provide single boolean or vector of booleans)
#' @param rename.Multi.to.NA Whether to rename cells labelled as "Multi" by scGate to "NA"
#' @param ncores The number of cores to use
#' @param progressbar Whether to show a progressbar or not
#'
#' @importFrom Hotelling clr
#' @importFrom parallelly availableCores
#' @importFrom stringr str_remove_all
#' @return Cell type compositions as a list of data.frames containing cell counts, relative abundance (freq) and clr-transformed freq (freq_clr), respectively.
#' @export calc_CTcomp
#'
#' @examples
#' devtools::install_github('satijalab/seurat-data')
#' library(SeuratData)
#' options(timeout = max(300, getOption("timeout")))
#' InstallData("panc8")
#' data("panc8")
#' panc8 = UpdateSeuratObject(object = panc8)
#' # Calculate overall composition
#' celltype.compositions.overall <- calc_CTcomp(object = panc8, annot.cols = "celltype")
#'
#' # Calculate sample-wise composition
#' celltype.compositions.sample_wise <- calc_CTcomp(object = panc8, annot.cols = "celltype", split.by = "orig.ident")
calc_CTcomp <- function(object = NULL,
                        dir = NULL,
                        split.by = NULL,
                        annot.cols = "scGate_multi",
                        min.cells = 10,
                        useNA = FALSE,
                        rename.Multi.to.NA = TRUE,
                        ncores = parallelly::availableCores() - 2, progressbar = T) {
  if (!is.null(object) & !is.null(dir)) {
    stop(paste("Cannot provide object and dir"))
  }

  if (length(annot.cols) == 1) {
    useNA <- ifelse(useNA == TRUE, "ifany", "no")
  } else {
    if (length(useNA) == 1) {
      useNA <- ifelse(useNA == TRUE, "ifany", "no")
      useNA <- rep(useNA, length(annot.cols))
    } else if (length(useNA) == length(annot.cols)) {
      x <- c()
      for (i in useNA) {
        x <- c(x, ifelse(useNA == TRUE, "ifany", "no"))
      }
      useNA <- x
    } else {stop("useNA has not the same length as annot.cols")}
  }

  # Either:
  # - The input is a single Seurat object
  # - The input is a list of Seurat objects
  # - The input is directory containing multiple .rds files (each a single object)
  if (isS4(object)) { # Input is a single Seurat object
    celltype.compositions <- list()
    for (i in 1:length(annot.cols)) {
      if (length(object@meta.data[[annot.cols[i]]]) == 0) {
        stop(paste("Annotation metadata column", annot.cols[i],"could not be found"))
      }

      if (rename.Multi.to.NA) {
        object@meta.data[[annot.cols[i]]][object@meta.data[[annot.cols[i]]] == "Multi"] <- NA
      }

      if (is.null(split.by)) {
        comp_table <- table(object@meta.data[[annot.cols[i]]],
                            useNA = useNA[i])
        comp_table_freq <- prop.table(comp_table+1) * 100 # To get percentage
      } else {
        if (length(object@meta.data[[split.by]]) == 0) {
          stop(paste("Sample column could not be found"))
        }
        comp_table <- table(object@meta.data[[annot.cols[i]]],
                            object@meta.data[[split.by]],
                            useNA = useNA[i])
        comp_table_freq <- prop.table(comp_table+1, margin = 2) * 100 # To get percentage
      }

      if (sum(comp_table) < min.cells) {
        warning(paste("There are less than", min.cells,
                      "cells detected. This is too few to calculate a reasonable celltype composition, If needed, set parameter min.cells = 0."))
        next
      }

      ## clr-transform
      comp_table_clr <- Hotelling::clr(t(comp_table_freq))

      ## Append
      celltype.compositions[[annot.cols[i]]][["cell_counts"]] <- as.data.frame.matrix(t(comp_table))
      celltype.compositions[[annot.cols[i]]][["freq"]] <- as.data.frame.matrix(t(comp_table_freq))
      celltype.compositions[[annot.cols[i]]][["freq_clr"]] <- as.data.frame.matrix(comp_table_clr)
    }
  } else { # Input is a list of Seurat objects or a directory containing multiple .rds files
    if (is.list(object)) {
      object_names <- names(object)
    } else if (!is.null(dir)) {
      files <- list.files(dir)
      files <- files[endsWith(files, '.rds')]
      object_names <- stringr::str_remove_all(files, '.rds')
    }

    # For each object calculate the comp_table (cell_counts) for each celltype (annot.col)
    # This is done per sample in order to load every sample only once, if a directory is provided as input
    for(j in 1:length(object_names)) {
      if (is.list(object)) {
        obj <- object[[j]]
      }
      if (!is.null(dir)) {
        obj <- readRDS(file.path(dir, files[j]))
      }

      if (j == 1) { # For the first sample, create a comp_tables list
        comp_tables <- list()
        for (i in 1:length(annot.cols)) {
          if (rename.Multi.to.NA) {
            obj@meta.data[[annot.cols[i]]][obj@meta.data[[annot.cols[i]]] == "Multi"] <- NA
          }
          if (length(obj@meta.data[[annot.cols[i]]]) == 0) {
            stop(paste("Annotation metadata column", annot.cols[i],"could not be found in", object_names[j]))
          }
          if (is.null(split.by)) {
            comp_tables[[annot.cols[i]]] <- as.data.frame.matrix(t(table(obj@meta.data[[annot.cols[i]]],
                                                                         useNA = useNA[i])))
            rownames(comp_tables[[annot.cols[i]]]) <- object_names[j]
          } else {
            stop(paste("Cannot define split.by if multiple objects are provided"))
          }
        }
      } else { # For subsequent samples, append to the comp_tables list
        for (i in 1:length(annot.cols)) {
          if (rename.Multi.to.NA) {
            obj@meta.data[[annot.cols[i]]][obj@meta.data[[annot.cols[i]]] == "Multi"] <- NA
          }
          if (length(obj@meta.data[[annot.cols[i]]]) == 0) {
            stop(paste("Annotation metadata column", annot.cols[i],"could not be found in", object_names[j]))
          }
          if (rename.Multi.to.NA) {
            obj@meta.data[[annot.cols[i]]][obj@meta.data[[annot.cols[i]]] == "Multi"] <- NA
          }
          if (is.null(split.by)) {
            x <- as.data.frame.matrix(t(table(obj@meta.data[[annot.cols[i]]],
                                              useNA = useNA[i])))
            rownames(x) <- object_names[j]
            comp_tables[[annot.cols[i]]] <- as.data.frame.matrix(t(transform(merge(t(x),t(comp_tables[[annot.cols[i]]]),by=0,all=TRUE), row.names=Row.names, Row.names=NULL)))
          } else {
            stop(paste("Cannot define split.by if multiple objects are provided"))
          }
        }
      }
    }

    celltype.compositions <- list()

    for (i in 1:length(annot.cols)) {
      comp_table <- comp_tables[[i]]
      comp_table[is.na(comp_table)] <- 0

      low_count_samples <- rownames(comp_table)[rowSums(comp_table) < min.cells]

      comp_table_freq <- as.data.frame(prop.table(as.matrix(comp_table)+1, margin = 1) * 100) # To get percentage

      if (length(low_count_samples) >= 1) {
        warning(paste("There are less than", min.cells, annot.cols[i],
                      "cells detected in sample(s)", low_count_samples,
                      ". For this sample(s), no celltype composition was calculated. If needed, set parameter min.cells = 0.\n"))
        comp_table_freq <- comp_table_freq[!rownames(comp_table_freq) %in% low_count_samples, ]
      }

      ## clr-transform
      comp_table_clr <- Hotelling::clr(comp_table_freq)

      ## Sort
      comp_table <- as.data.frame.matrix(comp_table)
      comp_table_clr <- as.data.frame.matrix(comp_table_clr)
      comp_table <- comp_table[order(row.names(comp_table)), ]
      comp_table_freq <- comp_table_freq[order(row.names(comp_table_freq)), ]
      comp_table_clr <- comp_table_clr[order(row.names(comp_table_clr)), ]

      ## Append
      celltype.compositions[[annot.cols[i]]][["cell_counts"]] <- comp_table
      celltype.compositions[[annot.cols[i]]][["freq"]] <- comp_table_freq
      celltype.compositions[[annot.cols[i]]][["freq_clr"]] <- comp_table_clr
    }
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
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom parallelly availableCores
#' @return A list of Seurat objects saved to disk as separate .rds files
#' @export save_objs
#'
#' @examples
#' save_objs(obj.list, "./output/samples")
save_objs <- function(obj.list,
                      dir,
                      ncores = parallelly::availableCores() - 2, progressbar = T){
  BiocParallel::bplapply(
    X = obj.list,
    BPPARAM =  BiocParallel::MulticoreParam(workers = ncores, progressbar = progressbar),
    FUN = function(x) {
      sample_name <- unique(x$Sample)
      file_name <- file.path(dir, sprintf("%s.rds", sample_name))
      saveRDS(x, file_name)
      })
}


#' Read all .rds files in a directory and return list of Seurat objects, in parallel
#'
#' @param dir File location
#' @param file.list A list of files (with full pathname)
#' @param ncores Number of CPU cores to use for parallelization
#' @param progressbar Whether to show a progressbar or not
#'
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom stringr str_remove_all
#' @importFrom parallelly availableCores
#' @return A list of Seurat objects read into R
#' @export read_objs
#'
#' @examples
#' obj.list <- read_objs("./output/samples")
read_objs <- function(dir = NULL,
                      file.list = NULL,
                      ncores = parallelly::availableCores() - 2, progressbar = T){
  if (!is.null(dir) & is.null(file.list)) {
    file_names <- list.files(dir)
    file_paths <- file.path(dir, file_names)
  } else if (is.null(dir) & !is.null(file.list)){
    file_paths <- file.list[endsWith(files, '.rds')]
    file_names <- gsub("^.*/", "", file_paths)
  }
  obj.list <- BiocParallel::bplapply(
    X = file_paths,
    BPPARAM =  BiocParallel::MulticoreParam(workers = ncores, progressbar = progressbar),
    FUN = function(x) {
      readRDS(file.path(x))
      })
  names(obj.list) <- stringr::str_remove_all(file_names, '.rds')
  return(obj.list)
}
