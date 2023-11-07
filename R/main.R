

#' Annotate cell types in parallel
#'
#' @param object A seurat object or a list of seurat objects
#' @param dir The directory where your sample .rds files(seurat objects) are located
#' @param scGate.model The scGate model to use (use get_scGateDB() to get a list of available models)
#' @param ref.maps A named list of the ProjecTILs reference maps to use
#' @param split.by A Seurat object metadata column to split by (e.g. sample names)
#' @param remerge When setting split.by, if remerge = TRUE one object will be returned. If remerge = FALSE a list of objects will be returned.
#' @param additional.signatures Adding UCell additional signatures to compute on each cell.
#' @param return.Seurat Whether to return a Hit S4 object or add the data into Seurat object metadata
#' @param ncores The number of cores to use
#' @param progressbar Whether to show a progressbar or not
#'
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom parallelly availableCores
#' @importFrom scGate scGate
#' @importFrom ProjecTILs ProjecTILs.classifier
#' @importFrom Seurat SplitObject
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
#'

annotate_cells <- function(object = NULL,
                           dir = NULL,
                           scGate.model = NULL,
                           ref.maps = NULL,
                           split.by = NULL,
                           remerge = TRUE,
                           additional.signatures = NULL,
                           return.Seurat = FALSE,
                           ncores = parallelly::availableCores() - 2,
                           progressbar = TRUE
                           ){

  if (is.null(object) & is.null(dir)) {
    stop("Please provide either a Seurat object or a directory with rds files")
  }

  if (!is.null(object) & !is.null(dir)) {
    stop("Cannot provide both object and directory")
  }

  if (!is.null(ref.maps)) {
    if (!is.list(ref.maps)) {
      stop("Please provide ref.maps as named list, containing the reference map(s), even if it is just one. The name will be used for the metadata column containing the cell type annotations")
    }
  }

  if (is.null(ref.maps) & is.null(scGate.model)) {
    stop("Please provide at least either an scGate model or reference map for cell type classification")
  }

  # split object into a list if indicated
  if (!is.null(split.by)) {
    object <- Seurat::SplitObject(object, split.by = split.by)

  }

  if (!is.null(dir)) {
    files <- list.files(dir,
                        pattern = "[.]rds$")
    # get ncores per file
    internal_cores <- floor(ncores / length(files))
    if(internal_cores == 0){internal_cores <- 1}

    # when running from directory only return seurat object
    return.Seurat <- TRUE

  } else if (!is.null(object)){
    # Define internal ncores for scGate and projectils functions
    if(is.list(object)){
      internal_cores <- floor(ncores / length(object))
      if(internal_cores == 0){internal_cores <- 1}
    } else {
      internal_cores <- ncores
    }
  }





  # Either:
  # - The input is a single Seurat object
  # - The input is a list of Seurat objects
  # - The input is directory containing multiple .rds files (each a single object)

  param <- BiocParallel::MulticoreParam(workers = ncores, progressbar = progressbar)


  # Run scGate, if model is provided
  if (!is.null(scGate.model)) {
    message("Running scGate\n")
    # If input is a single Seurat object
    if (isS4(object)) {
      if(class(object) != "Seurat"){
        stop("Not Seurat object included, cannot be processed.\n")
      }
      object <- scGate::scGate(object,
                               model = scGate.model,
                               additional.signatures = additional.signatures,
                               ncores = internal_cores)
    # If input is object list
    } else if (is.list(object)) {
      object <- BiocParallel::bplapply(
        X = object,
        BPPARAM = param,
        FUN = function(x) {
          if(class(x) != "Seurat"){
            stop("Not Seurat object included, cannot be processed.\n")
          }
          x <- scGate::scGate(x,
                              model=scGate.model,
                              additional.signatures = additional.signatures,
                              ncores = internal_cores)
        }
      )
    } else if (!is.null(dir)) { # If input is directory
      BiocParallel::bplapply(
        X = files,
        BPPARAM = param,
        FUN = function(file) {
          path <- file.path(dir, file)
          x <- readRDS(path)
          if(class(x) != "Seurat"){
            stop("Not Seurat object included, cannot be processed.\n")
          }
          x <- scGate::scGate(x, model=models.TME, ncores = 1)
            x <- scGate::scGate(x,
                                model= scGate.model,
                                additional.signatures = additional.signatures,
                                ncores = internal_cores)
          saveRDS(x, path)
        }
      )
    }



    message("Finished scGate\n")
  } else {
    message("Not running coarse cell type classification as no scGate model was indicated.\n")
  }


  # Run ProjecTILs if ref.maps is provided
  if (!is.null(ref.maps)) {
    for (i in 1:length(ref.maps)) {
      ref.map.name <- names(ref.maps)[i]
      message(paste("Running ProjecTILs with map:  ", ref.map.name, "\n"))

      # If input is a single Seurat object
      if (isS4(object)) {
        if(class(object) != "Seurat"){
          stop("Not Seurat object included, cannot be processed.\n")
        }
        object <- ProjecTILs::ProjecTILs.classifier(object,
                                                    ref.maps[[i]],
                                                    ncores = internal_cores)
        object <- addProjectilsClassification(object,
                                              ref.map.name = ref.map.name)

        # If input is object list
      } else if (is.list(object)) {
        BiocParallel::bplapply(
          X = object,
          BPPARAM = param,
          FUN = function(x) {
            if(class(x) != "Seurat"){
              stop("Not Seurat object included, cannot be processed.\n")
            }
            x <- ProjecTILs::ProjecTILs.classifier(x,
                                                   ref.maps[[i]],
                                                   ncores = internal_cores)
            x <- addProjectilsClassification(x, ref.map.name = ref.map.name)
            saveRDS(x, path)
          }
        )

        # If input is directory
      } else if (!is.null(dir)) {
        BiocParallel::bplapply(
          X = files,
          BPPARAM = param,
          FUN = function(file) {
            path <- file.path(dir, file)
            x <- readRDS(path)
            if(class(x) != "Seurat"){
              stop("Not Seurat object included, cannot be processed.\n")
            }
            x <- ProjecTILs::ProjecTILs.classifier(x,
                                                   ref.maps[[i]],
                                                   ncores = internal_cores)
            x <- addProjectilsClassification(x, ref.map.name = ref.map.name)
            saveRDS(x, path)
          }
        )
      }
      print(paste("Finished ProjecTILs with map:  ", ref.map.name))
    }
  } else {
    message("Not running reference mapping as no reference maps were indicated.\n")
  }


  if (!is.null(split.by)) {
    if (remerge) {
      object <- merge(object[[1]],object[2:length(object)])

      # Add consensus projectils annotation
      object <- classificationConsensus(object)
      }
    } else {
      object <- lapply(object, classificationConsensus)
    }


  if (!is.null(object)) {
    if(return.Seurat){
      return(object)
    } else {

      return(hit)
    }
  }
}





#' Calculate cell type composition
#'
#' @param object A seurat object or a list of seurat objects with HiTME classification on its metadata
#' @param dir Directory containing the sample .rds files
#' @param split.by A Seurat object metadata column to split by (e.g. sample names)
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
                        ncores = parallelly::availableCores() - 2,
                        progressbar = T) {

  if (is.null(object) & is.null(dir)) {
    stop("Please provide either a Seurat object or a directory with rds files")
  }

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
        stop(paste("Annotation metadata column", annot.cols[i],"could not be found in this Seurat object"))
      }

      if (rename.Multi.to.NA) {
        object@meta.data[[annot.cols[i]]][object@meta.data[[annot.cols[i]]] == "Multi"] <- NA
      }

      if (is.null(split.by)) {
        comp_table <- table(object@meta.data[[annot.cols[i]]],
                            useNA = useNA[i])
        comp_table_freq <- prop.table(comp_table) * 100 # To get percentage
        freq_clr <- prop.table(comp_table + 1) * 100
      } else {
        if (length(object@meta.data[[split.by]]) == 0) {
          stop(paste("Split.by column could not be found"))
        }
        comp_table <- table(object@meta.data[[annot.cols[i]]],
                            object@meta.data[[split.by]],
                            useNA = useNA[i])
        comp_table_freq <- prop.table(comp_table, margin = 2) * 100 # To get percentage
        freq_clr <- prop.table(comp_table + 1, margin = 2) * 100
      }

      if (sum(comp_table) < min.cells) {
        warning(paste("There are less than", min.cells,
                      "cells detected. This is too few to calculate a reasonable celltype composition, If needed, set parameter min.cells = 0."))
        next
      }

      ## clr-transform
      comp_table_clr <- Hotelling::clr(t(freq_clr))

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
