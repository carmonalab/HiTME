

#' Run HitME: classify cells
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
#' @param bparam A [BiocParallel::bpparam()] object that tells Run.HiTME how to parallelize. If provided, it overrides the `ncores` parameter.
#' @param progressbar Whether to show a progressbar or not
#'
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom parallelly availableCores
#' @importFrom scGate scGate
#' @importFrom SignatuR SignatuR
#' @importFrom ProjecTILs ProjecTILs.classifier
#' @importFrom Seurat SplitObject
#'
#' @return No return in R. The input files will be overwritten.
#' @export Run.HiTME
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
#' Run.HiTME(dir = path_data,
#'                scGate.model = models.TME,
#'                ref.maps = ref.maps,
#'                ncores = 6)
#'

Run.HiTME <- function(object = NULL,
                       dir = NULL,
                       scGate.model = "default",
                       scGate.model.branch = c("master", "dev"),
                       ref.maps = NULL,
                       split.by = NULL,
                       group.by = c("layer1","layer2"),
                       remerge = TRUE,
                       species = "human",
                       additional.signatures = "default",
                       return.Seurat = FALSE,
                       bparam = NULL,
                       ncores = parallelly::availableCores() - 2,
                       progressbar = TRUE,
                       ...){

  if (is.null(object) & is.null(dir)) {
    stop("Please provide either a Seurat object or a directory with rds files")
  }

  if (!is.null(object) & !is.null(dir)) {
    stop("Cannot provide both object and directory")
  }

  if (!is.null(ref.maps)) {
    if (!is.list(ref.maps)) {
      stop("Please provide ref.maps as named list, containing the reference map(s), even if it is just one. The name will be used for the metadata column containing the cell type annotations")
    } else if(!all(unlist(lapply(ref.maps, function(x){class(x) == "Seurat"})))) {
      stop("One or more reference maps are not a Seurat object")
    }
  }

  # split object into a list if indicated
  if (!is.null(split.by)) {
    if(!split.by %in% names(object@meta.data)){
      stop(paste(split.by, "is not a metadata column in this Seurat object"))
    }
    object <- Seurat::SplitObject(object, split.by = split.by)

  }

  if (!is.null(dir)) {
    files <- list.files(dir,
                        pattern = "[.]rds$")

    # when running from directory only return seurat object
    return.Seurat <- TRUE
  }

  # adapt species
  if(is.null(species)){
    stop("Please provide human or mouse as species")
  }
  species <- tolower(species)
  if(!species == "human"){
    if(grepl("homo|sap|huma", species)){
      species <- "human"
    }
  } else if (!species == "mouse"){
    if(grepl("mice|mus", species)){
      species <- "mouse"
    }
  } else {
    stop("Only supported species are human and mouse")
  }

  # if object is unique turn into a list
  if(!is.list(object)){
    object <- list(object)
  }

  if(ncores >= parallelly::availableCores()){
    ncores <- parallelly::availableCores() - 1
    message("Using all or more cores available in this computer, reducing number of cores to ", ncores)
  }

  # Reduce number of cores if file is huge
  if(sum(unlist(lapply(object, ncol)))>=30000){
    ncores <- 2
    message("Huge Seurat object, reducing number of cores to avoid memory issues to", ncores)
  }


  # set paralelization parameters
  if (is.null(bparam)) {
    if (ncores>1) {
      param <- BiocParallel::MulticoreParam(workers =  ncores, progressbar = progressbar)
    } else {
      param <- SerialParam()
    }
  } else {
    param <- bparam
  }

  # Either:
  # - The input is a single Seurat object
  # - The input is a list of Seurat objects
  # - The input is directory containing multiple .rds files (each a single object)

  # Run scGate, if model is provided
  if (!is.null(scGate.model)) {
    message("Running scGate\n")


    # Retrieve default scGate models if default
    if(!is.list(scGate.model) && tolower(scGate.model) == "default"){
      if(species == "human"){
        scGate.model <- scGate::get_scGateDB(branch = scGate.model.branch,
                                             force_update = T)[[species]][["TME_HiRes"]]
      } else if(species == "mouse"){
        scGate.model <- scGate::get_scGateDB(branch = scGate.model.branch,
                                             force_update = T)[[species]][["generic"]]
      }
      message(" - Running scGate model for ", paste(names(scGate.model), collapse = ", "), "\n")
    }

    # based on species get SignatuR signatures
    if(!is.list(additional.signatures) && tolower(additional.signatures) == "default"){
      if(species == "human"){
        sig.species <- "Hs"
      } else if(species == "mouse"){
        sig.species <- "Mm"
      }
      additional.signatures <- SignatuR::GetSignature(SignatuR[[sig.species]][["Programs"]])
      message(" - Adding additional signatures: ", paste(names(additional.signatures), collapse = ", "), "\n")
    }

  ## Run scGate
  if (is.list(object)) {
      object <- lapply(
        X = object,
        FUN = function(x) {
          if(class(x) != "Seurat"){
            stop("Not Seurat object included, cannot be processed.\n")
          }
          x <- scGate::scGate(x,
                              model=scGate.model,
                              additional.signatures = additional.signatures,
                              BPPARAM = param)
          # convert multi to NA
          x@meta.data <-
            x@meta.data %>% dplyr::mutate(scGate_multi = ifelse(scGate_multi == "Multi",
                                                            NA, scGate_multi))

          x@misc[["layer1_param"]] <- list()
          x@misc[["layer1_param"]][["layer1_models"]] <- names(scGate.model)
          x@misc[["layer1_param"]][["additional.signatures"]] <- names(additional.signatures)
          return(x)
        }
      )
    } else if (!is.null(dir)) { # If input is directory
      lapply(
        X = files,
        FUN = function(file) {
          path <- file.path(dir, file)
          x <- readRDS(path)
          if(class(x) != "Seurat"){
            stop("Not Seurat object included, cannot be processed.\n")
          }
          x <- scGate::scGate(x,
                              model= scGate.model,
                              additional.signatures = additional.signatures,
                              BPPARAM = param)

          # convert multi to NA
          x@meta.data <-
            x@meta.data %>% dplyr::mutate(scGate_multi = ifelse(scGate_multi == "Multi",
                                                            NA, scGate_multi))

          x@misc[["layer1_param"]] <- list()
          x@misc[["layer1_param"]][["layer1_models"]] <- names(scGate.model)
          x@misc[["layer1_param"]][["additional.signatures"]] <- names(additional.signatures)
          return(x)

          saveRDS(x, path)
        }
      )
    }
    # close paralel cores
    BiocParallel::bpstop(param)



    message("Finished scGate\n####################################################\n")
  } else {
    message("Not running coarse cell type classification as no scGate model was indicated.\n")
  }


  # Run ProjecTILs if ref.maps is provided
  if (!is.null(ref.maps)) {
    # Run each map in paralel

      if (is.list(object)) {
        object <- lapply(
          X = object,
          FUN = function(x) {
            if(class(x) != "Seurat"){
              stop("Not Seurat object included, cannot be processed.\n")
            }
            x <- ProjecTILs.classifier.multi(x,
                                             ref.maps = ref.maps,
                                             bparam = param)
            return(x)
          }
        )

        # If input is directory
      } else if (!is.null(dir)) {
        lapply(
          X = files,
          FUN = function(file) {
            path <- file.path(dir, file)
            x <- readRDS(path)
            if(class(x) != "Seurat"){
              stop("Not Seurat object included, cannot be processed.\n")
            }
            x <- ProjecTILs.classifier.multi(x,
                                             ref.maps = ref.maps,
                                             bparam = param)
            saveRDS(x, path)
          }
        )
      }
    # close paralel cores
    BiocParallel::bpstop(param)
    message("Finished Projectils\n####################################################\n")
  } else {
    message("Not running reference mapping as no reference maps were indicated.\n")
  }




  if (is.list(object)) {
    if (remerge && length(object)>1) {
        object <- merge(object[[1]],object[2:length(object)])
        object <- list(object)
    }
  } else {
    object <- list(object)
  }

  if(!return.Seurat){
      hit <- lapply(object, get.HiTObject,
                            group.by = group.by)
  }

  # if list is of 1, return object not list

  if (!is.null(object)) {
    if(return.Seurat){
      if(length(object)==1){
        object <- object[[1]]
      }
      return(object)
    } else {
      if(length(hit)==1){
        hit <- hit[[1]]
      }
      return(hit)
    }
  }

}



# Function get.HiTObject
#' Get a HiT object summarizing cell type classification and aggregating profiles.
#'
#'
#'
#' @param object A seurat object
#' @param group.by List with one or multiple Seurat object metadata columns with cell type predictions to group by (e.g. layer 1 cell type classification)
#' @param name.additional.signatures Names of additional signatures as found in object metadata to take into account
#' @param ... Additional parameters for \link{get.aggregated.profile}, \link{get.celltype.composition}, and \link{get.aggregated.signature} functions.
#' @import Seurat
#' @import SeuratObject
#' @return HiT object summarizing cell type classification and aggregating profiles.
#' @export get.HiTObject
#'


get.HiTObject <- function(object,
                          group.by = list("layer1" = c("scGate_multi"),
                                          "layer2" = c("funciontal.cluster")
                                          ),
                          name.additional.signatures = NULL,
                          ...){


  if (is.null(object)) {
    stop("Please provide either a Seurat object")
  } else if(class(object) != "Seurat"){
    stop("Not Seurat object included, cannot be processed.\n")
  }

  for(v in seq_along(group.by)){
    if(is.null(names(group.by[v]))){
      names(group.by[v]) <- paste0("layer", v)
    }
  }

  # if only group.by is indicated use the same for aggregation and composition
  if(!exists("group.by.aggregated")){
    group.by.aggregated <- group.by
  }

  if(!exists("group.by.composition")){
    group.by.composition <- group.by
  }

  if(is.null(c(group.by.aggregated, group.by.composition))){
    stop("Please provide at least one grouping variable")
  }

  ## Build S4 objects to store data
  setClass(Class = "HiT",
           slots = list(
             metadata = "data.frame",
             predictions = "list",
             composition = "list",
             aggregated_profile = "list",
             version = "list"
           )
  )

    # make list of list for layers for predictions slot
  pred.list <- list()
  for(a in names(group.by)){
    pred.list[[a]] <- list()
    pred.list[[a]] <- object@meta.data[,group.by[[a]], drop = F]
  }



  #run extended if object got scGate info
  # extract values from misc slot from object
  if(any(group.by == "scGate_multi")){
    layer.scgate <- names(group.by[group.by == "scGate_multi"])
    if(is.null(name.additional.signatures)){
      name.additional.signatures <- object@misc$layer1_param$additional.signatures
    }
  scgate.models <- object@misc$layer1_param$layer1_models


  if(!is.null(name.additional.signatures)){
  sig <- grep(paste(additional.signatures, collapse = "|"),
              names(object@meta.data), value = T)
  sig.df <- object@meta.data[, sig]
  pred.list[[layer.scgate]][["additional_signatures"]] <- sig.df
  }

  if(!is.null(scgate.models)){
  scgate <- grep(paste(scgate.models, collapse = "$|"),
                 names(object@meta.data), value = T)
  scgate.df <- object@meta.data[, scgate]
  pred.list[[layer.scgate]][["scGate_is.pure"]] <- scgate.df
  }

  if(!is.null(name.additional.signatures) && !is.null(scgate.models)){
  ucell <- names(object@meta.data)[!names(object@meta.data) %in% c(sig, scgate)] %>%
            grep("_UCell$", ., value = T)

  ucell.df <- object@meta.data[, ucell]

  pred.list[[layer.scgate]][["UCell_scores"]] <- ucell.df
  }
  }

  #run extended if object got projectils info
  # extract values from misc slot from object
  if(any(group.by == "functional.cluster")){
    layer.pt <- names(group.by[group.by == "functional.cluster"])
    pred.list[[layer.pt]][["functional.cluster"]] <- object@meta.data[,"functional.cluster", drop = F]
    pred.list[[layer.pt]][["functional.cluster.conf"]] <- object@meta.data[,"functional.cluster.conf", drop = F]

  }


  # Compute proportions
  comp.prop <- get.celltype.composition(object,
                                        group.by.composition = group.by.composition,
                                        ...)
  # Compute avg expression
  avg.expr <- get.aggregated.profile(object,
                                     group.by.aggregated = group.by.aggregated,
                                     ...)

  aggr.signature <- get.aggregated.signature(object,
                                             group.by.aggregated = group.by.aggregated,
                                             ...)


  hit <- new("HiT",
             metadata = object@meta.data,
             predictions = pred.list,
             aggregated_profile = list("Gene expression" = avg.expr,
                                       "Signatures" = aggr.signature),
             composition = comp.prop
  )


  return(hit)


}


#' Calculate cell type composition
#'
#' @param object A seurat object or a list of seurat objects
#' @param dir Directory containing the sample .rds files
#' @param split.by A Seurat object metadata column to split by (e.g. sample names)
#' @param group.by.composition The Seurat object metadata column(s) containing celltype annotations (provide as character vector, containing the metadata column name(s))
#' @param min.cells Set a minimum threshold for number of cells to calculate relative abundance (e.g. less than 10 cells -> no relative abundnace will be calculated)
#' @param useNA Whether to include not annotated cells or not (labelled as "NA" in the group.by.composition). Can be defined separately for each group.by.composition (provide single boolean or vector of booleans)
#' @param clr_zero_impute_perc To calculate the clr-transformed relative abundance ("clr_freq"), zero values are not allowed and need to be imputed (e.g. by adding a pseudo cell count). Instead of adding a pseudo cell count of flat +1, here a pseudo cell count of +1% of the total cell count will be added to all cell types, to better take into consideration the relative abundance ratios (e.g. adding +1 cell to a total cell count of 10 cells would have a different, i.e. much larger effect, than adding +1 to 1000 cells).
#' @param ncores The number of cores to use
#' @param progressbar Whether to show a progressbar or not
#'
#' @importFrom Hotelling clr
#' @importFrom parallelly availableCores
#' @importFrom stringr str_remove_all
#' @return Cell type compositions as a list of data.frames containing cell counts, relative abundance (freq) and clr-transformed freq (freq_clr), respectively.
#' @export get.celltype.composition
#'
#' @examples
#' devtools::install_github('satijalab/seurat-data')
#' library(SeuratData)
#' options(timeout = max(300, getOption("timeout")))
#' InstallData("panc8")
#' data("panc8")
#' panc8 = UpdateSeuratObject(object = panc8)
#' # Calculate overall composition
#' celltype.compositions.overall <- get.celltype.composition(object = panc8, group.by.composition = "celltype")
#'
#' # Calculate sample-wise composition
#' celltype.compositions.sample_wise <- get.celltype.composition(object = panc8, group.by.composition = "celltype", split.by = "orig.ident")
get.celltype.composition <- function(object = NULL,
                                     dir = NULL,
                                     split.by = NULL,
                                     group.by.composition = NULL,
                                     min.cells = 10,
                                     useNA = FALSE,
                                     clr_zero_impute_perc = 1) {

  if (!is.null(object) & !is.null(dir)) {
    stop(paste("Cannot provide object and dir"))
  }

  if (length(group.by.composition) == 1) {
    useNA <- ifelse(useNA == TRUE, "ifany", "no")
  } else {
    if (length(useNA) == 1) {
      useNA <- ifelse(useNA == TRUE, "ifany", "no")
      useNA <- rep(useNA, length(group.by.composition))
    } else if (length(useNA) == length(group.by.composition)) {
      x <- c()
      for (i in useNA) {
        x <- c(x, ifelse(useNA == TRUE, "ifany", "no"))
      }
      useNA <- x
    } else {stop("useNA has not the same length as group.by.composition")}
  }

  celltype.compositions <- list()

  # Either:
  # - The input is a single Seurat object
  # - The input is a list of Seurat objects
  # - The input is directory containing multiple .rds files (each a single object)
  if (isS4(object)) { # Input is a single Seurat object
    for (i in 1:length(group.by.composition)) {
      if (length(object@meta.data[[group.by.composition[[i]]]]) == 0) {
        stop(paste("Annotation metadata column", group.by.composition[[i]],"could not be found"))
      }

      if (is.null(split.by)) {
        comp_table <- table(object@meta.data[[group.by.composition[[i]]]],
                            useNA = useNA[i])
        comp_table_freq <- prop.table(comp_table) * 100 # To get percentage
      } else {

        if (length(object@meta.data[[split.by]]) == 0) {
          stop(paste("Sample column could not be found. If providing multiple columns, indicate a vector."))
        }

        comp_table <- table(object@meta.data[[group.by.composition[[i]]]],
                            object@meta.data[[split.by]],
                            useNA = useNA[i])

        comp_table_freq <- prop.table(comp_table, margin = 2) * 100 # To get percentage
      }

      if (sum(comp_table) < min.cells) {
        warning(paste("There are less than", min.cells,
                      "cells detected. This is too few to calculate a reasonable celltype composition, If needed, set parameter min.cells = 0."))
        next
      }

      ## clr-transform
      comp_table_clr <- Hotelling::clr(t(comp_table_freq + clr_zero_impute_perc))

      ## Append
      celltype.compositions[[group.by.composition[[i]]]][["cell_counts"]] <- as.data.frame.matrix(t(comp_table))
      celltype.compositions[[group.by.composition[[i]]]][["freq"]] <- as.data.frame.matrix(t(comp_table_freq))
      celltype.compositions[[group.by.composition[[i]]]][["freq_clr"]] <- as.data.frame.matrix(comp_table_clr)
    }
  } else { # Input is a list of Seurat objects or a directory containing multiple .rds files
    comp_tables <- list()

    if (is.list(object)) {
      object_names <- names(object)
    } else if (!is.null(dir)) {
      files <- list.files(dir, pattern = "[.]rds$")
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

      for (i in 1:length(group.by.composition)) {
        if (rename.Multi.to.NA) {
          obj@meta.data[[group.by.composition[[i]]]][obj@meta.data[[group.by.composition[[i]]]] == "Multi"] <- NA
        }
        if (length(obj@meta.data[[group.by.composition[[i]]]]) == 0) {
          stop(paste("Annotation metadata column", group.by.composition[[i]],"could not be found in", object_names[j]))
        }
        if (rename.Multi.to.NA) {
          obj@meta.data[[group.by.composition[[i]]]][obj@meta.data[[group.by.composition[[i]]]] == "Multi"] <- NA
        }
        if (!is.null(split.by)) {
          stop(paste("Cannot define split.by if multiple objects are provided"))
        } else {
          if (j == 1) { # For the first sample, simply add table to list
            comp_tables[[group.by.composition[[i]]]] <- as.data.frame.matrix(t(table(obj@meta.data[[group.by.composition[[i]]]],
                                                                                   useNA = useNA[i])))
            rownames(comp_tables[[group.by.composition[[i]]]]) <- object_names[j]
          } else { # For subsequent samples, need to merge tables
            comp_table_to_append <- as.data.frame.matrix(t(table(obj@meta.data[[group.by.composition[[i]]]],
                                                                 useNA = useNA[i])))
            rownames(comp_table_to_append) <- object_names[j]
            comp_tables_merged <- merge(t(comp_table_to_append),
                                        t(comp_tables[[group.by.composition[[i]]]]),
                                        by=0, all=TRUE)
            comp_tables[[group.by.composition[[i]]]] <- as.data.frame.matrix(t(transform(comp_tables_merged, row.names=Row.names, Row.names=NULL)))
          }
        }
      }
    }

    for (i in 1:length(group.by.composition)) {
      comp_table <- comp_tables[[i]]
      comp_table[is.na(comp_table)] <- 0

      low_count_samples <- rownames(comp_table)[rowSums(comp_table) < min.cells]

      comp_table_freq <- as.data.frame(prop.table(as.matrix(comp_table), margin = 1) * 100) # To get percentage

      if (length(low_count_samples) >= 1) {
        warning(paste("There are less than", min.cells, group.by.composition[[i]],
                      "cells detected in sample(s)", low_count_samples,
                      ". For this sample(s), no celltype composition was calculated. If needed, set parameter min.cells = 0.\n"))
        comp_table_freq <- comp_table_freq[!rownames(comp_table_freq) %in% low_count_samples, ]
      }

      ## clr-transform
      comp_table_clr <- as.data.frame.matrix(Hotelling::clr(comp_table_freq + clr_zero_impute_perc))

      ## Sort
      comp_table <- comp_table[order(row.names(comp_table)), ]
      comp_table_freq <- comp_table_freq[order(row.names(comp_table_freq)), ]
      comp_table_clr <- comp_table_clr[order(row.names(comp_table_clr)), ]

      ## Append
      celltype.compositions[[group.by.composition[[i]]]][["cell_counts"]] <- comp_table
      celltype.compositions[[group.by.composition[[i]]]][["freq"]] <- comp_table_freq
      celltype.compositions[[group.by.composition[[i]]]][["freq_clr"]] <- comp_table_clr
    }
  }
  return(celltype.compositions)
}

### Function to compute average expression

get.aggregated.profile <- function(object,
                                   group.by.aggregated = NULL,
                                   gene.filter = NULL,
                                   GO_accession = NULL,
                                   assay = "RNA",
                                   slot = "data") {
  if (is.null(object)) {
    stop("Please provide a Seurat object")
    if(class(object) != "Seurat"){
      stop("Please provide a Seurat object")
    }
  }

  if(is.null(group.by.aggregated)){
    message("No grouping provided, grouping by consensus of Projectils\n")
    group.by.aggregated <- "functional.cluster"
  } else if (!is.vector(group.by.aggregated)) {
    stop("Please provide one or various grouping variables as a vector")
  } else {
    group.by.aggregated <- unique(c(group.by.aggregated, "functional.cluster"))
  }

  if(any(!group.by.aggregated %in% names(object@meta.data))){
    g.missing <- group.by.aggregated[!group.by.aggregated %in% names(object@meta.data)]
    warning(paste(g.missing, "not in object metadata, it is not computed for average expression.\n"))
    group.by.aggregated <- group.by.aggregated[group.by.aggregated %in% names(object@meta.data)]
  }

  if(!is.null(gene.filter)){
    if(!is.list(gene.filter)){
      stop("Please add additional subsetting list of genes as a named list format")
    }}


  gene.filter.list <- list()

  # Dorothea transcription factors
  # data("entire_database", package = "dorothea")
  # gene.filter.list[["Dorothea_Transcription_Factors"]] <- entire_database$tf %>% unique()

  # Ribosomal genes
  gene.filter.list[["Ribosomal"]] <- SignatuR::GetSignature(SignatuR$Hs$Compartments$Ribo)

  # Add list genes from GEO accessions
  gene.filter.list <- c(gene.filter.list, GetGOList(GO_accession))

  # Add defined list of genes to filter
  for(a in seq_along(gene.filter)){
    if(is.null(names(gene.filter[[a]]))){
      i.name <- paste0("GeneList_", a)
    } else {
      i.name <- names(gene.filter[[a]])
    }
    gene.filter.list[[i.name]] <- gene.filter[[a]]
  }


  avg.exp <- list()

  # loop over different grouping

  for(i in group.by.aggregated){
    avg.exp[[i]] <- list()
    # only return avg expression for existing groups

    all.genes  <-
      Seurat::AverageExpression(object,
                                group.by = i,
                                assays = assay,
                                slot = slot)[[assay]]

    avg.exp[[i]][["All.genes"]] <- all.genes

    for(e in names(gene.filter.list)){
      keep <- gene.filter.list[[e]][gene.filter.list[[e]] %in% rownames(all.genes)]
      avg.exp[[i]][[e]] <- all.genes[keep,]
    }


  }

  return(avg.exp)
}





get.aggregated.signature <- function(object,
                                     group.by.aggregated = NULL,
                                     name.additional.signatures = NULL){

  if(is.null(group.by.aggregated)){
    stop("Please provide groping variable for aggregation")
  }

  if(is.null(name.additional.signatures)){
    name.additional.signatures <- object@misc$layer1_param$additional.signatures
  }

  if(is.null(name.additional.signatures)){
    stop("Please provide additional signatures")
  }

  if(!any(grepl(paste(name.additional.signatures, collapse = "|"),
                names(object@meta.data)))){
    stop("No additional signatures found in this object metadata")
  }

  add.sig.cols <- grep(paste(name.additional.signatures, collapse = "|"),
                       names(object@meta.data), value = T)

  aggr.sig <- list()

  for(e in group.by.aggregated){
    aggr.sig[[e]] <- object@meta.data %>%
      group_by(.data[[e]]) %>%
      summarize_at(add.sig.cols, mean, na.rm = T)
  }

  return(aggr.sig)

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




