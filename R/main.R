#' Classify cells using scGate and ProjecTILs.
#'
#' @param object A seurat object or a list of seurat objects
#' @param scGate.model The scGate model to use (use get_scGateDB() to get a list of available models), by default fetch the HiTME models: \code{ scGate::get_scGateDB(branch = scGate.model.branch)[[species]][["HiTME"]]}.
#' @param scGate.model.branch From which branch Run.HiTME fetch the scGate models, by default models are retrieved from \code{master} branch.
#' @param multi.asNA How to label cells that are "Pure" for multiple annotations: "Multi" (FALSE) or NA (TRUE)
#' @param additional.signatures UCell additional signatures to compute on each cell, by default \code{SignatuR} programs are included.
#' @param ref.maps A named list of the ProjecTILs reference maps to use. They ought to be Seurat objects. It is recommended to add in reference object slost \code{misc} the identifier connecting to layer 1 classification (scGate): \code{ref.map@misc$layer1_link}
#' @param split.by A Seurat object metadata column to split by (e.g. sample names).
#' @param layer1_link Column of metadata linking layer1 prediction (e.g. scGate ~ CellOntology_ID) in order to perform subsetting for second layer classification.
#' @param return.Seurat Whether to return a Hit object or add the data into Seurat object metadata
#' @param group.by If return.Seurat = F, variables to be used to summarize HiTME classification data in HiT object. See \link{get.HiTObject} for more information.
#' @param useNA Whether to include not annotated cells or not (labelled as "NA" in the group.by.composition) when summarizing into HiT object. Can be "no", "ifany", or "always". See \code{?table} and \link{get.celltype.composition} for more information. Default is "ifany".
#' @param remerge When setting split.by, if remerge = TRUE one object will be returned. If remerge = FALSE a list of objects will be returned.
#' @param species Define species used for to get correct gene signature list with SignatuR
#' @param ncores The number of cores to use, by default all available cores minus 2 are used.
#' @param bparam A \code{BiocParallel::bpparam()} object that tells Run.HiTME how to parallelize. If provided, it overrides the `ncores` parameter.
#' @param progressbar Whether to show a progressbar or not
#'
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom parallelly availableCores
#' @importFrom dplyr mutate filter %>%
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom scGate scGate get_scGateDB
#' @import SignatuR
#' @import scGate
#' @import ProjecTILs
#' @importFrom Seurat SplitObject
#' @importFrom data.table rbindlist setDT
#'
#' @return Seurat object with additional metadata showing cell type classification, or HiT object summarizing such classification.
#' @export Run.HiTME
#'
#' @examples
#' library(ProjecTILs)
#'
#' path_data <- file.path("~", "Dropbox", "CSI", "Standardized_SingleCell_Datasets", "ZhangY_2022_34653365", "output", "test")
#'
#' # Define scGate model
#' scGate_models_DB <- scGate::get_scGateDB(branch = "master", verbose = TRUE, force_update = TRUE)
#' models.TME <- scGate_models_DB$human$TME_HiRes
#'
#' # Load ProjecTILs reference maps
#' path_ref <- "~/Dropbox/CSI/reference_atlases"
#' ref.maps <- list(CD8 = load.reference.map(file.path(path_ref, "CD8T_human_ref_v1.rds")),
#'                  CD4 = load.reference.map(file.path(path_ref, "CD4T_human_ref_v2.rds")),
#'                  DC = load.reference.map(file.path(path_ref, "DC_human_ref_v1.rds")),
#'                  MoMac = load.reference.map(file.path(path_ref, "MoMac_human_v1.rds")))
#' Run.HiTME(object = seurat.object,
#'           scGate.model = models.TME,
#'           ref.maps = ref.maps,
#'           ncores = 6)
#'

Run.HiTME <- function(object = NULL,
                      scGate.model = "default",
                      scGate.model.branch = c("master", "dev"),
                      multi.asNA = TRUE,
                      additional.signatures = "default",
                      ref.maps = NULL,
                      split.by = NULL,
                      layer1_link = "CellOntology_ID",
                      return.Seurat = TRUE,
                      group.by = list("layer1" = c("scGate_multi"),
                                      "layer2" = c("functional.cluster")
                      ),
                      useNA = FALSE,
                      remerge = TRUE,
                      species = "human",
                      bparam = NULL,
                      ncores = parallelly::availableCores() - 2,
                      progressbar = TRUE) {

  if (is.null(object)) {
    stop("Please provide a Seurat object or a list of them")
  }

  if (is.list(object) &&
      !is.null(split.by)) {
    stop("split.by only supported for a single Seurat object, not a list.\n
         Merge list before running HiTME")
  }

  # split object into a list if indicated
  if (!is.null(split.by)) {
    if (is.list(object)) {
      stop("Split.by argument not supported when providing a list of Seurat objects. Set split.by = NULL or merge list.")
    }
    if (!split.by %in%
        names(object@meta.data)) {
      stop(paste("split.by argument: ", split.by, " is not a metadata column in this Seurat object"))
    }
    object <- Seurat::SplitObject(object, split.by = split.by)

  } else {
    # if not applying split.by not merge by default
    remerge <-  FALSE
  }

  if (!return.Seurat &&
      is.null(group.by)) {
    warning("If setting return Seurat as FALSE, HiT summarized object will be returned. Need to indicate group.by variable indicating cell type classification\n
         e.g. group.by = list(\"layer1\" = c(\"scGate_multi\"),\"layer2\" = c(\"functional.cluster\"))\n
            Returning Seurat object, not HiT object.")
    return.Seurat = TRUE
  }

  # adapt species
  if (is.null(species)) {
    stop("Please provide human or mouse as species")
  }
  species <- tolower(species)
  if (!species == "human") {
    if (grepl("homo|sap|huma", species)) {
      species <- "human"
    }
  } else if (!species == "mouse") {
    if (grepl("mice|mus", species)) {
      species <- "mouse"
    }
  } else {
    stop("Only supported species are human and mouse")
  }

  # if object is unique turn into a list
  if (!is.list(object)) {
    remerge <-  FALSE
    object <- list(object)
  }


  # Warning to reduce number of cores if file is huge
  if (any(lapply(object, ncol))>=30000) {
    warning("Huge Seurat object, consider reducing number of cores to avoid memory issues")
  }


  # set parallelization parameters
  param <- set_parallel_params(ncores = ncores,
                               bparam = bparam,
                               progressbar = progressbar)


  # Get default additional signatures
  # based on species get SignatuR signatures
  if (!is.null(additional.signatures)) {
    # convert to list if not
    if (!is.list(additional.signatures)) {
      additional.signatures <- list(additional.signatures)
    }

    if (length(additional.signatures) == 1 &&
        tolower(additional.signatures) == "default") {
      if (species == "human") {
        sig.species <- "Hs"
      } else if (species == "mouse") {
        sig.species <- "Mm"
      }
      additional.signatures <- SignatuR::GetSignature(SignatuR::SignatuR[[sig.species]][["Programs"]])
      selected.SignatuR.programs <- c("IFN", "HeatShock", "cellCycle.G1S", "cellCycle.G2M")
      additional.signatures <- additional.signatures[selected.SignatuR.programs]
      message(" - Adding additional signatures: ", paste(names(additional.signatures), collapse = ", "), "\n")
    }
  } else {
    for (v in seq_along(additional.signatures)) {
      if (is.null(names(additional.signatures)[[v]]) || is.na(names(additional.signatures)[[v]])) {
        names(additional.signatures)[[v]] <- paste0("Additional_signature", v)
      }
    }
  }


  # Run scGate, if model is provided
  if (!is.null(scGate.model)) {
    message("## Running scGate\n")

    if (!is.list(scGate.model)) {
      scGate.model <- list(scGate.model)
    }

    # Retrieve default scGate models if default
    if (length(scGate.model) == 1 &&
        tolower(scGate.model) == "default") {
      scGate.model.branch <- scGate.model.branch[1]
      if (species == "human") {
        scGate.model <- scGate::get_scGateDB(branch = scGate.model.branch,
                                             force_update = FALSE)[[species]][["HiTME"]]
      } else if (species == "mouse") {
        scGate.model <- scGate::get_scGateDB(branch = scGate.model.branch,
                                             force_update = FALSE)[[species]][["HiTME"]]
      }
      message(" - Running scGate model for ", paste(names(scGate.model), collapse = ", "), "\n")
    }



    ## Run scGate
    object <- lapply(
      X = object,
      function(x) {
        if (!inherits(x, "Seurat")) {
          stop("Not Seurat object included, cannot be processed.\n")
        }
        x <- scGate::scGate(x,
                            model=scGate.model,
                            additional.signatures = additional.signatures,
                            BPPARAM = param,
                            multi.asNA = multi.asNA)
        return(x)
      }
    )

    message("Finished scGate\n####################################################\n")
  } else {
    message("Not running coarse cell type classification as no scGate model was indicated.\n")
  }

  # Instance if we want to run additional signatures but not scGate
  if (is.null(scGate.model) &&
      !is.null(additional.signatures)) {
    message("Running additional Signatures but not Running scGate classification\n")
    object <- lapply(
      X = object,
      function(x) {
        if (!inherits(x, "Seurat")) {
          stop("Not Seurat object included, cannot be processed.\n")
        }
        x <- UCell::AddModuleScore_UCell(x,
                                         features = additional.signatures,
                                         BPPARAM = param)
        return(x)
      }
    )
  }



  # Run ProjecTILs if ref.maps is provided
  if (!is.null(ref.maps)) {
    # convert ref.maps to list if not
    if (!is.list(ref.maps)) {
      ref.maps <- list(ref.maps)
    }

    # give name to ref.maps list if not present
    for (v in seq_along(ref.maps)) {
      if (is.null(names(ref.maps)[[v]]) || is.na(names(ref.maps)[[v]])) {
        names(ref.maps)[[v]] <- paste0("Map_", v)
      }
    }

    # check that all ref maps are Seurat objects
    if (suppressWarnings(!all(lapply(ref.maps, function(x) {inherits(x, "Seurat")})))) {
      message("Some or all reference maps are not a Seurat object, please provide reference maps as Seurat objects.\nNot running Projectils.")
    } else {

      message("## Running Projectils\n")

      object <- lapply(
        X = object,
        function(x) {
          if (!inherits(x, "Seurat")) {
            stop("Not Seurat object included, cannot be processed.\n")
          }
          x <- ProjecTILs.classifier.multi(x,
                                           ref.maps = ref.maps,
                                           bparam = param,
                                           layer1_link = layer1_link)
          # Check if ProjecTIL columns were added
          # If not, add NA columns
          ProjecTILs_cols <- c("functional.cluster", "functional.cluster.conf")
          if (!any(ProjecTILs_cols %in% names(x@meta.data))) {
            x@meta.data[ProjecTILs_cols[!ProjecTILs_cols %in% names(x@meta.data)]] <- NA
          }
          return(x)
        }
      )
      message("Finished Projectils\n####################################################\n")
    }

  } else {
    message("Not running reference mapping as no reference maps were indicated.\n")
  }


  # Remerge if indicated
  if (remerge && length(object)>1) {
    object <- merge(object[[1]],object[2:length(object)])
    object <- list(object)
  }


  # add layer_links metadata levels
  # misc slots are removed when merging objects
  object <- lapply(object,
                   function(x) {
                     # add misc slot, removed when merging
                     x@misc[["layer1_param"]] <- list()
                     x@misc[["layer1_param"]][["scGate_models"]] <- names(scGate.model)
                     x@misc[["layer1_param"]][["additional.signatures"]] <- names(additional.signatures)

                     if (!is.null(scGate.model)) {
                       levs <- names(scGate.model)
                       if (!multi.asNA) {
                         levs <- c(levs, "Multi")
                       }
                       x$scGate_multi <- factor(x$scGate_multi,
                                                levels = levs)
                     }

                     if (!is.null(ref.maps)) {
                       # get ref.maps all cells types
                       if ("functional.cluster" %in% names(x@meta.data)) {
                         all.levels <- lapply(ref.maps, function(x) {
                           unique(x$functional.cluster)
                         })
                         x$functional.cluster <- factor(x$functional.cluster,
                                                        levels = unlist(all.levels))

                         # add each level to misc
                         names(all.levels) <- lapply(ref.maps, function(x) {
                           x@misc$layer1_link
                         })
                       } else {
                         all.levels <- NULL
                       }
                       x@misc[["layer2_param"]][["functional.cluster"]][["levels2_per_levels1"]] <- all.levels
                       # all refs indicated in the function
                       x@misc[["layer2_param"]][["functional.cluster"]][["References_user_specified"]] <- names(ref.maps)
                     }
                     return(x)
                   }
  )


  # if group.by parameters are not present in metadata, return Seurat
  if (!return.Seurat) {
    if (suppressWarnings(!all(lapply(object, function(x) {any(names(x@meta.data) %in% group.by)})))) {
      message("None of the 'group.by' variables found, returning Seurat object not HiT object.")
      return.Seurat <- TRUE
    } else {
      message("\nBuilding HiT object\n")

      hit <- BiocParallel::bplapply(
        X = object,
        BPPARAM = param,
        function(x) {
          x <- get.HiTObject(x,
                             group.by = group.by,
                             useNA = useNA
          )
        }
      )

    }
  }


  # if list is of 1, return object not list

  if (return.Seurat) {
    if (length(object)==1) {
      object <- object[[1]]
    }
    return(object)
  } else {
    if (length(hit)==1) {
      hit <- hit[[1]]
    }
    return(hit)
  }


}



# Function get.HiTObject
#' Get a HiT object summarizing cell type classification and aggregated profiles.
#'
#'
#'
#' @param object A Seurat object or a list of Seurat objects
#' @param group.by List or vector with one or multiple Seurat object metadata columns with cell type predictions to group by (e.g. layer 1 cell type classification)
#' @param split.by A Seurat object metadata column to split by compositional data (e.g. sample names) \link{get.celltype.composition}
#' @param min.cells.composition Parameter for internal \link{get.celltype.composition}. Minimum number of cells annotated in a group.by parameter to render the cell type composition. If the number of cells annotated for a certain group.by parameter is less, compositional data will not be rendered. Default value is 10.
#' @param min.cells.aggregated Parameter for internal \link{get.aggregated.profile} and \link{get.aggregated.signature}. Minimum number of cells per sample and cell type to aggregate data for pseudobulk or aggregated signatures. Aggregated data for cell types with less than that value will not be returned in aggregated data. Default value is 10.
#' @param name.additional.signatures Names of additional signatures as found in object metadata to take into account.
#' @param useNA logical whether to return aggregated profile for NA (undefined) cell types, default is FALSE.
#' @param clr_zero_impute_perc Parameter for internal \link{get.celltype.composition}.
#' @param layers_links Parameter for internal \link{get.celltype.composition}
#' @param gene.filter List of genes to subset for aggregated expression. Parameter for internal \link{get.aggregated.profile}.
#' @param assay Parameter for internal \link{get.aggregated.profile}.
#' @param layer1_link Column of metadata linking layer1 prediction (e.g. scGate ~ CellOntology_ID) in order to perform subsetting for second layer classification.
#'
#' @importFrom methods setClass new
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom parallelly availableCores
#'
#' @import SeuratObject
#'
#' @return HiT object summarizing cell type classification and aggregated profiles.
#' @export get.HiTObject
#'


get.HiTObject <- function(object,
                          group.by = list("layer1" = c("scGate_multi"),
                                          "layer2" = c("functional.cluster")
                          ),
                          split.by = NULL,
                          min.cells.composition = 10,
                          min.cells.aggregated = 10,
                          name.additional.signatures = NULL,
                          useNA = FALSE,
                          clr_zero_impute_perc = 1,
                          layers_links = c("scGate_multi" = "functional.cluster"),
                          gene.filter = NULL,
                          assay = "RNA",
                          layer1_link = "CellOntology_ID",

                          ncores = parallelly::availableCores() - 2,
                          bparam = NULL,
                          progressbar = TRUE) {

  args <- list(group.by,
               split.by,
               min.cells.composition,
               min.cells.aggregated,
               name.additional.signatures,
               useNA,
               clr_zero_impute_perc,
               layers_links,
               gene.filter,
               assay,
               layer1_link)

  # if object is a single Seurat object, turn into a list
  if (!is.list(object)) {
    hit.list <- get.HiTObject.helper(object, group.by,
                                     split.by,
                                     min.cells.composition,
                                     min.cells.aggregated,
                                     name.additional.signatures,
                                     useNA,
                                     clr_zero_impute_perc,
                                     layers_links,
                                     gene.filter,
                                     assay,
                                     layer1_link)
  } else {
    # set parallelization parameters
    param <- set_parallel_params(ncores = ncores,
                                 bparam = bparam,
                                 progressbar = progressbar)

    hit.list <- BiocParallel::bplapply(X = object,
                                       BPPARAM = param,
                                       function(x) {
                                         do.call(get.HiTObject.helper,
                                                 c(x, args)
                                         )
                                       })
  }

  return(hit.list)
}

get.HiTObject.helper <- function(object,
                                 group.by = list("layer1" = c("scGate_multi"),
                                                 "layer2" = c("functional.cluster")
                                 ),
                                 split.by = NULL,
                                 min.cells.composition = 10,
                                 min.cells.aggregated = 10,
                                 name.additional.signatures = NULL,
                                 useNA = FALSE,
                                 clr_zero_impute_perc = 1,
                                 layers_links = c("scGate_multi" = "functional.cluster"),
                                 gene.filter = NULL,
                                 assay = "RNA",
                                 layer1_link = "CellOntology_ID") {

  if (is.null(object)) {
    stop("Please provide a Seurat object")
  } else if (!inherits(object, "Seurat")) {
    stop("Not Seurat object included, cannot be processed.\n")
  }


  if (is.null(group.by)) {
    stop("Please provide at least one grouping variable")
  }

  if (!is.list(group.by)) {
    group.by <- as.list(group.by)
  }

  if (!any(group.by %in% names(object@meta.data))) {
    stop("Group.by variables ", paste(group.by, collapse = ", "), " not found in metadata")
  } else if (!all(group.by %in% names(object@meta.data))) {
    group.by <- group.by[group.by %in% names(object@meta.data)]
    message("Only found ", paste(group.by, collapse = ", ") , " as grouping variable for HiT Object.")
  }

  for (v in seq_along(group.by)) {
    if (is.null(names(group.by)[[v]]) || is.na(names(group.by)[[v]])) {
      names(group.by)[[v]] <- paste0("layer", v)
    }
  }


  # Make list of list for layers for predictions slot
  pred.list <- list()
  for (a in names(group.by)) {
    pred.list[[a]] <- list()
    pred.list[[a]][[group.by[[a]]]] <- object@meta.data[,group.by[[a]], drop = FALSE]
  }



  # Run extended if object got scGate info
  # Extract values from misc slot from object
  if (any(group.by == "scGate_multi")) {
    layer.scgate <- names(group.by[group.by == "scGate_multi"])
    if (is.null(name.additional.signatures)) {
      name.additional.signatures <- object@misc$layer1_param$additional.signatures
    }
    scgate.models <- object@misc$layer1_param$scGate_models


    if (!is.null(name.additional.signatures)) {
      sig <- grep(paste(name.additional.signatures, collapse = "|"),
                  names(object@meta.data), value = TRUE)
      sig.df <- object@meta.data[, sig]
      pred.list[[layer.scgate]][["additional_signatures"]] <- sig.df
    }

    if (!is.null(scgate.models)) {
      scgate <- grep(paste(paste0(scgate.models, "$"), collapse = "|"),
                     names(object@meta.data), value = TRUE)
      scgate.df <- object@meta.data[, scgate]
      pred.list[[layer.scgate]][["scGate_is.pure"]] <- scgate.df
    }

    if (!is.null(name.additional.signatures) && !is.null(scgate.models)) {
      ucell <- names(object@meta.data)[!names(object@meta.data) %in% c(sig, scgate)] %>%
        grep("_UCell$", ., value = TRUE)

      ucell.df <- object@meta.data[, ucell]

      pred.list[[layer.scgate]][["UCell_scores"]] <- ucell.df
    }
  }

  # Run extended if object got ProjecTILs info
  # Extract values from misc slot from object
  if (any(group.by == "functional.cluster")) {
    layer.pt <- names(group.by[group.by == "functional.cluster"])
    pred.list[[layer.pt]][["functional.cluster"]] <- object@meta.data[,"functional.cluster", drop = FALSE]
    pred.list[[layer.pt]][["functional.cluster.conf"]] <- object@meta.data[,"functional.cluster.conf", drop = FALSE]
  }


  # Compute proportions
  # message("\nComputing cell type composition...\n")

  comp.prop <- get.celltype.composition(object,
                                        group.by.composition = group.by,
                                        min.cells.composition = min.cells.composition,
                                        split.by = split.by,
                                        useNA = useNA,
                                        clr_zero_impute_perc = clr_zero_impute_perc,
                                        layer1_link = layer1_link,
                                        layers_links = layers_links)
  # Compute avg expression
  # message("\nComputing aggregated profile...\n")

  avg.expr <- get.aggregated.profile(object,
                                     group.by.aggregated = group.by,
                                     gene.filter = gene.filter,
                                     min.cells.aggregated = min.cells.aggregated,
                                     assay = assay,
                                     useNA = useNA)

  aggr.signature <- get.aggregated.signature(object,
                                             group.by.aggregated = group.by,
                                             min.cells.aggregated = min.cells.aggregated,
                                             name.additional.signatures = name.additional.signatures,
                                             useNA = useNA)

  hit <- methods::new("HiT",
                      metadata = object@meta.data,
                      predictions = pred.list,
                      aggregated_profile = list("pseudobulk" = avg.expr,
                                                "signatures" = aggr.signature),
                      composition = comp.prop
  )
  return(hit)
}


#' Calculate cell type composition or frequencies
#'
#' @param object A seurat object
#' @param group.by.composition The Seurat object metadata column(s) containing celltype annotations (provide as character vector, containing the metadata column name(s))
#' @param split.by A Seurat object metadata column to split by (e.g. sample names)
#' @param min.cells.composition Set a minimum threshold for number of cells to calculate relative abundance (e.g. less than 10 cells -> no relative abundance will be calculated)
#' @param useNA Whether to include not annotated cells or not (labelled as "NA" in the group.by.composition). Can be defined separately for each group.by.composition (provide single boolean or vector of booleans)
#' @param layers_links Named vector indicating the relation of multiple layer classification, default is \code{c("scGate_multi", "functional.cluster")}. Name of the vector element ought to be layer1 and element layer2. This parameter is used to compute relative compositional data of layer2 within layer1 classes.
#' @param clr_zero_impute_perc To calculate the clr-transformed relative abundance ("clr_freq"), zero values are not allowed and need to be imputed (e.g. by adding a pseudo cell count). Instead of adding a pseudo cell count of flat +1, here a pseudo cell count of +1% of the total cell count will be added to all cell types, to better take into consideration the relative abundance ratios (e.g. adding +1 cell to a total cell count of 10 cells would have a different, i.e. much larger effect, than adding +1 to 1000 cells).
#' @param layer1_link Column of metadata linking layer1 prediction (e.g. scGate ~ CellOntology_ID) in order to perform subsetting for second layer classification.

#' @importFrom Hotelling clr
#' @importFrom dplyr group_by summarize filter ungroup mutate select left_join n coalesce bind_rows across all_of
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom rrapply rrapply
#'
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
#'
get.celltype.composition <- function(object = NULL,
                                     group.by.composition = NULL,
                                     split.by = NULL,
                                     min.cells.composition = 10,
                                     useNA = FALSE,
                                     layers_links = c("scGate_multi" = "functional.cluster"),
                                     clr_zero_impute_perc = 1,
                                     layer1_link = "CellOntology_ID") {

  if (is.null(object)) {
    stop("Please provide a Seurat object or metadata as dataframe")
  }

  # input can be a Seurat object or a dataframe containing its meta.data
  # convert object to metadata if seurat object is provided
  if (inherits(object, "Seurat")) {
    meta.data <- object@meta.data
    if (is.null(meta.data)) {
      stop("No metadata found in this Seurat object")
    }
  } else {
    stop("Not Seurat object or dataframe included, cannot be processed.\n")
  }

  if (is.null(group.by.composition)) {
    stop("Please specificy a group.by.composition variable")
  }

  # Assess wheter split.by variable is in metadata
  if (!is.null(split.by) &&
      !split.by %in% names(meta.data)) {
    stop("Split.by variable not found in meta.data!\n")
  }

  if (length(useNA) != 1 &&
      length(useNA) != length(group.by.composition)) {
    stop("useNA variable must be of length 1 or the same length as group.by.composition (group.by)")
  }

  # convert group.by.composition to list if not
  if (!is.list(group.by.composition)) {
    group.by.composition <- as.list(group.by.composition)
  }


  # Rename group.by.composition if not indicated
  for (v in seq_along(group.by.composition)) {
    if (is.null(names(group.by.composition)[[v]]) ||
        is.na(names(group.by.composition)[[v]])) {
      names(group.by.composition)[[v]] <- group.by.composition[[v]]
    }

  }



  # Keep only grouping variables in metadata
  if (!any(group.by.composition %in% names(meta.data))) {
    stop("Group.by variables ", paste(group.by.composition, collapse = ", "), " not found in metadata")
  } else if (!all(group.by.composition %in% names(meta.data))) {
    group.present <- group.by.composition %in% names(meta.data)
    # accommodate useNA vector
    if (length(useNA) == length(group.by.composition)) {
      useNA <- useNA[group.present]

    }
    # keep only grouping variables found in metadata
    group.by.composition <- group.by.composition[group.present]
    message("Only found ", paste(group.by.composition, collapse = ", ") , " as grouping variables.")
  }

  # Factorize cell type to keep all in compositional data
  meta.data <- meta.data %>%
    dplyr::mutate_at(unlist(group.by.composition), as.factor)

  # evaluate which layers are included in group.by.composition and layers_links
  for (l in seq_along(layers_links)) {
    if (!all(c(names(layers_links[l]), layers_links[l]) %in%
             group.by.composition)) {
      message("Not computing relative proportion values of layer2 within layer1 for ",
              paste0(names(layers_links[l]), " -> ", layers_links[l]))
    }
  }


  # Rep useNa parameters according to group.by.composition variable
  if (length(useNA) == 1) {
    useNA <- rep(useNA, length(group.by.composition))
  }

  # list to store compositions
  celltype.compositions <- list()


  for (i in seq_along(group.by.composition)) {

    suppressMessages(
      {
        ctable <- compositional_data(data = meta.data,
                                     split.by = split.by,
                                     group.by.1 = group.by.composition[[i]],
                                     useNA = useNA[i],
                                     clr_zero_impute_perc = clr_zero_impute_perc)

        # stop if very few cells
        if (sum(ctable$cell_counts) < min.cells.composition) {
          # Return empty data.frame
          ctable <- ctable[0,]
        }
        # get proportion relative to layer1 types
        if (group.by.composition[[i]] %in% layers_links) {
          # stop if very few cells
          if (sum(ctable$cell_counts) < min.cells.composition) {
            # Return empty list
            ctable <- list()
          } else {
            lay <- layers_links[which(layers_links == group.by.composition[[i]])]
            if (lay %in% names(object@misc$layer2_param)) {
              meta.split <- split(meta.data,
                                  meta.data[[names(lay)]])
              # switch list names to cell ontology ID
              cellonto_dic <- lapply(meta.split, function(x) {
                nam <- unique(x[[names(lay)]])
                val <- unique(x[[layer1_link]])
                names(val) <- nam
                return(val)
              }) %>%
                unname() %>%
                unlist()

              levs <- object@misc$layer2_param[[lay]]$levels2_per_levels1
              names(levs) <- names(cellonto_dic[match(names(levs), cellonto_dic)])

              # If a celltype was not detected, drop it
              levs <- levs[!is.na(names(levs))]

              # Filter only layer1 cells with representation in layer2
              meta.split <- meta.split[names(levs)]

              # add factor to group by layer1
              for (f in names(meta.split)) {
                meta.split[[f]]$functional.cluster <- factor(meta.split[[f]]$functional.cluster,
                                                             levels = levs[[f]])
              }

              ctable.split <- lapply(meta.split,
                                     compositional_data,
                                     split.by = split.by,
                                     group.by.1 = group.by.composition[[i]],
                                     useNA = useNA[i])

              # If not enough cells, return empty dataframe
              ctable.split <- lapply(ctable.split,
                                     function (x) {
                                       if (sum(x$cell_counts) < min.cells.composition) {
                                         x[0,]
                                       } else {x}
                                     })

              ctable <- c(ctable.split)
            }
          }
        }
      })

    ## Append
    celltype.compositions[[names(group.by.composition)[i]]] <- ctable
  }

  return(celltype.compositions)

}



#' Compute aggregated gene expression
#'
#' Function to compute aggregated expression (pseudobulk, i.e. sum counts per ident), and average expression by indicated cell type or grouping variable.
#'
#'
#' @param object A seurat object or a list of seurat objects
#' @param group.by.aggregated The Seurat object metadata column(s) containing celltype annotations
#' @param min.cells.aggregated Minimum number of cells per sample and cell type to aggregate data for pseudobulk or aggregated signatures. Aggregated data for cell types with less than that value will not be returned in aggregated data. Default value is 10.
#' @param assay Assay to retrieve information. By default "RNA".
#' @param useNA logical whether to return aggregated profile for NA (undefined) cell types, default is FALSE.
#' @param ... Extra parameters for internal Seurat functions: AverageExpression, AggregateExpression, FindVariableFeatures

#' @importFrom Seurat AverageExpression AggregateExpression FindVariableFeatures
#' @importFrom dplyr filter
#' @return Average and aggregated expression as a list of matrices for all genes and indicated gene lists filtering.
#' @export get.aggregated.profile


get.aggregated.profile <- function(object,
                                   group.by.aggregated = NULL,
                                   min.cells.aggregated = 10,
                                   assay = "RNA",
                                   useNA = FALSE,
                                   ...) {

  if (is.null(object)) {
    stop("Please provide a Seurat object")
    if (!inherits(object, "Seurat")) {
      stop("Please provide a Seurat object")
    }
  }

  if (is.null(group.by.aggregated)) {
    stop("Please specificy a group.by.aggregated variable")
  }

  # convert group.by.aggregated to list if not
  if (!is.list(group.by.aggregated)) {
    group.by.aggregated <- as.list(group.by.aggregated)
  }

  # Rename group.by.aggregated if not indicated
  for (v in seq_along(group.by.aggregated)) {
    if (is.null(names(group.by.aggregated)[[v]]) || is.na(names(group.by.aggregated)[[v]])) {
      names(group.by.aggregated)[[v]] <- group.by.aggregated[[v]]
    }
  }

  if (!any(group.by.aggregated %in% names(object@meta.data))) {
    stop("Group.by variables ", paste(group.by.aggregated, collapse = ", "), " not found in metadata")
  } else if (!all(group.by.aggregated %in% names(object@meta.data))) {
    group.by.aggregated <- group.by.aggregated[group.by.aggregated %in% names(object@meta.data)]
    message("Only found ", paste(group.by.aggregated, collapse = ", ") , " as grouping variables.")
  }


  avg.exp <- list()

  # loop over different grouping
  for (i in names(group.by.aggregated)) {

    # compute pseudobulk
    suppressWarnings({

      # Calculate total (sample) pseudobulks
      if (i == names(group.by.aggregated)[1]) {

        # Calculate pseudobulk for ALL cells in sample
        avg.exp[[i]] <- object@assays[["RNA"]]["counts"]
        row_names <- row.names(avg.exp[[i]])
        avg.exp[[i]] <- Matrix::Matrix(rowSums(avg.exp[[i]]))
        row.names(avg.exp[[i]]) <- row_names
        colnames(avg.exp[[i]]) <- "all"

        # Calculate pseudobulk for only annotated cells in sample

        # Subset to remove not annotated cells
        object <- object[, which(!is.na(object[[group.by.aggregated[[i]]]]))]

        # Calculate pseudobulk of all annotated cells
        mat <- object@assays[["RNA"]]["counts"]
        row_names <- row.names(mat)
        mat <- Matrix::Matrix(rowSums(mat))
        row.names(mat) <- row_names
        colnames(mat) <- "all.annotated_only"

        avg.exp[[i]] <- cbind(avg.exp[[i]], mat)
      }

      # remove from aggregated data cell with less than min.cells.aggregated
      cnts <- compositional_data(object@meta.data,
                                 group.by.1 = group.by.aggregated[[i]],
                                 only.counts = TRUE,
                                 useNA = useNA)

      if (nrow(cnts) == 0) {
        next
      }

      # Remove cell types with less than min.cells.aggregated
      keep <- cnts[cnts[["cell_counts"]] > min.cells.aggregated, 1] %>%
        unlist()
      object@meta.data$fltr <- object@meta.data[[group.by.aggregated[[i]]]]
      object <- object[, as.character(object$fltr) %in% keep]

      # Handle not annotated cells being labelled as NA
      # object@meta.data[[group.by.aggregated[[i]]]] <- as.character(object@meta.data[[group.by.aggregated[[i]]]])
      # object@meta.data[[group.by.aggregated[[i]]]][is.na(object@meta.data[[group.by.aggregated[[i]]]])] <- "NA"
      # object@meta.data[[group.by.aggregated[[i]]]] <- as.factor(object@meta.data[[group.by.aggregated[[i]]]])

      if (length(unique(object@meta.data[[group.by.aggregated[[i]]]])) >= 2) {
        mat <-
          Seurat::AggregateExpression(object,
                                      group.by = group.by.aggregated[[i]],
                                      assays = assay,
                                      verbose = FALSE,
                                      ...)[[assay]]
        avg.exp[[i]] <- cbind(avg.exp[[i]], mat)
      } else {
        # Handle case if there is only one cell type
        col_name <- as.character(unique(object@meta.data[[group.by.aggregated[[i]]]]))
        colnames(avg.exp[[i]]) <- col_name
      }
    })

    if (useNA == FALSE) {
      # Drop NA column
      avg.exp[[i]] <- avg.exp[[i]][, !colnames(avg.exp[[i]]) %in% "NA", drop = FALSE]
    }
  }
  return(avg.exp)
}


#' Compute aggregated additional signatures by cell type
#'
#' Function to compute aggregated signatures of predicted cell types.
#'
#'
#' @param object A seurat object or metadata data frame.
#' @param group.by.aggregated The Seurat object metadata column(s) containing celltype annotations (idents).
#' @param min.cells.aggregated Minimum number of cells per sample and cell type to aggregate data for pseudobulk or aggregated signatures. Aggregated data for cell types with less than that value will not be returned in aggregated data. Default value is 10.
#' @param name.additional.signatures Names of additional signatures to compute the aggregation per cell type.
#' @param fun Function to aggregate the signature, e.g. mean or sum.
#' @param useNA logical whether to return aggregated signatures for NA (undefined) cell types, default is FALSE.

#' @importFrom dplyr group_by summarize_at filter
#' @return Aggregated signature score for each indicated cell type grouping Results is NULL of not additional signatures are indicated or present in metadata.
#' @export get.aggregated.signature


get.aggregated.signature <- function(object,
                                     group.by.aggregated = NULL,
                                     min.cells.aggregated = 10,
                                     name.additional.signatures = NULL,
                                     fun = mean,
                                     useNA = FALSE) {

  # input can be a Seurat object or a dataframe containing its meta.data
  # convert object to metadata if seurat object is provided
  if (inherits(object, "Seurat")) {
    meta.data <- object@meta.data
    if (is.null(meta.data)) {
      stop("No metadata found in this Seurat object")
    }
  } else if (inherits(object, "data.frame")) {
    meta.data <- object
  } else {
    stop("Not Seurat object or dataframe included, cannot be processed.\n")
  }


  if (is.null(group.by.aggregated)) {
    stop("Please specificy a group.by.aggregated variable")
  }

  # convert group.by.aggregated to list if not
  if (!is.list(group.by.aggregated)) {
    group.by.aggregated <- as.list(group.by.aggregated)
  }

  # Rename group.by.aggregated if not indicated
  for (v in seq_along(group.by.aggregated)) {
    if (is.null(names(group.by.aggregated)[[v]]) ||
        is.na(names(group.by.aggregated)[[v]])) {
      names(group.by.aggregated)[[v]] <- group.by.aggregated[[v]]
    }
  }

  if (!any(group.by.aggregated %in% names(meta.data))) {
    stop("Group.by variables ", paste(group.by.aggregated, collapse = ", "), " not found in metadata")
  } else if (!all(group.by.aggregated %in% names(meta.data))) {
    group.by.aggregated <- group.by.aggregated[group.by.aggregated %in% names(meta.data)]
    message("Only found ", paste(group.by.aggregated, collapse = ", ") , " as grouping variable.")
  }

  if (is.null(name.additional.signatures)) {
    name.additional.signatures <- object@misc$layer1_param$additional.signatures
  }

  if (is.null(name.additional.signatures)) {
    message("No additional signatures indicated. Returning NULL")
    aggr.sig <- NULL
  } else {

    if (!any(grepl(paste(name.additional.signatures, collapse = "|"),
                   names(meta.data)))) {
      stop("No additional signatures found in this object metadata")
    }

    add.sig.cols <- grep(paste(name.additional.signatures, collapse = "|"),
                         names(meta.data), value = TRUE)

    if (length(add.sig.cols) > length(name.additional.signatures)) {
      for (i in name.additional.signatures) {
        if (sum(grep(i, names(meta.data))) > 1) {
          meta_cols <- names(meta.data)[grep(i, names(meta.data))]
          message(paste("Signatue", i, "was found in multiple metadata columns:", meta_cols))
        }
      }
      stop("The name of at least one signature provided was found in multiple metadata columns.
           Please give them a more unique name, e.g. by appending '_signature' to the name")
    }

    aggr.sig <- list()

    for (e in names(group.by.aggregated)) {
      # remove from aggregated data cell with less than min.cells.aggregated
      cnts <- compositional_data(meta.data,
                                 group.by.1 = group.by.aggregated[[e]],
                                 only.counts = TRUE,
                                 useNA = useNA)

      keep <- cnts[cnts[["cell_counts"]] > min.cells.aggregated, 1] %>%
        unlist()


      aggr.sig[[e]] <- meta.data %>%
        dplyr::filter(.data[[group.by.aggregated[[e]]]] %in% keep) %>%
        dplyr::group_by(.data[[group.by.aggregated[[e]]]]) %>%
        dplyr::summarize_at(add.sig.cols, fun, na.rm = TRUE)

      # filter out NA if useNA=F
      if (!useNA) {
        aggr.sig[[e]] <- aggr.sig[[e]] %>%
          dplyr::filter(!is.na(.data[[group.by.aggregated[[e]]]]))
      }

      colnames(aggr.sig[[e]])[1] <- "celltype"
    }
  }
  return(aggr.sig)
}


#' Retrieve list of genes from GO accessions
#'
#' Function to fetch GO accessions gene list from biomaRt.
#'
#'
#' @param GO_accession Vector of GO accession ID (e.g. \code{c("GO:0003700","GO:0005125")}). Vector could be named for clarity in the output.
#' @param species species to retrieve the genes, for now suppported human and mice.
#' @param host BioMart host to connect. Default is \link{https://dec2021.archive.ensembl.org/}, as it fails less.

#' @importFrom dplyr filter
#' @importFrom biomaRt useMart getBM

#' @return List of genes for each GO accession requested. If named vector is provided, lists output are named as GO:accession.ID_named (e.g. "GO:0004950_cytokine_receptor_activity")
#' @export get.GOList


get.GOList <- function(GO_accession = NULL,
                       species = "Human",
                       host = "https://dec2021.archive.ensembl.org/") {

  if (is.null(GO_accession)) {
    stop("Please provide at least one GO accession ID (e.g. GO:0003700")
  } else {
    if (length(names(GO_accession)) == 0) {
      names(GO_accession) <- GO_accession
    }
  }

  species <- tolower(species)


  # adapt species
  if (is.null(species)) {
    stop("Please provide human or mouse as species")
  }
  species <- tolower(species)
  if (grepl("homo|sapi|huma", species)) {
    dataset <- "hsapiens_gene_ensembl"
  } else if (grepl("mice|mus", species)) {
    dataset <- "mmusculus_gene_ensembl"
  } else {
    stop("Only supported species are human and mouse")
  }
  # Load default GO accession in case cannot connect to biomaRt
  # gene.data.bm <- data("GO_accession_default.RData")

  # Retrieve genes for each GO
  gene.data.bm <-
    tryCatch(
      {

        ensembl = biomaRt::useMart("ensembl",
                                   dataset = dataset,
                                   host = host)

        gene.data.bm <- biomaRt::getBM(attributes=c('hgnc_symbol',
                                                    'go_id'),
                                       filters = 'go',
                                       values = GO_accession,
                                       mart = ensembl)
        gene.data.bm
      },
      error = function(e) {
        message("GO accession from biomaRt not possible")
        message(e)
        NULL
      }
    )

  if (!is.null(gene.data.bm)) {
    gene.data <- gene.data.bm %>%
      dplyr::filter(go_id %in% GO_accession) %>%
      dplyr::filter(hgnc_symbol != "") %>%
      split(., as.factor(.$go_id)) %>%
      lapply(., function(x) x[,1])

    if (length(names(gene.data)) != length(GO_accession)) {
      warning("Additional GO accession provided not found in GO database")
    }
    names(gene.data) <- paste0(names(gene.data), "_",
                               names(GO_accession)[GO_accession %in% names(gene.data)]) %>%
      gsub(" ", "_", .)
  } else {
    gene.data <- NULL
  }
  return(gene.data)
}


#' Merge HiTObjects
#'
#' @param hit.object List of HiTObjects
#' @param group.by If only merging for certain layers of annotation is intended, layers names can be indicated here as vector. Otherwise all layers present in all HiT object will be merged.
#' @param metadata.vars Variables to keep as metadata. (Default: NULL, keeping unique metadata columns per sample, dropping single-cell metadata)
#' @param pseudobulk.matrix Paramater to determine whther obtain the pseudobulk matrix as a single matrix (\code{"unique"}), or as one matrix for each cell type in the layer (\code{"list"})
#' @param ncores The number of cores to use, by default, all available cores - 2.
#' @param bparam A \code{BiocParallel::bpparam()} object that tells how to parallelize. If provided, it overrides the `ncores` parameter.
#' @param progressbar Whether to show a progressbar or not
#' @param verbose Whether to show optional messages or not

#' @importFrom dplyr mutate mutate_if filter %>% mutate_all full_join
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom parallelly availableCores
#' @importFrom data.table rbindlist
#' @importFrom methods slot

#' @return Merged HiTObject
#' @export merge.HiTObjects
#'

merge.HiTObjects <- function(hit.object = NULL,
                             group.by = NULL,
                             metadata.vars = NULL,
                             pseudobulk.matrix = "list",
                             ncores = parallelly::availableCores() - 2,
                             bparam = NULL,
                             progressbar = FALSE,
                             verbose = FALSE) {

  if (is.null(hit.object) ||
      !is.list(hit.object)) {
    stop("Please provide a list of HiT object containing more than one sample")
    if (suppressWarnings(!all(lapply(hit.object, function(x) {inherits(x, "HiT")})))) {
      stop("Not all components of the list are HiT objects.")
    }
  }

  # TODO: delete? -------------------------------
  # # give name to list of hit objects
  for (v in seq_along(hit.object)) {
    if (is.null(names(hit.object)[v]) ||
        is.na(names(hit.object)[v])) {
      names(hit.object)[v] <- paste0("Sample", v)
    }
  }

  # fetch which layers are included in all HiT objects
  layers_in_hit <- lapply(hit.object, function(x) {
    a <- names(x@composition)
    b <- names(x@aggregated_profile$pseudobulk)
    c <- names(x@aggregated_profile$signatures)
    u <- unique(c(a,b,c))
    return(u)
  })

  # if group.by is not NULL retrieve the layers present in HiT object
  if (is.null(group.by)) {
    group.by <- unique(unlist(layers_in_hit))
  }


  if (suppressWarnings(!all(lapply(layers_in_hit, function(x) {any(group.by %in% x)})))) {
    stop("None of the supplied HiT object contain at least one of these layers: ", paste(group.by, collapse = ", "),
         ".\nPlease make sure to indicate the name of the layer.")
  } else {
    present <- sapply(group.by,
                      function(char) {
                        unlist(lapply(layers_in_hit,
                                      function(df) {char %in% df}
                        ))
                      }
    )

    present.sum <- colSums(present)

    for (l in names(group.by)) {
      message("** ", l, " present in ", present.sum[[l]], " / ", length(hit.object), " HiT objects.")
    }
  }

  if (!is.null(metadata.vars)) {
    if (verbose) {message("\n#### Metadata ####")}
    if (suppressWarnings(!all(lapply(hit.object, function(x) {any(metadata.vars %in% names(x@metadata))})))) {
      message("Not all supplied HiT object contain ", paste(metadata.vars, collapse = ", "),
              "metadata elements in their metadata")
    }
    if (verbose) {
      in.md <- sapply(metadata.vars,
                      function(char) {
                        unlist(lapply(hit.object,
                                      function(df) {char %in% colnames(df@metadata)}
                        ))
                      }
      )
      in.md.sum <- colSums(in.md)

      for (l in metadata.vars) {
        message("** ", l, " present in ",
                in.md.sum[[l]], " / ", length(hit.object),
                " HiT objects.")
      }
    }
  }

  # set parallelization parameters
  param <- set_parallel_params(ncores = ncores,
                               bparam = bparam,
                               progressbar = progressbar)

  # Join data from all HiT objects

  # Metadata
  if (is.null(metadata.vars)) {
    # Per HiTObject from input, get metadata column names which have all the same value.
    # E.g. each HiTObject is from a specific sample, so each HiTObject has metadata column "Sample" which contains all the same value (e.g. "sample1" for sample 1, "sample2" for sample 2, etc.)
    # Other metadata columns, e.g. scGate_multi can be dropped, as the merged HiTObject is a summary per sample, so single-cell metadata columns don't make sense.
    all.names <- list()
    for (x in names(hit.object)) {
      all.names[[x]] <- names(hit.object[[x]]@metadata)
    }
    all.names_in.all <- Reduce(intersect, all.names)

    umc <- list()
    for (x in names(hit.object)) {
      umc[[x]] <- apply(hit.object[[x]]@metadata[all.names_in.all], 2, function(y) length(unique(y))) == 1
    }
    umc <- umc %>% Reduce("&", .)
    metadata.vars <- names(umc)[umc == TRUE]
  }
  metadata <- lapply(names(hit.object),
                     function (x) {hit.object[[x]]@metadata[1, metadata.vars] %>%
                         mutate(hitme.sample = x)
                     }
  )
  metadata <- data.table::rbindlist(metadata, use.names=TRUE, fill=TRUE)

  comp.prop <- list()
  avg.expr <- list()
  aggr.signature <- list()

  for (gb in group.by) {

    layer_present <- row.names(present)[present[,gb]]

    # Composition
    message("Merging compositions of " , gb, "...")

    type <- "composition"

    # assess that all HiT objects in the list contain the composition data for that layer
    is_df_check <- hit.object[layer_present] %>%
      lapply(methods::slot, name = type) %>%
      lapply("[[", gb) %>%
      lapply(is.data.frame) %>%
      unlist()
    is_list_check <- hit.object[layer_present] %>%
      lapply(methods::slot, name = type) %>%
      lapply("[[", gb) %>%
      lapply(function(x) {is.list(x) &!inherits(x, "data.frame")}) %>%
      unlist()

    if (all(is_df_check)) {
      df <- BiocParallel::bplapply(X = layer_present,
                                   BPPARAM = param,
                                   function(x) {
                                     hit.object[[x]]@composition[[gb]] %>%
                                       mutate(hitme.sample = x)
                                   })
      df <- data.table::rbindlist(df, use.names=TRUE, fill=TRUE)
      comp.prop[[gb]] <- df

    } else if (all(is_list_check)) {
      gb_sublevel_unique_names <- hit.object[layer_present] %>%
        lapply(methods::slot, name = type) %>%
        lapply("[[", gb) %>%
        lapply(names) %>%
        unlist() %>%
        unique()
      for (i in gb_sublevel_unique_names) {
        df <- BiocParallel::bplapply(X = layer_present,
                                     BPPARAM = param,
                                     function(x) {
                                       if (!is.null(hit.object[[x]]@composition[[gb]][[i]])) {
                                         hit.object[[x]]@composition[[gb]][[i]] %>%
                                           mutate(hitme.sample = x)
                                       }
                                     })
        df <- data.table::rbindlist(df, use.names=TRUE, fill=TRUE)
        comp.prop[[gb]][[i]] <- df
      }
    }


    # Aggregated_profile
    message("Merging aggregated profiles of " , gb, "...")

    type <- "pseudobulk"

    # Aggregated_profile Pseudobulk
    celltypes <- lapply(names(hit.object),
                        function(x) {
                          colnames(hit.object[[x]]@aggregated_profile[[type]][[gb]])
                        }) %>%
      unlist() %>% unique()

    for (ct in celltypes) {
      ct_present <- lapply(names(hit.object),
                           function(x) {
                             ct %in% colnames(hit.object[[x]]@aggregated_profile[[type]][[gb]]) }) %>%
        unlist()
      layer_ct_present <- names(hit.object)[(names(hit.object) %in% layer_present) & ct_present]
      df <- BiocParallel::bplapply(X = layer_ct_present,
                                   BPPARAM = param,
                                   function(x) {
                                     hit.object[[x]]@aggregated_profile[[type]][[gb]][, ct]
                                   })
      df <- do.call(cbind, df)
      df <- Matrix::Matrix(df, sparse = TRUE)
      colnames(df) <- layer_ct_present
      avg.expr[[gb]][[ct]] <- df
    }

    # join each matrix per celltype into a single matrix changing the colnames to accommodate sample source
    if (tolower(pseudobulk.matrix) == "unique") {
      avg.expr[[gb]] <- lapply(names(avg.expr[[gb]]),
                               function(x) {
                                 mat <- avg.expr[[gb]][[x]]
                                 colnames(mat) <- paste(x, colnames(mat), sep = "__")
                                 mat <- mat %>% as.data.frame() %>%
                                   tibble::rownames_to_column("gene")
                               }) %>%
        reduce(full_join, by = "gene") %>%
        # convert NA to 0
        mutate_if (is.numeric, ~ifelse(is.na(.), 0, .)) %>%
        tibble::column_to_rownames("gene") %>%
        as.matrix() %>%
        Matrix::Matrix(., sparse = TRUE)
    }


    # Aggregated_profile Signatures

    type <- "signatures"

    df <- BiocParallel::bplapply(X = layer_present,
                                 BPPARAM = param,
                                 function(x) {
                                   if (!is.null(hit.object[[x]]@aggregated_profile[[type]][[gb]])) {
                                     hit.object[[x]]@aggregated_profile[[type]][[gb]] %>%
                                       mutate(hitme.sample = x)
                                   }
                                 })
    df <- data.table::rbindlist(df,
                                use.names = TRUE,
                                fill = TRUE)
    aggr.signature[[gb]] <- df
  }

  hit <- methods::new("HiT",
                      metadata = metadata,
                      composition = comp.prop,
                      aggregated_profile = list("pseudobulk" = avg.expr,
                                                "signatures" = aggr.signature)
  )
  return(hit)
}


#' Get metrics of cell types pseudobulk clustering
#'
#' @param hit.object A Hit class object obtained with \link{get.HiTObject} or pseudobulk raw count matrix (\code{HitObject@aggregated_profile$pseudobulk$layer1})
#' @param cluster.by Vector indicating the variable for clustering, default is celltype (for the annotation) and sample
#' @param cluster.by.drop.na Whether to keep (FALSE) or drop (TRUE) NAs present in cluster.by column.
#' @param batching Vector indicating the variable for batching to allow calculating scores per batch, to account for batch effect
#' @param scores Scores to compute clustering of samples based on cell type prediction.
#' @param modularity.k Number of k-nearest neighbours to use to build graph for the calculation of the modularity score
#' @param dist.method Method to compute distance between celltypes, default is euclidean.
#' @param hclust.method Hierarchical clustering method for hclust. Options are: "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC). (See hclust package for details)
#' @param ntests Number of shuffling events to calculate p-value for scores
#' @param seed Set seed for random shuffling events to calculate p-value for scores
#' @param pval.combine Method for combining p-values if calculated using batching. Default is "zmethod" weighted (Stouffer's) Z-method, weighting by sample size per batch. Alternatively, Fisher's method can be used with "fisher".
#' @param pca_comps_labs_invisible Parameter to pass to fviz_pca defining which labels to plot for composition PCA, see documentation of fviz_pca for details.
#' @param pca_pb_labs_invisible Parameter to pass to fviz_pca defining which labels to plot for pseudobulk PCA, see documentation of fviz_pca for details.
#' @param pca_sig_labs_invisible Parameter to pass to fviz_pca defining which labels to plot for signatures PCA, see documentation of fviz_pca for details.
#' @param ndim Number of dimensions to be use for PCA clustering metrics. Default is 10.
#' @param nVarGenes Number of variable genes to assess samples. Default is 500.
#' @param gene.filter Named list of genes to subset their aggregated expression. Default is \code{"HVG"} or \code{NULL} would indicate "Highly variable genes", number of genes can be set at \code{nVarGenes}. If "default_filter" is indicated transcription factor (GO:0003700), cytokines (GO:0005125), cytokine receptors (GO:0004896), chemokines (GO:0008009), and chemokines receptors (GO:0004950) are subsetted out of the total genes and accounted for. For additional GO accession gene list you may use the function \link{getGOList} to fetch them from biomaRt.
#' @param black.list List of genes to discard from clustering, if "default" object "default_black_list" object is used. Alternative black listed genes can be provided as a vector or list.
#' @param ncores The number of cores to use, by default, all available cores - 2.
#' @param bparam A \code{BiocParallel::bpparam()} object that tells how to parallelize. If provided, it overrides the `ncores` parameter.
#' @param progressbar Whether to show a progressbar or not

#' @importFrom tidyr separate pivot_wider
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom parallelly availableCores
#' @importFrom data.table rbindlist
#' @importFrom dplyr mutate mutate_if filter %>% coalesce mutate_all full_join row_number
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom tibble rownames_to_column column_to_rownames remove_rownames
#' @importFrom caret nearZeroVar
#' @importFrom ggplot2 aes geom_point guides theme geom_col labs guide_legend annotate theme_bw ggtitle geom_ribbon element_text
#' @importFrom DESeq2 DESeqDataSetFromMatrix vst estimateSizeFactors
#' @importFrom MatrixGenerics rowVars rowMins
#' @importFrom SummarizedExperiment assay
#' @importFrom BiocGenerics counts
#' @importFrom ggdendro ggdendrogram
#' @importFrom stringr str_to_title
#' @importFrom stats prcomp na.omit formula pnorm t.test
#' @importFrom factoextra fviz_pca
#' @importFrom scran buildKNNGraph
#' @importFrom igraph modularity set_vertex_attr layout_nicely V
#' @importFrom patchwork wrap_plots plot_layout plot_annotation
#' @importFrom ggraph ggraph geom_edge_link geom_node_point
#' @importFrom metap sumlog sumz

#' @return Metrics of cell types pseudobulk clustering
#' @export get.cluster.score
#'


get.cluster.score <- function(hit.object = NULL,
                              cluster.by = NULL,
                              cluster.by.drop.na = TRUE,
                              batching = NULL,
                              scores = c("Silhouette_isolated", "Silhouette", "Modularity"),
                              modularity.k = 3,
                              dist.method = "euclidean",
                              hclust.method = "complete",

                              # For scores p-value calculation
                              ntests = 100, # number of shuffling events
                              seed = 22, # seed for random shuffling
                              pval.combine.method = "weighted_zmethod",

                              # For PCA
                              pca_comps_labs_invisible = c("quali"),
                              pca_pb_labs_invisible = c("var", "quali"),
                              pca_sig_labs_invisible = c("quali"),

                              # Pseudobulk params
                              ndim = 10,
                              nVarGenes = 500,
                              gene.filter = "HVG",
                              black.list = NULL,

                              ncores = round(parallelly::availableCores() / 3), # to reduce memory load
                              bparam = NULL,
                              progressbar = TRUE) {

  if (is.null(hit.object) ||
      !inherits(hit.object, "HiT")) {
    stop("Please provide a Hit class object or a count matrix.")
  }
  if (is.null(cluster.by)) {
    stop("Please provide a metadata column name to cluster by.")
  }

  if (!any(cluster.by %in% names(hit.object@metadata))) {
    stop("Group.by variables ", paste(cluster.by, collapse = ", "), " not found in metadata")
  } else if (!all(cluster.by %in% names(hit.object@metadata))) {
    cluster.by <- cluster.by[cluster.by %in% names(hit.object@metadata)]
    message("Only found ", paste(cluster.by, collapse = ", ") , " as grouping variable for HiT Object.")
  }

  # Need to replace special characters
  colnames(hit.object@metadata) <- make.names(colnames(hit.object@metadata))
  cluster.by <- make.names(cluster.by)

  for (i in cluster.by) {
    if (length(unique(hit.object@metadata[[i]])) == 1) {
      stop("All values are the same in cluster.by ", i, ". Please provide a metadata column with at least two different groups.")
    }
  }

  # Convert NA to factor level or drop
  for (i in cluster.by) {
    if (cluster.by.drop.na) {
      hit.object@metadata <- hit.object@metadata[!is.na(hit.object@metadata[[i]]), ]
    } else {
      hit.object@metadata[[i]] <- as.character(hit.object@metadata[[i]])
      hit.object@metadata[[i]][is.na(hit.object@metadata[[i]])] <- "NA"
    }
    # cluster.by column must be factor
    hit.object@metadata[[i]] <- as.factor(hit.object@metadata[[i]])
  }

  scores <- str_to_title(tolower(scores))

  # set parallelization parameters
  param <- set_parallel_params(ncores = ncores,
                               bparam = bparam,
                               progressbar = progressbar)

  # empty list to fill in the loop
  results <- list()


  # Process data ###############################################
  for (cluster_col in cluster.by) {
    message("Processing ", cluster_col)

    ## Process celltype composition ###############################################
    message("\nProcessing cell type composition\n")

    type <- "composition"

    comp_layers <- names(hit.object@composition)

    for (layer in comp_layers) {
      if (inherits(hit.object@composition[[layer]], "data.frame")) {
        mat <- hit.object@composition[[layer]][, c("celltype", "clr", "hitme.sample"), with = FALSE]
        if (cluster.by.drop.na) {
          mat <- mat %>% filter(hitme.sample %in% hit.object@metadata[["hitme.sample"]])
        }
        mat <- mat %>%
          tidyr::pivot_wider(names_from = hitme.sample,
                             values_from = clr) %>%
          stats::na.omit() %>%
          tibble::column_to_rownames(var = "celltype") %>%
          scale(center = TRUE,
                scale = FALSE)

        if (is.null(batching))  {
          cluster_labels <- hit.object@metadata %>%
            dplyr::filter(hitme.sample %in% colnames(mat)) %>%
            .[[cluster_col]]

          results[[cluster_col]][[type]][[layer]] <-
            get.scores(matrix = mat,
                       cluster_labels = cluster_labels,
                       scores = scores,
                       modularity.k = modularity.k,
                       dist.method = dist.method,
                       ntests = ntests,
                       seed = seed,
                       title = paste(cluster_col,
                                     stringr::str_to_title(type),
                                     layer),
                       invisible = pca_comps_labs_invisible)
        }

        if (!is.null(batching))  {
          for (b_var in batching) {
            if (b_var != cluster_col) {
              b_var_res_summary <- list()
              for (b in unique(hit.object@metadata[[b_var]])) {
                meta <- hit.object@metadata %>%
                  dplyr::filter(get(b_var) == b) %>%
                  dplyr::filter(hitme.sample %in% colnames(mat))

                cluster_labels <- meta[[cluster_col]]

                m <- mat[ , colnames(mat) %in% as.character(meta[["hitme.sample"]])] %>%
                  scale(center = TRUE,
                        scale = FALSE)

                results[[cluster_col]][[type]][[layer]][[b_var]][[b]] <-
                  get.scores(matrix = m,
                             cluster_labels = cluster_labels,
                             scores = scores,
                             modularity.k = modularity.k,
                             dist.method = dist.method,
                             ntests = ntests,
                             seed = seed,
                             title = paste(cluster_col,
                                           stringr::str_to_title(type),
                                           layer),
                             invisible = pca_comps_labs_invisible)
                for (score in scores) {
                  b_var_res_summary[[score]][["summary"]] <- c(
                    b_var_res_summary[[score]][["summary"]],
                    results[[cluster_col]][[type]][[layer]][[b_var]][[b]][["Scores"]][[score]][["summary"]])
                  b_var_res_summary[[score]][["n"]] <- c(
                    b_var_res_summary[[score]][["n"]],
                    results[[cluster_col]][[type]][[layer]][[b_var]][[b]][["Scores"]][[score]][["n"]])
                  b_var_res_summary[[score]][["p_value"]] <- c(
                    b_var_res_summary[[score]][["p_value"]],
                    results[[cluster_col]][[type]][[layer]][[b_var]][[b]][["Scores"]][[score]][["p_value"]])
                }
              }

              for (score in scores) {
                b_var_res_summary[[score]][["summary"]] <-
                  mean(b_var_res_summary[[score]][["summary"]])

                # Requires the number of samples per batch, so run before summing n
                p_values_not_na <- stats::na.omit(b_var_res_summary[[score]][["p_value"]])
                p_values_not_na_len <- length(p_values_not_na)
                if (p_values_not_na_len > 1){
                  b_var_res_summary[[score]] <- combine_pvals(b_var_res_summary[[score]],
                                                              pval.combine.method)
                } else if (p_values_not_na_len == 1) {
                  b_var_res_summary[[score]][["p_value"]] <- p_values_not_na
                } else {
                  b_var_res_summary[[score]][["p_value"]] <- NULL
                  b_var_res_summary[[score]][["summary"]] <- NULL
                }

                b_var_res_summary[[score]][["n"]] <-
                  sum(b_var_res_summary[[score]][["n"]])
              }
              results[[cluster_col]][[type]][[layer]][[b_var]][["all"]][["Scores"]] <- b_var_res_summary
            }
          }
        }

      } else if (is.list(hit.object@composition[[layer]])) {
        results[[cluster_col]][[type]][[layer]] <-
          BiocParallel::bplapply(
            X = names(hit.object@composition[[layer]]),
            BPPARAM = param,
            function(i){
              mat <- hit.object@composition[[layer]][[i]][, c("celltype", "clr", "hitme.sample"), with = F]
              if (cluster.by.drop.na) {
                mat <- mat %>% filter(hitme.sample %in% hit.object@metadata[["hitme.sample"]])
              }
              mat <- mat %>%
                tidyr::pivot_wider(names_from = hitme.sample,
                                   values_from = clr) %>%
                tibble::column_to_rownames(var = "celltype")

              mat <- scale(mat, center = TRUE,
                           scale = FALSE)

              if (is.null(batching))  {
                cluster_labels <- hit.object@metadata %>%
                  dplyr::filter(hitme.sample %in% colnames(mat)) %>%
                  .[[cluster_col]]

                if (nrow(mat) > 1) {
                  res <- get.scores(matrix = mat,
                                    cluster_labels = cluster_labels,
                                    scores = scores,
                                    modularity.k = modularity.k,
                                    dist.method = dist.method,
                                    ntests = ntests,
                                    seed = seed,
                                    title = paste(cluster_col,
                                                  stringr::str_to_title(type),
                                                  layer,
                                                  i),
                                    invisible = pca_comps_labs_invisible)
                  return(res)
                } else {
                  return(NULL)
                }
              }

              if (!is.null(batching)) {
                res <- list()
                for (b_var in batching) {
                  if (b_var != cluster_col) {
                    b_var_res_summary <- list()
                    for (b in unique(hit.object@metadata[[b_var]])) {
                      meta <- hit.object@metadata %>%
                        dplyr::filter(get(b_var) == b) %>%
                        dplyr::filter(hitme.sample %in% colnames(mat))

                      cluster_labels <- meta[[cluster_col]]

                      m <- mat[ , colnames(mat) %in% as.character(meta[["hitme.sample"]])] %>%
                        scale(center = TRUE, scale = FALSE)

                      if (nrow(m) > 1) {
                        res[[b_var]][[b]] <-
                          get.scores(matrix = m,
                                     cluster_labels = cluster_labels,
                                     scores = scores,
                                     modularity.k = modularity.k,
                                     dist.method = dist.method,
                                     ntests = ntests,
                                     seed = seed,
                                     title = paste(cluster_col,
                                                   stringr::str_to_title(type),
                                                   layer,
                                                   i),
                                     invisible = pca_comps_labs_invisible)
                      }

                      for (score in scores) {
                        b_var_res_summary[[score]][["summary"]] <- c(
                          b_var_res_summary[[score]][["summary"]],
                          res[[b_var]][[b]][["Scores"]][[score]][["summary"]])
                        b_var_res_summary[[score]][["n"]] <- c(
                          b_var_res_summary[[score]][["n"]],
                          res[[b_var]][[b]][["Scores"]][[score]][["n"]])
                        b_var_res_summary[[score]][["p_value"]] <- c(
                          b_var_res_summary[[score]][["p_value"]],
                          res[[b_var]][[b]][["Scores"]][[score]][["p_value"]])
                      }
                    }

                    for (score in scores) {
                      b_var_res_summary[[score]][["summary"]] <-
                        mean(b_var_res_summary[[score]][["summary"]])

                      # Requires the number of samples per batch, so run before summing n
                      p_values_not_na <- stats::na.omit(b_var_res_summary[[score]][["p_value"]])
                      p_values_not_na_len <- length(p_values_not_na)
                      if (p_values_not_na_len > 1) {
                        b_var_res_summary[[score]] <- combine_pvals(b_var_res_summary[[score]],
                                                                    pval.combine.method)
                      } else if (p_values_not_na_len == 1) {
                        b_var_res_summary[[score]][["p_value"]] <- p_values_not_na
                      } else {
                        b_var_res_summary[[score]][["p_value"]] <- NULL
                        b_var_res_summary[[score]][["summary"]] <- NULL
                      }

                      b_var_res_summary[[score]][["n"]] <-
                        sum(b_var_res_summary[[score]][["n"]])
                    }
                    res[[b_var]][["all"]][["Scores"]] <- b_var_res_summary
                  }
                }
                if(length(res) == 0) {
                  return(NULL)
                } else {
                  return(res)
                }
              }
            }
          )
        names(results[[cluster_col]][[type]][[layer]]) <-
          names(hit.object@composition[[layer]])
      }
    }


    ## Process pseudobulk ###############################################
    message("\nProcessing Pseudobulks\n")

    type <- "pseudobulk"

    pb_layers <- names(hit.object@aggregated_profile[[type]])

    for (layer in pb_layers) {
      results[[cluster_col]][[type]][[layer]] <- BiocParallel::bplapply(
        X = names(hit.object@aggregated_profile[[type]][[layer]]),
        BPPARAM = param,
        function(i){
          mat <- hit.object@aggregated_profile[[type]][[layer]][[i]]
          meta <- hit.object@metadata %>%
            dplyr::filter(hitme.sample %in% colnames(mat))
          cluster_labels <- meta[[cluster_col]]
          if (cluster.by.drop.na) {
            mat <- mat[, colnames(mat) %in% as.character(meta[["hitme.sample"]])]
          }

          if (is.null(batching))  {
            if (length(unique(cluster_labels)) > 1) {
              mat <- preproc_pseudobulk(matrix = mat,
                                        metadata = meta,
                                        cluster.by = cluster_col,
                                        nVarGenes = nVarGenes,
                                        gene.filter = gene.filter,
                                        black.list = black.list)

              res <- get.scores(matrix = mat,
                                cluster_labels = cluster_labels,
                                scores = scores,
                                modularity.k = modularity.k,
                                dist.method = dist.method,
                                ntests = ntests,
                                seed = seed,
                                title = paste(cluster_col,
                                              stringr::str_to_title(type),
                                              layer,
                                              i),
                                invisible = pca_pb_labs_invisible)
              return(res)
            } else {
              return(NULL)
            }
          }

          if (!is.null(batching)) {
            res <- list()
            for (b_var in batching) {
              if (b_var != cluster_col) {
                b_var_res_summary <- list()
                for (b in unique(hit.object@metadata[[b_var]])) {
                  met <- hit.object@metadata %>%
                    dplyr::filter(get(b_var) == b) %>%
                    dplyr::filter(hitme.sample %in% colnames(mat))

                  cluster_labels <- met[[cluster_col]]

                  m <- mat[ , colnames(mat) %in% met[["hitme.sample"]]]

                  if (length(unique(cluster_labels)) > 1) {
                    m <- preproc_pseudobulk(matrix = m,
                                            metadata = met,
                                            cluster.by = cluster_col,
                                            nVarGenes = nVarGenes,
                                            gene.filter = gene.filter,
                                            black.list = black.list)

                    if (nrow(m) > 1) {
                      res[[b_var]][[b]] <-
                        get.scores(matrix = m,
                                   cluster_labels = cluster_labels,
                                   scores = scores,
                                   modularity.k = modularity.k,
                                   dist.method = dist.method,
                                   ntests = ntests,
                                   seed = seed,
                                   title = paste(cluster_col,
                                                 stringr::str_to_title(type),
                                                 layer,
                                                 i),
                                   invisible = pca_pb_labs_invisible)
                    }
                  }
                  for (score in scores) {
                    b_var_res_summary[[score]][["summary"]] <- c(
                      b_var_res_summary[[score]][["summary"]],
                      res[[b_var]][[b]][["Scores"]][[score]][["summary"]])
                    b_var_res_summary[[score]][["n"]] <- c(
                      b_var_res_summary[[score]][["n"]],
                      res[[b_var]][[b]][["Scores"]][[score]][["n"]])
                    b_var_res_summary[[score]][["p_value"]] <- c(
                      b_var_res_summary[[score]][["p_value"]],
                      res[[b_var]][[b]][["Scores"]][[score]][["p_value"]])
                  }
                }

                for (score in scores) {
                  b_var_res_summary[[score]][["summary"]] <-
                    mean(b_var_res_summary[[score]][["summary"]])

                  # Requires the number of samples per batch, so run before summing n
                  p_values_not_na <- stats::na.omit(b_var_res_summary[[score]][["p_value"]])
                  p_values_not_na_len <- length(p_values_not_na)
                  if (p_values_not_na_len > 1) {
                    b_var_res_summary[[score]] <- combine_pvals(b_var_res_summary[[score]],
                                                                pval.combine.method)
                  } else if (p_values_not_na_len == 1) {
                    b_var_res_summary[[score]][["p_value"]] <- p_values_not_na
                  } else {
                    b_var_res_summary[[score]][["p_value"]] <- NULL
                    b_var_res_summary[[score]][["summary"]] <- NULL
                  }

                  b_var_res_summary[[score]][["n"]] <-
                    sum(b_var_res_summary[[score]][["n"]])
                }
                res[[b_var]][["all"]][["Scores"]] <- b_var_res_summary
              }
            }
            if(length(res) == 0) {
              return(NULL)
            } else {
              return(res)
            }
          }
        }
      )
      names(results[[cluster_col]][[type]][[layer]]) <-
        names(hit.object@aggregated_profile[[type]][[layer]])
    }


    ## Process signatures ###############################################
    message("\nProcessing Signatures\n")

    type <- "signatures"

    comp_layers <- names(hit.object@aggregated_profile[[type]])

    if (!is.null(comp_layers)) {
      for (layer in comp_layers) {
        cols <- colnames(hit.object@aggregated_profile[[type]][[layer]])
        signatures <- cols[! cols %in% c("celltype", "hitme.sample")]

        results[[cluster_col]][[type]][[layer]] <- BiocParallel::bplapply(
          X = signatures,
          BPPARAM = param,
          function(i){
            mat <- hit.object@aggregated_profile[[type]][[layer]][, c("celltype", i, "hitme.sample"), with = FALSE]
            if (cluster.by.drop.na) {
              mat <- mat %>% filter(hitme.sample %in% hit.object@metadata[["hitme.sample"]])
            }
            mat <- mat %>%
              tidyr::pivot_wider(names_from = hitme.sample, values_from = i) %>%
              tibble::column_to_rownames(var = "celltype") %>%
              replace(is.na(.), 0)

            if (is.null(batching))  {
              mat <- scale(mat,
                           center = TRUE,
                           scale = TRUE)

              # Drop columns containing NAN, caused by scaling zero variance columns
              mat <- mat[ , colSums(is.nan(mat)) == 0]

              cluster_labels <- hit.object@metadata %>%
                filter(hitme.sample %in% colnames(mat)) %>%
                .[[cluster_col]]

              res <- get.scores(matrix = mat,
                                cluster_labels = cluster_labels,
                                scores = scores,
                                modularity.k = modularity.k,
                                dist.method = dist.method,
                                ntests = ntests,
                                seed = seed,
                                title = paste(cluster_col,
                                              stringr::str_to_title(type),
                                              layer,
                                              i),
                                invisible = pca_sig_labs_invisible)
              return(res)
            }

            if (!is.null(batching)) {
              res <- list()
              for (b_var in batching) {
                if (b_var != cluster_col) {
                  b_var_res_summary <- list()
                  for (b in unique(hit.object@metadata[[b_var]])) {
                    meta <- hit.object@metadata %>%
                      dplyr::filter(get(b_var) == b) %>%
                      dplyr::filter(hitme.sample %in% colnames(mat))

                    cluster_labels <- meta[[cluster_col]]

                    m <- mat[ , colnames(mat) %in% as.character(meta[["hitme.sample"]])] %>%
                      scale(center = TRUE,
                            scale = TRUE)

                    # Drop columns containing NAN, caused by scaling zero variance columns
                    m <- m[ , colSums(is.nan(m)) == 0]

                    res[[b_var]][[b]] <-
                      get.scores(matrix = m,
                                 cluster_labels = cluster_labels,
                                 scores = scores,
                                 modularity.k = modularity.k,
                                 dist.method = dist.method,
                                 ntests = ntests,
                                 seed = seed,
                                 title = paste(cluster_col,
                                               stringr::str_to_title(type),
                                               layer,
                                               i),
                                 invisible = pca_sig_labs_invisible)
                    for (score in scores) {
                      b_var_res_summary[[score]][["summary"]] <- c(
                        b_var_res_summary[[score]][["summary"]],
                        res[[b_var]][[b]][["Scores"]][[score]][["summary"]])
                      b_var_res_summary[[score]][["n"]] <- c(
                        b_var_res_summary[[score]][["n"]],
                        res[[b_var]][[b]][["Scores"]][[score]][["n"]])
                      b_var_res_summary[[score]][["p_value"]] <- c(
                        b_var_res_summary[[score]][["p_value"]],
                        res[[b_var]][[b]][["Scores"]][[score]][["p_value"]])
                    }
                  }

                  for (score in scores) {
                    b_var_res_summary[[score]][["summary"]] <-
                      mean(b_var_res_summary[[score]][["summary"]])

                    # Requires the number of samples per batch, so run before summing n
                    p_values_not_na <- stats::na.omit(b_var_res_summary[[score]][["p_value"]])
                    p_values_not_na_len <- length(p_values_not_na)
                    if (p_values_not_na_len > 1) {
                      b_var_res_summary[[score]] <- combine_pvals(b_var_res_summary[[score]],
                                                                  pval.combine.method)
                    } else if (p_values_not_na_len == 1) {
                      b_var_res_summary[[score]][["p_value"]] <- p_values_not_na
                    } else {
                      b_var_res_summary[[score]][["p_value"]] <- NULL
                      b_var_res_summary[[score]][["summary"]] <- NULL
                    }

                    b_var_res_summary[[score]][["n"]] <-
                      sum(b_var_res_summary[[score]][["n"]])
                  }
                  res[[b_var]][["all"]][["Scores"]] <- b_var_res_summary
                }
              }
              if(length(res) == 0) {
                return(NULL)
              } else {
                return(res)
              }
            }
          }
        )
        names(results[[cluster_col]][[type]][[layer]]) <- signatures
      }
    }
  }

  # Save user set parameters for summarize.cluster.scores
  results[["params"]][["cluster.by"]] <- cluster.by
  results[["params"]][["batching"]] <- batching
  results[["params"]][["scores"]] <- scores

  return(results)
}




#' Summarize scores and plot heat map
#'
#' @param data Output from get.cluster.scores
#' @param topN Integer indicating number of topN highest scoring (most discriminating) features ranked for each score
#' @param create_plots Boolean indicating whether to create and show plots or not
#' @param p.adjustment Whether to adjust p-value columns or not
#' @param p.adjust.method Method for adjusting p-values (see stats::p.adjust for methods)
#' @param p.value_cutoff p-value (mean of all p-value columns) cutoff to filter out non-significant results
#'
#' @importFrom dplyr mutate filter
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
#' @importFrom stats p.adjust na.omit quantile
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom pheatmap pheatmap
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom ggplotify as.ggplot
#'
#' @return Average silhouette widths per grouping variable. Optionally, a heatmap plot for visualization
#' @export summarize.cluster.scores
#'


summarize.cluster.scores <- function(data = NULL,
                                     topN = NULL,
                                     create_plots = TRUE,
                                     p.adjustment = TRUE,
                                     p.adjust.method = "fdr",
                                     p.value_cutoff = 0.05) {

  show.variables <- c(".summary", ".n", ".p_value")

  if (is.null(data)) {
    message("Please provide scores object (output from get.cluster.scores)")
  }

  cluster.by <- data[["params"]][["cluster.by"]]
  batching <- data[["params"]][["batching"]]
  scores <- data[["params"]][["scores"]]

  data[["params"]] <- NULL # Remove params

  data_conts <- unlist(data)
  df_pars <- list()

  for (par in show.variables) {
    data_conts_temp <- data_conts[endsWith(names(data_conts), par)]

    if (!is.null(batching)) {
      data_conts_temp <- data_conts_temp[grepl(".all.Scores.",
                                               names(data_conts_temp))]
      names(data_conts_temp) <- gsub("all.", "", names(data_conts_temp))
    }

    if (length(data_conts_temp) == 0) {
      stop("Error happened at ", show.variables, ". No values left after filtering.")
    } else {
      target_string <- paste0(scores, collapse = "|")
      target_string <- paste0("(", target_string, "+)")
      data_conts_names_split <- strsplit(gsub(target_string,"~\\1",
                                              names(data_conts_temp)), "~")

      df <- data.frame(t(do.call(cbind, data_conts_names_split))) %>%
        dplyr::mutate(X3 = sapply(data_conts_temp, "[[", 1)) %>%
        tidyr::pivot_wider(names_from = "X1",
                           values_from = "X3") %>%
        tibble::column_to_rownames(var = "X2") %>%
        t() %>%
        as.data.frame()

      row.names(df) <- gsub(".Scores.",
                            "",
                            row.names(df))

      df_pars[[par]] <- df
    }
  }

  # Check if all columns with sample numbers are equal
  df <- as.matrix(df_pars[[".n"]])
  all_n_cols_equal <- all(apply(df, 2, identical, df[,1]))
  if (all_n_cols_equal) {
    df_pars[[".n"]] <- df_pars[[".n"]][, 1, drop = FALSE]
    colnames(df_pars[[".n"]]) <- "n"
  }

  df <- merge(df_pars[[show.variables[1]]],
              df_pars[[show.variables[2]]],
              by = "row.names",
              all = TRUE)
  row.names(df) <- df$Row.names
  df$Row.names <- NULL

  df <- merge(df,
              df_pars[[show.variables[3]]],
              by = "row.names",
              all = TRUE)
  row.names(df) <- df$Row.names
  df$Row.names <- NULL

  df_cluster.by_list <- list()
  for (c in cluster.by) {
    df_cluster.by_list[[c]] <- df %>%
      dplyr::filter(row.names(.) %>%
               startsWith(c))
    row.names(df_cluster.by_list[[c]]) <- gsub(paste0(c, "."),
                                               "",
                                               row.names(df_cluster.by_list[[c]]))
    # Remove batching variable name
    if (!is.null(batching)) {
      batch_names <- batching[!batching %in% c]
      batch_names <- batch_names %>%
        paste0(".", .) %>%
        paste(collapse = "|")
      row.names(df_cluster.by_list[[c]]) <- gsub(batch_names,
                                                 "",
                                                 row.names(df_cluster.by_list[[c]]))
    }

    # Adjust p-values and filter
    if (p.adjustment) {
      p_val_cols <- which(grepl(".p_value",
                                colnames(df_cluster.by_list[[c]])))
      for (i in p_val_cols){
        df_cluster.by_list[[c]][, i] <-
          stats::p.adjust(df_cluster.by_list[[c]][, i],
                          method = p.adjust.method)
      }

      df_cluster.by_list[[c]] <-
        df_cluster.by_list[[c]][rowMeans(df_cluster.by_list[[c]][, p_val_cols]) <= p.value_cutoff, ]

      df_cluster.by_list[[c]] <- stats::na.omit(df_cluster.by_list[[c]])

      if (nrow(df_cluster.by_list[[c]]) == 0) {
        df_cluster.by_list[[c]] <- NULL
        message(paste("For ", c, " no separation was found after p-value cutoff. You can try to set it higher."))

        next
      }
    }

    # Get topN scored results
    if (!is.null(topN) &&
        is.numeric(topN) &&
        length(topN) == 1) {
      summary_score_cols <- which(grepl(".summary",
                                        colnames(df_cluster.by_list[[c]])))
      topN_rownames <- c()
      for (i in summary_score_cols) {
        topN_rownames_i <- row.names(df_cluster.by_list[[c]])[
          order(df_cluster.by_list[[c]][, i], decreasing = TRUE)][1:topN]
        topN_rownames <- c(topN_rownames,
                           topN_rownames_i)
      }
      df_cluster.by_list[[c]] <- df_cluster.by_list[[c]][unique(topN_rownames), ]
    }
  }

  # Remove NULL elements
  df_cluster.by_list <- Filter(Negate(is.null), df_cluster.by_list)

  # Check if all df_cluster.by_list were NULL
  if (length(df_cluster.by_list) == 0) {

    message("No significant separation found for cluster.by provided")

  } else {

    if (create_plots) {
      df_list <- list()
      plot_list <- list()

      for (n in names(df_cluster.by_list)) {
        df <- stats::na.omit(df_cluster.by_list[[n]])

        n_breaks <- min(100, dim(df)[1] * dim(df)[2])
        quantiles <- seq(0, 1, 1/n_breaks)

        # Check that variance is not zero
        nonzero_var_cols <- unlist(lapply(df, function(x) !length(unique(x))==1))
        breaks <- df[, nonzero_var_cols] %>%
          scale() %>%
          stats::quantile(., quantiles) %>%
          unique()

        color_breaks <- length(breaks)

        color <- grDevices::colorRampPalette(
          rev(RColorBrewer::brewer.pal(n = 7,
                                       name = "RdYlBu")))(color_breaks)

        if (nrow(df) > 1) {
          scale <- "column"
        } else {
          scale <- "none"
        }

        df_pval_invert <- df
        df_pval <- df_pval_invert[, grepl(".p_value", names(df))]
        df_pval <- 1 / (df_pval + 1e-16)
        df_pval_invert[, grepl(".p_value", names(df))] <- df_pval

        plot_list[[n]] <- pheatmap::pheatmap(df_pval_invert,
                                             main = n,
                                             angle_col = 315,
                                             scale = scale,
                                             display_numbers = round(df, 3),
                                             number_color = "black",
                                             color = color,
                                             breaks = breaks,
                                             legend_breaks = 0,
                                             legend_labels = "",
                                             cluster_cols = FALSE,
                                             cluster_rows = FALSE,
                                             silent = TRUE)[[4]]
      }
      df_cluster.by_list[["plots"]][["plot_list"]] <- plot_list

      g <- gridExtra::grid.arrange(
        gridExtra::arrangeGrob(grobs = plot_list,
                               ncol=length(plot_list))
      )
      df_cluster.by_list[["plots"]][["summary_plot"]] <- ggplotify::as.ggplot(g)
    }

    return(df_cluster.by_list)
  }
}




#' Render plots summarizing celltype proportions and distribution in samples
#'
#'
#' @param obj.list List of Seurat objects
#' @param annot.col Metadata column(s) containing the cell type annotations
#' @param bottom.mar Adjustable bottom margin for long sample names

#' @importFrom stats setNames

#' @return Get percentage of not annotated cells per sample and plot it.
#' @export nas.per.sample
#'

nas.per.sample <- function (obj.list = NULL,
                            annot.col = c("scGate_multi"),
                            return.plot = TRUE,
                            bottom.mar = 10.2) {
  if (is.null(obj.list) &
      !is.list(obj.list) &
      !all(lapply(obj.list, inherits, "Seurat"))) {
    stop("Please provide a list of seurat objects")
  }

  na_list <- list()

  for (col in annot.col) {
    na_perc_per_sample <- c()
    for (i in names(obj.list)) {
      if (col %in% names(obj.list[[i]]@meta.data)) {
        percs <- prop.table(table(obj.list[[i]]@meta.data[[col]], useNA = "ifany"))*100
        nas_perc <- unname(percs[is.na(names(percs))])
        na_perc_per_sample <- c(na_perc_per_sample, stats::setNames(nas_perc, i))
      } else {
        stop(paste(col, " not found in obj.list item ", i))
      }
    }
    par(mar = c(bottom.mar, 4.1, 4.1, 2.1))
    if (return.plot) {
      barplot(na_perc_per_sample,
              main = col,
              ylab = "Percentage of NA values",
              las=2)
    }
    na_list[[col]] <- na_perc_per_sample
  }
  return(na_list)
}




#' Render bar plots summarizing celltype proportions and distribution in samples and groups
#'
#'
#' @param hit.object A Hit class object (typically after applying merge.HiTObjects onto a list of HiTObjects)
#' @param layer Default "layer1" if you have one cell type annotation layer in your hit.object. Alternatively "layer2" etc. if you have multiple layers of annotation depths.
#' @param return.plot.to.var Optionally, you can save the ggplots to a variable if you would like to further modify and adapt the plots on your own.
#' @param facet.by This allows you to pass a metadata column name present in your hit.object$metadata to show your samples in facets with ggplot facet_grid, for example by "condition".

#' @importFrom ggplot2 ggplot aes geom_bar theme element_text ggtitle facet_grid
#' @importFrom stats reformulate
#' @importFrom patchwork wrap_plots

#' @return Plotting function to show the cell type composition from HiTME object across different samples.
#' @export composition.barplot
#'

composition.barplot <- function (hit.object = NULL,
                                 sample.col = NULL,
                                 layer = "layer1",
                                 return.plot.to.var = FALSE,
                                 facet.by = NULL) {

  # Need to replace special characters
  colnames(hit.object@metadata) <- make.names(colnames(hit.object@metadata))

  sample.col <- "hitme.sample"

  if (is.null(hit.object)) {
    stop("Please provide input hit.object")
  }
  if (!length(layer) == 1 || !is.character(layer)) {
    stop("Please provide one character string for layer parameter")
  }
  if (!is.null(facet.by)) {
    if (!is.character(facet.by)) {
      stop("Please provide a character string or a vector of character strings for the facet.by parameter")
    }
    facet.by <- make.names(facet.by)
    facet.by.in.colnames <- facet.by %in% names(hit.object@metadata)
    if (!all(facet.by.in.colnames)) {
      facet.by.not.in.colnames <- facet.by[!facet.by.in.colnames]
      stop(paste("facet.by ", facet.by.not.in.colnames, " not found in hit.object@metadata column names"))
    }
  }


  comps <- hit.object@composition[[layer]]
  meta <- hit.object@metadata

  if (!is.null(facet.by)) {
    facet.by_reformulate <- reformulate(facet.by)
  }

  if (is.data.frame(comps)) {
    comp <- merge(comps, meta[, c(sample.col, facet.by), drop=FALSE], by = sample.col)

    p <- ggplot(comp, aes(x = hitme.sample, y = freq, fill = celltype)) +
      geom_bar(stat = "identity") +
      theme(axis.text.x = element_text(angle = 45, hjust=1))

    if (!is.null(facet.by)) {

      p <- p + facet_grid(facet.by_reformulate,
                          space  = "free",
                          scales = "free")
    }

    if (return.plot.to.var) {
      return(p)
    } else {
      print(p)
    }
  }

  else {
    p_list <- list()
    for (ct in names(comps)) {
      comp <- hit.object@composition[[layer]][[ct]]
      comp <- merge(comp, meta[, c(sample.col, facet.by), drop=FALSE], by = sample.col)

      p_list[["plot_list"]][[ct]] <- ggplot(comp, aes(x = hitme.sample, y = freq, fill = celltype)) +
        geom_bar(stat = "identity") +
        ggtitle(ct) +
        theme(axis.text.x = element_text(angle = 45, hjust=1)) +
        ggtitle(ct)

      if (!is.null(facet.by)) {
        p_list[["plot_list"]][[ct]] <- p_list[["plot_list"]][[ct]] +
          facet_grid(facet.by_reformulate,
                     space  = "free",
                     scales = "free")
      }

    }
    p_list[["arranged_plots"]] <- patchwork::wrap_plots(p_list[["plot_list"]])

    if (return.plot.to.var) {
      return(p_list)
    } else {
      print(p_list[["arranged_plots"]])
    }
  }
}


#' Render box plots summarizing celltype proportions and distribution in samples and groups
#'
#'
#' @param hit.object A Hit class object (typically after applying merge.HiTObjects onto a list of HiTObjects)
#' @param plot.var Column in the hit.object$composition: either "freq" for cell type relative abundance in percent or "clr" (for centered log-ratio transformed). Default: "clr" as it is better suited for statistical analysis and is better able to also show low abundant cell types.
#' @param layer Default "layer1" if you have one cell type annotation layer in your hit.object. Alternatively "layer2" etc. if you have multiple layers of annotation depths.
#' @param return.plot.to.var Optionally, you can save the ggplots to a variable if you would like to further modify and adapt the plots on your own.
#' @param group.by This allows you to pass a metadata column name present in your hit.object$metadata to show your samples in groups, for example by "condition".
#' @param facet.by This allows you to pass a metadata column name present in your hit.object$metadata to show your samples in facets with ggplot facet_grid, for example by "condition".
#' @param pval.method Specify method how to calculate the p-value. Wilcoxon test is recommended as compositional data might not fullfill the assumption of a gaussian distribution. For alternatives, see documentation of ggpubr::stat_pwc.
#' @param p.adjust.method Method for adjusting p-values (see ggpubr::stat_pwc for available methods)
#' @param palette Choose a palette of your liking. For available palettes, see ggsci package. Default: "lancet"
#' @param legend.position Where to put the legend. Possible options: "top", "right", "bottom", "left"

#' @importFrom ggplot2 ggplot aes geom_boxplot theme element_text ggtitle facet_grid position_jitterdodge
#' @importFrom ggpubr stat_pwc ggboxplot
#' @importFrom stats reformulate
#' @importFrom patchwork wrap_plots

#' @return Plotting function to show the cell type composition from HiTME object across different samples.
#' @export composition.boxplot
#'

composition.boxplot <- function (hit.object = NULL,
                                 plot.var = "clr",
                                 layer = "layer1",
                                 return.plot.to.var = FALSE,
                                 group.by = NULL,
                                 facet.by = NULL,
                                 pval.method = "wilcox.test",
                                 p.adjust.method = "BH",
                                 palette = "lancet",
                                 legend.position = "right") {

  # Need to replace special characters
  colnames(hit.object@metadata) <- make.names(colnames(hit.object@metadata))

  sample.col <- "hitme.sample"

  if (is.null(hit.object)) {
    stop("Please provide input hit.object")
  }
  if (!length(plot.var) == 1 ||
      !is.character(plot.var) ||
      !plot.var %in% c("freq", "clr")) {
    stop("Please provide one character string for plot.var, either 'freq' or 'clr'")
  }
  if (!length(layer) == 1 || !is.character(layer)) {
    stop("Please provide one character string for layer parameter")
  }
  if (!is.null(group.by)) {
    if (!length(group.by) == 1 || !is.character(group.by)) {
      stop("Please provide one character string for the group.by parameter")
    }
    group.by <- make.names(group.by)
    group.by.gg <- sym(group.by)
    nr_of_boxplots <- hit.object@metadata[[group.by]] %>%
      unique() %>%
      length() %>%
      "*"(1.5) %>%
      round()
    hit.object@metadata[group.by] <- lapply(hit.object@metadata[group.by], as.factor)
  }
  if (!is.null(facet.by)) {
    if (!is.character(facet.by)) {
      stop("Please provide a character string or a vector of character strings for the facet.by parameter")
    }
    facet.by <- make.names(facet.by)

    facet.by.in.colnames <- facet.by %in% names(hit.object@metadata)
    if (!all(facet.by.in.colnames)) {
      facet.by.not.in.colnames <- facet.by[!facet.by.in.colnames]
      stop(paste("facet.by ", facet.by.not.in.colnames, " not found in hit.object@metadata column names"))
    }
  }

  comps <- hit.object@composition[[layer]]
  meta <- hit.object@metadata

  plot.var.gg <- sym(plot.var)

  if (is.data.frame(comps)) {
    comp <- merge(comps, meta[, c(sample.col, group.by, facet.by), drop=FALSE], by = sample.col)

    # Need to check if group.by is NULL
    # Due to a presumed bug, if group.by is passed as variable to ggboxplot, even if it is assigned NULL, it throws an error
    if (is.null(group.by)) {
      p <- ggboxplot(comp,
                     x = "celltype",
                     y = plot.var,
                     outlier.shape = NA,
                     palette = palette,
                     facet.by = facet.by,
                     legend = legend.position) +
        geom_jitter(width = 0.2, size = 1)
    } else {
      p <- ggboxplot(comp,
                     x = "celltype",
                     y = plot.var,
                     color = group.by,
                     outlier.shape = NA,
                     palette = palette,
                     facet.by = facet.by,
                     legend = legend.position) +
        geom_jitter(mapping = aes(color = !!group.by.gg), position=position_jitterdodge(jitter.width = 1/nr_of_boxplots), size = 1) +
        stat_pwc(aes(group = !!group.by.gg),
                 label = "p.signif",
                 method = pval.method,
                 p.adjust.method = p.adjust.method,
                 p.adjust.by = "panel",
                 tip.length = 0,
                 hide.ns = TRUE)
    }

    p <- p +
      theme(axis.text.x = element_text(angle = 45, hjust=1)) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

    if (return.plot.to.var) {
      return(p)
    } else {
      print(p)
    }
  }

  else {
    p_list <- list()
    for (ct in names(comps)) {
      comp <- hit.object@composition[[layer]][[ct]]
      comp <- merge(comp, meta[, c(sample.col, group.by, facet.by), drop=FALSE], by = sample.col)

      # Need to check if group.by is NULL
      # Due to a presumed bug, if group.by is passed as variable to ggboxplot, even if it is assigned NULL, it throws an error
      if (is.null(group.by)) {
        p_list[["plot_list"]][[ct]] <- ggboxplot(comp,
                                                 x = "celltype",
                                                 y = plot.var,
                                                 outlier.shape = NA,
                                                 palette = palette,
                                                 facet.by = facet.by,
                                                 legend = legend.position) +
          geom_jitter(width = 0.2, size = 1)
      } else {
        p_list[["plot_list"]][[ct]] <- ggboxplot(comp,
                                                 x = "celltype",
                                                 y = plot.var,
                                                 color = group.by,
                                                 outlier.shape = NA,
                                                 palette = palette,
                                                 facet.by = facet.by,
                                                 legend = legend.position) +
          geom_jitter(mapping = aes(color = !!group.by.gg), position=position_jitterdodge(jitter.width = 1/nr_of_boxplots), size = 1) +
          stat_pwc(aes(group = !!group.by.gg),
                   label = "p.signif",
                   method = pval.method,
                   p.adjust.method = p.adjust.method,
                   p.adjust.by = "panel",
                   tip.length = 0,
                   hide.ns = TRUE)
      }

      p_list[["plot_list"]][[ct]] <- p_list[["plot_list"]][[ct]] +
        theme(axis.text.x = element_text(angle = 45, hjust=1)) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

    }
    p_list[["arranged_plots"]] <- patchwork::wrap_plots(p_list[["plot_list"]])

    if (return.plot.to.var) {
      return(p_list)
    } else {
      print(p_list[["arranged_plots"]])
    }
  }
}



#' Build confusion matrix-like plots between different cell type classification approaches
#'
#'
#' @param hit.object List of Hit class object
#' @param var.1 1st of grouping variables on the x-axis. If using \code{relative = TRUE} proportions will be normalized to this variable.
#' @param var.2 2nd of grouping variables on the x-axis.
#' @param relative Whether to show absolute number of cells, or relative number cells out of the 1st variable indicted. Default is FALSE.
#' @param useNA Whether to include not annotated cells or not (labelled as "NA"). Can be "no", "ifany", or "always". See \code{?table} for more information. Default is "ifany".
#' @param type Type of plot to render. Either "tile" (default) as a confusion matrices or "barplot".

#' @importFrom dplyr mutate group_by distinct
#' @importFrom ggplot2 aes geom_point geom_tile geom_col scale_fill_gradient labs geom_label theme
#' @importFrom cowplot theme_cowplot
#' @importFrom data.table rbindlist

#' @return Plots to evaluate the correspondence between different classification methods.
#' @export plot.confusion.matrix
#'

plot.confusion.matrix <- function(hit.object = NULL,
                                  var.1 = NULL,
                                  var.2 = NULL,
                                  relative = FALSE,
                                  useNA = "ifany",
                                  type = "tile") {

  if (is.null(hit.object)) {
    stop("Please provide a single one or a list of HiT object")
  }

  if (!is.list(hit.object)) {
    hit.object <- list(hit.object)
  }

  if (suppressWarnings(!all(lapply(hit.object, function(x) {inherits(x, "HiT")})))) {
    stop("Not all components of the list are HiT objects.")
  }

  #give name to list of hit objects
  for (v in seq_along(hit.object)) {
    if (is.null(names(hit.object)[[v]]) ||
        is.na(names(hit.object)[[v]])) {
      names(hit.object)[[v]] <- paste0("Sample", v)
    }
  }

  # join vars to plot
  vars <- c(var.1, var.2)

  if (any(is.null(vars))) {
    stop("Please provide 2 cell type classification labels common in all elements of the list of Hit objects: var.1 and var.2")
  }

  if (suppressWarnings(!all(lapply(hit.object,
                                   function(x) {
                                     any(vars %in% names(x@metadata))
                                   })))) {
    stop("Not all supplied HiT object contain ",
         paste(vars, collapse = ", "),
         " group.by elements in their metadata")
  }


  # build all table
  data <- lapply(names(hit.object), function(y) {
    sel <- names(hit.object[[y]]@metadata) %in% vars
    a <- hit.object[[y]]@metadata[,sel, drop = FALSE] %>%
      dplyr::mutate(hitme.sample = y)
    return(a)
  }) %>%
    data.table::rbindlist(fill = TRUE) %>%
    as.data.frame()

  # Get counts for every group
  pz <- table(data[[var.1]],
              data[[var.2]],
              useNA = useNA)
  if (relative) {
    pz <- pz %>% prop.table(margin = 1)
    legend <- "Relative counts"
  } else {
    legend <- "Absolute counts"
  }


  if (tolower(type) == "tile") {
    plot <- pz %>% as.data.frame() %>%
      ggplot2::ggplot(ggplot2::aes(Var1, Var2,
                                   fill = Freq)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient(name = legend,
                                   low = "#FFFFC8",
                                   high = "#7D0025") +
      ggplot2::labs(x = var.1,
                    y = var.2) +
      ggplot2::geom_label(ggplot2::aes(label = ifelse(Freq > 0,
                                                      round(Freq, 1),
                                                      NA)),
                          color = "white",
                          alpha = 0.6,
                          fill = "black") +
      cowplot::theme_cowplot()
  } else if (grepl("bar|col", type, ignore.case = TRUE)) {
    plot <- pz %>% as.data.frame() %>%
      ggplot2::ggplot(ggplot2::aes(Var1, Freq,
                                   fill = Var2)) +
      ggplot2::labs(x = var.1,
                    y = legend,
                    fill = var.2) +
      ggplot2::geom_col() +
      ggplot2::theme_bw()
  }

  plot <- plot +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                       hjust = 1,
                                                       vjust = 1))
  return(plot)
}




#' Plot gene expression along scGate model
#'
#'
#' @param object Seurat object or list of them, if list is provided one plot per Seurat object is returned
#' @param scGate.model scGate model to plot their gating genes
#' @param group.by Whether to group gene expression by a variable in metadata
#' @param split.by Whether to generate one plot for each grouping variable provided

#' @import scGate
#' @importFrom Seurat VlnPlot
#' @importFrom ggplot2 ylab ggtitle ggplot theme_void
#' @importFrom ggpubr ggarrange annotate_figure

#' @return Plots to evaluate the gene expression of scGate models
#' @export plot_gene_gating
#'

plot_gene_gating <- function(object = NULL,
                             scGate.model = NULL,
                             group.by = NULL,
                             split.by = NULL) {

  if (is.null(object)) {
    stop("Please provide a Seurat object or a list of them")
  }

  if (is.list(object) && !is.null(split.by)) {
    stop("split.by only supported for a single Seurat object, not a list.\n
         Merge list before running HiTME")
  }

  if (is.null(scGate.model)) {
    stop("Please provide a scGate model or list of them")
  }

  if (!is.list(scGate.model)) {
    scGate.model <- list("scGate_model" = scGate.model)
  }

  # split object into a list if indicated
  if (!is.null(split.by)) {
    if (is.list(object)) {
      stop("Split.by argument not supported when providing a list of
           Seurat objects. Set split.by = NULL or merge list.")
    }
    if (!split.by %in% names(object@meta.data)) {
      stop(paste("split.by argument: ",
                 split.by,
                 " is not a metadata column in this Seurat object"))
    }
    object <- Seurat::SplitObject(object, split.by = split.by)

  }

  if (!is.null(group.by)) {
    if (!group.by %in% names(object@meta.data)) {
      stop(paste("group.by argument: ",
                 group.by,
                 " is not a metadata column in this Seurat object"))
    }
  } else {
    group.by <- "orig.ident"
  }

  # if object is unique turn into a list
  if (!is.list(object)) {
    object <- list("object" = object)
  }

  plots <- list()
  # render plots
  model <- scGate:::table.to.model(scGate.model)

  suppressWarnings({
    for (ob in names(object)) {
      if (!inherits(object[[ob]], "Seurat")) {
        stop("Not Seurat object included, cannot be processed.\n")
      }

      # keep only genes expressed
      sc_names <- row.names(object[[ob]])[rowSums(object[[ob]]) > 0]
      # make plot list for each level
      pl.list <- list()

      for (a in names(model)) {
        m <- model[[a]] %>% unlist(recursive = FALSE)
        pl.sublist <- list()
        for (e in names(m)) {
          # remove - sign on negative markers
          feat <- m[[e]] %>% gsub("-", "", .)

          feat <- intersect(feat, sc_names)

          if (length(feat)>0) {
            # do not stack if only one gene is present
            stack <- ifelse(length(feat) > 1,
                            TRUE,
                            FALSE)

            pl.sublist[[e]] <-
              Seurat::VlnPlot(object[[ob]],
                              features = feat,
                              group.by = group.by,
                              stack = stack,
                              pt.size = 0,
                              flip = TRUE) +
              ggplot2::ggtitle(e) +
              Seurat::NoLegend() +
              ggplot2::xlab("") +
              {if (length(feat) == 1) {
                ggplot2::ylab(feat)
              }} +
              ggplot2::theme(plot.title = element_text(hjust = 0.5))

          }
        }

        pl.list[[a]] <- pl.sublist
      }

      # max number of plots
      max <- lapply(pl.list, length) %>%
        unlist() %>%
        max()

      for (p in names(pl.list)) {

        # add blank plots if needed
        while(length(pl.list[[p]])<max) {
          void <- ggplot2::ggplot() +
            ggplot2::theme_void()
          pl.list[[p]][[paste0("void", length(pl.list[[p]])+1)]] <- void
        }

        join.plot <- ggpubr::ggarrange(plotlist = pl.list[[p]],
                                       ncol= 1,
                                       nrow = length(pl.list[[p]]))
        join.plot <- ggpubr::annotate_figure(join.plot,
                                             top = ggpubr::text_grob(p, face = "bold",
                                                                     color = "darkblue",
                                                                     size = 26))
        pl.list[[p]] <- join.plot
      }


      plots[[ob]] <- ggpubr::ggarrange(plotlist = pl.list,
                                       nrow= 1,
                                       ncol = length(pl.list))
    }
  })
  return(plots)
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
                      ncores = parallelly::availableCores() - 2,
                      progressbar = TRUE) {

  BiocParallel::bplapply(
    X = names(obj.list),
    BPPARAM =  BiocParallel::MulticoreParam(workers = ncores,
                                            progressbar = progressbar),
    function(x) {
      file_name <- file.path(dir, sprintf("%s.rds", x))
      saveRDS(obj.list[[x]], file_name)
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
                      ncores = parallelly::availableCores() - 2,
                      progressbar = TRUE) {

  if (!is.null(dir) & is.null(file.list)) {
    file_names <- list.files(dir)
    file_paths <- file.path(dir, file_names)
  } else if (is.null(dir) & !is.null(file.list)) {
    file_paths <- file.list[endsWith(files, '.rds')]
    file_names <- gsub("^.*/", "", file_paths)
  }
  obj.list <- BiocParallel::bplapply(
    X = file_paths,
    BPPARAM =  BiocParallel::MulticoreParam(workers = ncores,
                                            progressbar = progressbar),
    function(x) {
      readRDS(file.path(x))
    })
  names(obj.list) <- stringr::str_remove_all(file_names, '.rds')
  return(obj.list)
}



#' Add metadata to an annotated Seurat object to then obtain a HiT object
#'
#' Add metadata in misc slot of a Seurat object previously annotated to run \link{get.HiTObject}.
#' This needed data is by default incorporated by \link{Run.HiTME}, however, in objects annotated otherwise or splitted or merged these slots are lost. This metadata is needed basically to be able to link layer1 and layer2.
#'
#'
#' @param object Annotated Seurat object
#' @param scGate.models names of the scGate models used for layer1 annotation. Usually this is a list, here names of the list should be provided.
#' @param add.signatures names of the signatures computed during HiTME annotation, provide names of the list of signatures.
#' @param layer2_annotation Metadata column with layer2 annotation data.
#' @param ref.maps Reference maps used for annotation on layer2 can be provided to extract their information.
#' @param layer2_levels If reference maps are not provided, a vector or list with the levels of the layer2 annotation can be provided. These must be named with corresponding layer1_link.


#' @return Reference compatible with HiTME annotation, linking layer1 and layer2.
#' @export add.HiTMetadata
#'
#' @examples
#' # scGate models to run
#' HiT_scGate_models <- all_models$human$HiTME
#' # additional signatures
#' additional.signatures <- GetSignature(SignatuR[[species]][["Programs"]])
#'
#' # add scGate_link to ref.maps, by default coarse cell type cell ontology ID
#'  layer1.links <- list(CD8 = "CL:0000625", CD4 = "CL:0000624", DC = "CL:0000451", MoMac = "CL:0000576_CL:0000235")
#'  for (a in names(ref.maps)) {
#'   ref.maps[[a]]@misc$layer1_link <- layer1.links[[a]]
#'  }
#'
#'  # Run HiTME
#'  object <- Run.HiTME(object,
#'   scGate.model = HiT_scGate_models,
#'   ref.maps = ref.maps)
#'
#'  # misc data should be contained in the resulting object, if lost run:
#'  object <- add.HiTMetadata(object = object,
#'               scGate.models = names(HiT_scGate_models),
#'               add.signatures = names(additional.signatures),
#'               ref.maps = ref.maps)
#'
#' # Alternatively to specificying reference maps used, a list named with layer1_link and levels for each refernece map can be used:
#' layer2_levels <- list("CL:0000625" = c("CD8.NaiveLike, CD8.Tex, CD8.Tpex...),
#'                         ...)

add.HiTMetadata <- function (object = NULL,
                             scGate.models = NULL,
                             add.signatures = NULL,
                             layer2_annotation = "functional.cluster",
                             ref.maps = NULL,
                             layer2_levels = NULL) {

  if(is.null(object)){
    stop("Please provide an annotated Seurat object")
  }

  # add misc slot, removed when merging or splitting
  object@misc[["layer1_param"]] <- list()
  object@misc[["layer1_param"]][["scGate_models"]] <- scGate.models
  object@misc[["layer1_param"]][["additional.signatures"]] <- add.signatures

  if (!is.null(scGate.models)) {
    object$scGate_multi <- factor(object$scGate_multi,
                                  levels = scGate.models)
  }

  if(is.null(layer2_levels) && is.null(ref.maps)){
    warning("Not adding misc metadata for layer2 as no reference maps or layer2_levels were indicated.")
  }

  if (!layer2_annotation %in% names(object@meta.data)) {
    stop(layer2_annotation, "not found in object metadata")
  }

  if(!is.null(ref.maps)){
    ref.maps <- as.list(ref.maps)

    layer2_levels <- lapply(ref.maps, function(x) {
      unique(x[[layer2_annotation]])
    })

    names(layer2_levels) <- lapply(ref.maps, function(y) {
      y@misc$layer1_link
    })
  }


  object@meta.data[[layer2_annotation]] <- factor(object@meta.data[[layer2_annotation]],
                                                  levels = unlist(layer2_levels))


  object@misc[["layer2_param"]][[layer2_annotation]][["levels2_per_levels1"]] <- layer2_levels

  return(object)
}



