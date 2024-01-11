

#' Classify cells using scGate and ProjecTILs.
#'
#' @param object A seurat object or a list of seurat objects
#' @param scGate.model The scGate model to use (use get_scGateDB() to get a list of available models), by default fetch the HiTME models: \code{ scGate::get_scGateDB(branch = scGate.model.branch)[[species]][["HiTME"]]}.
#' @param scGate.model.branch From which branch Run.HiTME fetch the scGate models, by default models are retrieved from \code{master} branch.
#' @param additional.signatures UCell additional signatures to compute on each cell, by default \code{SignatuR} programs are included.
#' @param ref.maps A named list of the ProjecTILs reference maps to use. They ought to be Seurat objects. It is recommended to add in reference object slost \code{misc} the identifier connecting to layer 1 classification (scGate): \code{ref.map@misc$layer1_link}
#' @param split.by A Seurat object metadata column to split by (e.g. sample names).
#' @param layer1_link Column of metadata linking layer1 prediction (e.g. scGate ~ CellOntology_ID) in order to perform subsetting for second layer classification.
#' @param return.Seurat Whether to return a Hit object or add the data into Seurat object metadata
#' @param group.by If return.Seurat = F, variables to be used to summarize HiTME classification data in HiT object. See \link{get.HiTObject} for more information.
#' @param useNA Whether to include not annotated cells or not (labelled as "NA" in the group.by.composition) when summarizing into HiT object. Can be "no", "ifany", or "always". See \code{?table} and \link{get.celltype.composition} for more information. Default is "ifany".
#' @param remerge When setting split.by, if remerge = TRUE one object will be returned. If remerge = FALSE a list of objects will be returned.
#' @param ncores The number of cores to use, by default all available cores minus 2 are used.
#' @param bparam A \code{BiocParallel::bpparam()} object that tells Run.HiTME how to parallelize. If provided, it overrides the `ncores` parameter.
#' @param progressbar Whether to show a progressbar or not
#'
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom parallelly availableCores
#' @importFrom dplyr mutate filter %>%
#' @importFrom tibble column_to_rownames
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
#' scGate_models_DB <- scGate::get_scGateDB(branch = "master", verbose = T, force_update = TRUE)
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
                       progressbar = TRUE){

  if (is.null(object)) {
    stop("Please provide a Seurat object or a list of them")
  }

  if(is.list(object) && !is.null(split.by)){
    stop("split.by only supported for a single Seurat object, not a list.\n
         Merge list before running HiTME")
  }

  # split object into a list if indicated
  if (!is.null(split.by)) {
    if(is.list(object)){
      stop("Split.by argument not supported when providing a list of Seurat objects. Set split.by = NULL or merge list.")
    }
    if(!split.by %in% names(object@meta.data)){
      stop(paste("split.by argument: ", split.by, " is not a metadata column in this Seurat object"))
    }
    object <- Seurat::SplitObject(object, split.by = split.by)

  } else {
    # if not applying split.by not merge by default
    remerge <-  FALSE
  }

  if(!return.Seurat && is.null(group.by)){
    warning("If setting return Seurat as FALSE, HiT summarized object will be returned. Need to indicate group.by variable indicating cell type classification\n
         e.g. group.by = list(\"layer1\" = c(\"scGate_multi\"),\"layer2\" = c(\"functional.cluster\"))\n
            Returning Seurat object, not HiT object.")
    return.Seurat = TRUE
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
    remerge <-  FALSE
    object <- list(object)
  }


  # Warning to reduce number of cores if file is huge
  if(any(lapply(object, ncol))>=30000){
    warning("Huge Seurat object, consider reducing number of cores to avoid memory issues")
  }


  # set paralelization parameters
  param <- set_paralel_params(ncores = ncores,
                               bparam = bparam,
                               progressbar = progressbar)


  # Get default additional signatures
  # based on species get SignatuR signatures
  if(!is.null(additional.signatures)){
    # convert to list if not
    if(!is.list(additional.signatures)){
      additional.signatures <- list(additional.signatures)
    }

    if(length(additional.signatures) == 1 && tolower(additional.signatures) == "default"){
      if(species == "human"){
        sig.species <- "Hs"
      } else if(species == "mouse"){
        sig.species <- "Mm"
      }
      additional.signatures <- SignatuR::GetSignature(SignatuR::SignatuR[[sig.species]][["Programs"]])
      message(" - Adding additional signatures: ", paste(names(additional.signatures), collapse = ", "), "\n")
    }
    } else {
      for(v in seq_along(additional.signatures)){
        if(is.null(names(additional.signatures)[[v]]) || is.na(names(additional.signatures)[[v]])){
          names(additional.signatures)[[v]] <- paste0("Additional_signature", v)
        }
      }
    }


  # Run scGate, if model is provided
  if(!is.null(scGate.model)) {
    message("## Running scGate\n")

    if(!is.list(scGate.model)){
      scGate.model <- list(scGate.model)
    }

    # Retrieve default scGate models if default
    if(length(scGate.model) == 1 && tolower(scGate.model) == "default"){
      scGate.model.branch <- scGate.model.branch[1]
      if(species == "human"){
        scGate.model <- scGate::get_scGateDB(branch = scGate.model.branch,
                                             force_update = F)[[species]][["HiTME"]]
      } else if(species == "mouse"){
        scGate.model <- scGate::get_scGateDB(branch = scGate.model.branch,
                                             force_update = F)[[species]][["HiTME"]]
      }
      message(" - Running scGate model for ", paste(names(scGate.model), collapse = ", "), "\n")
    }



  ## Run scGate
    object <- lapply(
              X = object,
              FUN = function(x) {
                    if(class(x) != "Seurat"){
                      stop("Not Seurat object included, cannot be processed.\n")
                    }
                    x <- scGate::scGate(x,
                                        model=scGate.model,
                                        additional.signatures = additional.signatures,
                                        BPPARAM = param,
                                        multi.asNA = T)
                    return(x)
              }
    )

    message("Finished scGate\n####################################################\n")
  } else {
    message("Not running coarse cell type classification as no scGate model was indicated.\n")
  }

    # Instance if we want to run additional signatures but not scGate
  if(is.null(scGate.model) && !is.null(additional.signatures)){

  message("Running additional Signatures but not Running scGate classification\n")
  object <- lapply(
    X = object,
    FUN = function(x) {
      if(class(x) != "Seurat"){
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
    if(!is.list(ref.maps)){
      ref.maps <- list(ref.maps)
    }

    # give name to ref.maps list if not present
    for(v in seq_along(ref.maps)){
      if(is.null(names(ref.maps)[[v]]) || is.na(names(ref.maps)[[v]])){
        names(ref.maps)[[v]] <- paste0("Map_", v)
      }
    }

    # check that all ref maps are Seurat objects
    if(suppressWarnings(!all(lapply(ref.maps, function(x){class(x) == "Seurat"})))){
      message("Some or all reference maps are not a Seurat object, please provide reference maps as Seurat objects.\nNot running Projectils.")
    } else {

      message("## Running Projectils\n")

        object <- lapply(
                  X = object,
                  FUN = function(x) {
                    if(class(x) != "Seurat"){
                      stop("Not Seurat object included, cannot be processed.\n")
                    }
                    x <- ProjecTILs.classifier.multi(x,
                                                     ref.maps = ref.maps,
                                                     bparam = param,
                                                     layer1_link = layer1_link)
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
                   function(x){
                     # add misc slot, removed when merging
                     x@misc[["layer1_param"]] <- list()
                     x@misc[["layer1_param"]][["scGate_models"]] <- names(scGate.model)
                     x@misc[["layer1_param"]][["additional.signatures"]] <- names(additional.signatures)

                     if(!is.null(scGate.model)){
                       x$scGate_multi <- factor(x$scGate_multi,
                                                levels = names(scGate.model))
                     }
                     if(!is.null(ref.maps)){
                       # get ref.maps all cells types
                       if("functional.cluster" %in% names(x@meta.data)){
                         all.levels <- lapply(ref.maps, function(x){
                           unique(x$functional.cluster)
                         })
                         x$functional.cluster <- factor(x$functional.cluster,
                                                        levels = unlist(all.levels))

                         # add each level to misc
                         names(all.levels) <- lapply(ref.maps, function(x){
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

                   })


  # if group.by parameters are not present in metadata, return Seurat
  if(!return.Seurat){
    if(suppressWarnings(!all(lapply(object, function(x){any(names(x@meta.data) %in% group.by)})))){
      message("None of the 'group.by' variables found, returning Seurat object not HiT object.")
      return.Seurat <- TRUE
    } else {
      message("\nBuilding HiT object\n")

      hit <- BiocParallel::bplapply(
        X = object,
        BPPARAM = param,
        FUN = function(x) {
          x <- get.HiTObject(x,
                             group.by = group.by,
                             useNA = useNA
                             )
        }
      )

    }
  }


  # if list is of 1, return object not list

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



# Function get.HiTObject
#' Get a HiT object summarizing cell type classification and aggregated profiles.
#'
#'
#'
#' @param object A Seurat object
#' @param group.by List with one or multiple Seurat object metadata columns with cell type predictions to group by (e.g. layer 1 cell type classification)
#' @param name.additional.signatures Names of additional signatures as found in object metadata to take into account.
#' @param useNA logical whether to return aggregated profile for NA (undefined) cell types, default is FALSE.
#' @param clr_zero_impute_perc Parameter for internal \link{get.celltype.composition}.
#' @param layers_links Parameter for internal \link{get.celltype.composition}
#' @param gene.filter List of genes to subset for aggregated expression. Parameter for internal \link{get.aggregated.profile}.
#' @param assay Parameter for internal \link{get.aggregated.profile}.
#' @param layer1_link Column of metadata linking layer1 prediction (e.g. scGate ~ CellOntology_ID) in order to perform subsetting for second layer classification.
#'
#' @importFrom methods setClass new
#' @import SeuratObject
#'
#' @return HiT object summarizing cell type classification and aggregated profiles.
#' @export get.HiTObject
#'


get.HiTObject <- function(object,
                          group.by = list("layer1" = c("scGate_multi"),
                                          "layer2" = c("functional.cluster")
                                          ),
                          name.additional.signatures = NULL,
                          useNA = FALSE,
                          clr_zero_impute_perc = 1,
                          layers_links = c("scGate_multi" = "functional.cluster"),
                          gene.filter = NULL,
                          assay = "RNA",
                          layer1_link = "CellOntology_ID"
                          ){


  if (is.null(object)) {
    stop("Please provide a Seurat object")
  } else if(class(object) != "Seurat"){
    stop("Not Seurat object included, cannot be processed.\n")
  }


  if(is.null(group.by)){
    stop("Please provide at least one grouping variable")
  }

  if(!is.list(group.by)){
    group.by <- as.list(group.by)
  }

  if(!any(group.by %in% names(object@meta.data))){
    stop("Group.by variables ", paste(group.by, collapse = ", "), " not found in metadata")
  } else if (!all(group.by %in% names(object@meta.data))){
    group.by <- group.by[group.by %in% names(object@meta.data)]
    message("Only found ", paste(group.by, collapse = ", ") , " as grouping variable for HiT Object.")
  }

  for(v in seq_along(group.by)){
    if(is.null(names(group.by)[[v]]) || is.na(names(group.by)[[v]])){
      names(group.by)[[v]] <- paste0("layer", v)
    }
  }


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
  sig <- grep(paste(name.additional.signatures, collapse = "|"),
              names(object@meta.data), value = T)
  sig.df <- object@meta.data[, sig]
  pred.list[[layer.scgate]][["additional_signatures"]] <- sig.df
  }

  if(!is.null(scgate.models)){
  scgate <- grep(paste(paste0(scgate.models, "$"), collapse = "|"),
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
  message("\nComputing cell type composition...\n")

  comp.prop <- get.celltype.composition(object,
                                        group.by.composition = group.by,
                                        useNA = useNA,
                                        clr_zero_impute_perc = clr_zero_impute_perc,
                                        layer1_link = layer1_link,
                                        layers_links = layers_links)
  # Compute avg expression
  message("\nComputing aggregated profile...\n")

  avg.expr <- get.aggregated.profile(object,
                                     group.by.aggregated = group.by,
                                     gene.filter = gene.filter,
                                     assay = assay,
                                     useNA = useNA)

  aggr.signature <- get.aggregated.signature(object,
                                             group.by.aggregated = group.by,
                                             name.additional.signatures = name.additional.signatures,
                                             useNA = useNA)


  hit <- methods::new("HiT",
                       metadata = object@meta.data,
                       predictions = pred.list,
                       aggregated_profile = list("Pseudobulk" = avg.expr,
                                                 "Signatures" = aggr.signature),
                       composition = comp.prop
                      )


  return(hit)


}


#' Calculate cell type composition or frequencies
#'
#' @param object A seurat object
#' @param group.by.composition The Seurat object metadata column(s) containing celltype annotations (provide as character vector, containing the metadata column name(s))
#' @param split.by A Seurat object metadata column to split by (e.g. sample names)
#' @param min.cells Set a minimum threshold for number of cells to calculate relative abundance (e.g. less than 10 cells -> no relative abundance will be calculated)
#' @param useNA Whether to include not annotated cells or not (labelled as "NA" in the group.by.composition). Can be defined separately for each group.by.composition (provide single boolean or vector of booleans)
#' @param layers_link Named vector indicating the relation of multiple layer classification, default is \code{c("scGate_multi", "functional.cluster")}. Name of the vector element ought to be layer1 and element layer2. This parameter is used to compute relative compositional data of layer2 within layer1 classes.
#' @param clr_zero_impute_perc To calculate the clr-transformed relative abundance ("clr_freq"), zero values are not allowed and need to be imputed (e.g. by adding a pseudo cell count). Instead of adding a pseudo cell count of flat +1, here a pseudo cell count of +1% of the total cell count will be added to all cell types, to better take into consideration the relative abundance ratios (e.g. adding +1 cell to a total cell count of 10 cells would have a different, i.e. much larger effect, than adding +1 to 1000 cells).
#' @param layer1_link Column of metadata linking layer1 prediction (e.g. scGate ~ CellOntology_ID) in order to perform subsetting for second layer classification.

#' @importFrom Hotelling clr
#' @importFrom dplyr group_by summarize filter ungroup mutate select left_join n
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
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
                                     min.cells = 10,
                                     useNA = FALSE,
                                     layers_links = c("scGate_multi" = "functional.cluster"),
                                     clr_zero_impute_perc = 1,
                                     layer1_link = "CellOntology_ID") {

  if (is.null(object)) {
    stop("Please provide a Seurat object or metadata as dataframe")
  }

  # input can be a Seurat object or a dataframe containing its meta.data
  # convert object to metadata if seurat object is provided
  if(class(object) == "Seurat"){
    meta.data <- object@meta.data
    if(is.null(meta.data)){
      stop("No metadata found in this Seurat object")
      }
  } else {
    stop("Not Seurat object or dataframe included, cannot be processed.\n")
  }

  if(is.null(group.by.composition)){
    stop("Please specificy a group.by.composition variable")
  }

  # Assess wheter split.by variable is in metadata
  if(!is.null(split.by) && !split.by %in% names(meta.data)){
    stop("Split.by variable not found in meta.data!\n")
  }

  if(length(useNA) != 1 && length(useNA) != length(group.by.composition)){
    stop("useNA variable must be of length 1 or the same length as group.by.composition (group.by)")
  }

  # convert group.by.composition to list if not
  if(!is.list(group.by.composition)){
    group.by.composition <- as.list(group.by.composition)
  }


  # Rename group.by.composition if not indicated
  for(v in seq_along(group.by.composition)){
    if(is.null(names(group.by.composition)[[v]]) || is.na(names(group.by.composition)[[v]])){
      names(group.by.composition)[[v]] <- group.by.composition[[v]]
    }

  }



  # Keep only grouping variables in metadata
  if(!any(group.by.composition %in% names(meta.data))){
    stop("Group.by variables ", paste(group.by.composition, collapse = ", "), " not found in metadata")
  } else if (!all(group.by.composition %in% names(meta.data))){
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
  for(l in seq_along(layers_links)){
    if(!all(c(names(layers_links[l]), layers_links[l]) %in%
           group.by.composition)){
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
      if (sum(ctable$cell_counts) < min.cells) {
        warning(paste("There are less than ", min.cells,
                      "cells detected for ", group.by.composition[[i]],
                      ". This is too few to calculate a reasonable celltype composition.\nIf needed, set parameter min.cells = 0."))
        next
      }

      # get proportion relative to layer1 types
      if(group.by.composition[[i]] %in% layers_links){
        lay <- layers_links[which(layers_links == group.by.composition[[i]])]
        if(lay %in% names(object@misc$layer2_param)){
          meta.split <- split(meta.data, meta.data[[names(lay)]])
          # switch list names to cell ontology ID
          cellonto_dic <- lapply(meta.split, function(x){
                          nam <- unique(x[[names(lay)]])
                          val <- unique(x[[layer1_link]])
                          names(val) <- nam
                          return(val)
            }) %>% unname() %>%  unlist()

          levs <- object@misc$layer2_param[[lay]]$levels2_per_levels1
          names(levs) <- names(cellonto_dic[match(names(levs), cellonto_dic)])


          # Filter only layer1 cells with representation in layer2
          meta.split <- meta.split[names(levs)]

          # add factor to group by layer1
          for(f in names(meta.split)){
            meta.split[[f]]$functional.cluster <- factor(meta.split[[f]]$functional.cluster,
                                                         levels = levs[[f]])
          }

          ctable.split <- lapply(meta.split, compositional_data,
                                 split.by = split.by,
                                 group.by.1 = group.by.composition[[i]],
                                 useNA = useNA[i])
          ctable <- c(all = list(ctable),
                         ctable.split)
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
#' @param assay Assay to retrieve information. By default "RNA".
#' @param useNA logical whether to return aggregated profile for NA (undefined) cell types, default is FALSE.
#' @param ... Extra parameters for internal Seurat functions: AverageExpression, AggregateExpression, FindVariableFeatures

#' @importFrom Seurat AverageExpression AggregateExpression FindVariableFeatures
#' @return Average and aggregated expression as a list of matrices for all genes and indicated gene lists filtering.
#' @export get.aggregated.profile


get.aggregated.profile <- function(object,
                                   group.by.aggregated = NULL,
                                   assay = "RNA",
                                   useNA = FALSE,
                                   ...) {

  if (is.null(object)) {
    stop("Please provide a Seurat object")
    if(class(object) != "Seurat"){
      stop("Please provide a Seurat object")
    }
  }

  if(is.null(group.by.aggregated)){
    stop("Please specificy a group.by.aggregated variable")
  }

  # convert group.by.aggregated to list if not
  if(!is.list(group.by.aggregated)){
    group.by.aggregated <- as.list(group.by.aggregated)
  }

  # Rename group.by.aggregated if not indicated
  # Rename group.by.aggregated if not indicated
  for(v in seq_along(group.by.aggregated)){
    if(is.null(names(group.by.aggregated)[[v]]) || is.na(names(group.by.aggregated)[[v]])){
      names(group.by.aggregated)[[v]] <- group.by.aggregated[[v]]
    }
  }

  if(!any(group.by.aggregated %in% names(object@meta.data))){
    stop("Group.by variables ", paste(group.by.aggregated, collapse = ", "), " not found in metadata")
  } else if (!all(group.by.aggregated %in% names(object@meta.data))){
    group.by.aggregated <- group.by.aggregated[group.by.aggregated %in% names(object@meta.data)]
    message("Only found ", paste(group.by.aggregated, collapse = ", ") , " as grouping variables.")
  }


  avg.exp <- list()

  # loop over different grouping

  for(i in names(group.by.aggregated)){

    # if useNA = TRUE, transform NA to character
    if(useNA){
      object@meta.data[is.na(object@meta.data[[group.by.aggregated[[i]]]]),
                       group.by.aggregated[[i]]] <- "NA"
    }

  # compute pseudobulk
    suppressWarnings(
      {
        avg.exp[[i]]  <-
          Seurat::AggregateExpression(object,
                                      group.by = group.by.aggregated[[i]],
                                      assays = assay,
                                      verbose = F,
                                      ...)[[assay]]
      })


    if(ncol(avg.exp[[i]]) == 1){
      for(av in names(avg.exp)){
        colnames(avg.exp[[i]]) <-
          unique(object@meta.data[!is.na(object@meta.data[[group.by.aggregated[[i]]]]),
                                  group.by.aggregated[[i]]])
      }
    }

      # add colnames if only one cell type is found
    if(ncol(avg.exp[[i]] ) == 1){
        colnames(avg.exp[[i]] ) <-
          unique(object@meta.data[!is.na(object@meta.data[[group.by.aggregated[[i]]]]),
                                  group.by.aggregated[[i]]])
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
#' @param name.additional.signatures Names of additional signatures to compute the aggregation per cell type.
#' @param fun Function to aggregate the signature, e.g. mean or sum.
#' @param useNA logical whether to return aggregated signatures for NA (undefined) cell types, default is FALSE.

#' @importFrom dplyr group_by summarize_at filter
#' @return Aggregated signature score for each indicated cell type grouping Results is NULL of not additional signatures are indicated or present in metadata.
#' @export get.aggregated.signature


get.aggregated.signature <- function(object,
                                     group.by.aggregated = NULL,
                                     name.additional.signatures = NULL,
                                     fun = mean,
                                     useNA = FALSE){


  # input can be a Seurat object or a dataframe containing its meta.data
  # convert object to metadata if seurat object is provided
  if(class(object) == "Seurat"){
    meta.data <- object@meta.data
    if(is.null(meta.data)){
      stop("No metadata found in this Seurat object")
    }
  } else if (class(object) == "data.frame"){
    meta.data <- object
  } else {
    stop("Not Seurat object or dataframe included, cannot be processed.\n")
  }


  if(is.null(group.by.aggregated)){
    stop("Please specificy a group.by.aggregated variable")
  }

  # convert group.by.aggregated to list if not
  if(!is.list(group.by.aggregated)){
    group.by.aggregated <- as.list(group.by.aggregated)
  }

  # Rename group.by.aggregated if not indicated
  for(v in seq_along(group.by.aggregated)){
    if(is.null(names(group.by.aggregated)[[v]]) || is.na(names(group.by.aggregated)[[v]])){
      names(group.by.aggregated)[[v]] <- group.by.aggregated[[v]]
    }
  }

  if(!any(group.by.aggregated %in% names(meta.data))){
    stop("Group.by variables ", paste(group.by.aggregated, collapse = ", "), " not found in metadata")
  } else if (!all(group.by.aggregated %in% names(meta.data))){
    group.by.aggregated <- group.by.aggregated[group.by.aggregated %in% names(meta.data)]
    message("Only found ", paste(group.by.aggregated, collapse = ", ") , " as grouping variable.")
  }

  if(is.null(name.additional.signatures)){
    name.additional.signatures <- object@misc$layer1_param$additional.signatures
  }

  if(is.null(name.additional.signatures)){
    message("No additional signatures indicated. Returning NULL")
    aggr.sig <- NULL
  } else {

    if(!any(grepl(paste(name.additional.signatures, collapse = "|"),
                  names(meta.data)))){
      stop("No additional signatures found in this object metadata")
    }

    add.sig.cols <- grep(paste(name.additional.signatures, collapse = "|"),
                         names(meta.data), value = T)

    aggr.sig <- list()

    for(e in names(group.by.aggregated)){
      aggr.sig[[e]] <- meta.data %>%
        dplyr::group_by(.data[[group.by.aggregated[[e]]]]) %>%
        dplyr::summarize_at(add.sig.cols, fun, na.rm = T)

        # filter out NA if useNA=F
      if(!useNA){
        aggr.sig[[e]] <- aggr.sig[[e]] %>%
                          dplyr::filter(!is.na(.data[[group.by.aggregated[[e]]]]))
      }

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
                       host = "https://dec2021.archive.ensembl.org/"
                      ){



  if(is.null(GO_accession)){
    stop("Please provide at least one GO accession ID (e.g. GO:0003700")
  } else {
    if(length(names(GO_accession)) == 0){
      names(GO_accession) <- GO_accession
    }
  }

  species <- tolower(species)


  # adapt species
  if(is.null(species)){
    stop("Please provide human or mouse as species")
  }
  species <- tolower(species)
  if(grepl("homo|sapi|huma", species)){
    dataset <- "hsapiens_gene_ensembl"
  } else if (grepl("mice|mus", species)){
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
      error = function(e){
        message("GO accession from biomaRt not possible")
        message(e)
        NULL
      }
    )

  if(!is.null(gene.data.bm)){
    gene.data <- gene.data.bm %>%
      dplyr::filter(go_id %in% GO_accession) %>%
      dplyr::filter(hgnc_symbol != "") %>%
      split(., as.factor(.$go_id)) %>%
      lapply(., function(x) x[,1])

    if(length(names(gene.data)) != length(GO_accession)){
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


#' Compute aggregation metrics of summarized data
#'
#'
#' @param object List of Hit class object
#' @param group.by Grouping variable for
#' @param min.cells Minimum number of cells to consider for aggregated expression.
#' @param metadata.vars Variables to show as metadata
#' @param score Score to compute clustering of samples based on cell type prediction, either Silhoutte or modularity.
#' @param dist.method Method to compute distance between celltypes, default euclidean.
#' @param ndim Number of dimensions to be use for PCA clustering metrics.
#' @param nVarGenes Number of variable genes to assess samples.
#' @param gene.filter Additional named list of genes to subset their aggregated expression, if "default" is indicated, all genes,
#' ribosomal genes, transcription factor (GO:0003700), cytokines (GO:0005125),
#'   cytokine receptors (GO:0004896), chemokines (GO:0008009), and chemokines receptors (GO:0004950) are subsetted and accounted for.
#' @param GO_accession Additional GO accessions to subset genes for aggregation, by default the list indicated in \code{gene.filter} are returned.
#' @param black.list List of genes to discard from clustering, if "default" object "default_black_list" object is used. Alternative black listed genes can be provided as a vector or list.
#' @param ncores The number of cores to use, by default, all available cores - 2.
#' @param bparam A \code{BiocParallel::bpparam()} object that tells how to parallelize. If provided, it overrides the `ncores` parameter.
#' @param progressbar Whether to show a progressbar or not

#' @importFrom dplyr mutate mutate_if filter %>% coalesce mutate_all full_join row_number
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom parallelly availableCores
#' @importFrom tibble rownames_to_column column_to_rownames remove_rownames
#' @importFrom caret nearZeroVar
#' @importFrom ggplot2 aes geom_point guides theme geom_col labs geom_hline guide_legend geom_vline theme_bw ggtitle
#' @importFrom DESeq2 DESeqDataSetFromMatrix vst estimateSizeFactors
#' @importFrom MatrixGenerics rowVars
#' @importFrom SummarizedExperiment assay
#' @importFrom BiocGenerics counts
#' @importFrom data.table rbindlist
#' @importFrom ggdendro ggdendrogram

#' @return Summarize of pseudobulk metrics for each celltype based on different predictions across all samples provided.
#' @export get.cluster.samples
#'

get.cluster.samples <- function(object = NULL,
                            group.by = list("layer1" = c("scGate_multi"),
                                            "layer2" = c("functional.cluster")
                            ),
                            min.cells = 10,
                            metadata.vars = NULL,
                            score = c("silhouette"),
                            dist.method = "euclidean",
                            ndim = 10,
                            nVarGenes = 500,
                            gene.filter = NULL,
                            GO_accession = NULL,
                            ntests = 100,
                            black.list = NULL,
                            ncores = parallelly::availableCores() - 2,
                            bparam = NULL,
                            progressbar = TRUE
                            ){


  if (is.null(object) || !is.list(object)) {
    stop("Please provide a list of HiT object containing more than one sample")
    if(suppressWarnings(!all(lapply(object, function(x){class(x) == "HiT"})))){
      stop("Not all components of the list are HiT objects.")
    }
  }

  #give name to list of hit objects
  for(v in seq_along(object)){
    if(is.null(names(object)[[v]]) || is.na(names(object)[[v]])){
      names(object)[[v]] <- paste0("Sample", v)
    }
  }


  if(is.null(group.by)){
    stop("Please provide at least one grouping variable for cell type classification common in all elements of the list of Hit objects")
  } else {
    if(!is.list(group.by)){
      group.by <- as.list(group.by)
    }
  # give name to list of grouping.by variables
    for(v in seq_along(group.by)){
      if(is.null(names(group.by)[[v]]) || is.na(names(group.by)[[v]])){
        names(group.by)[[v]] <- paste0("layer", v)
      }
    }
  }

  if(suppressWarnings(!all(lapply(object, function(x){any(group.by %in% names(x@metadata))})))) {
    stop("Not all supplied HiT object contain ", paste(group.by, collapse = ", "),
         "group.by elements in their metadata")
  } else {
    message("\n#### Group by ####")
    present <- sapply(group.by,
                      function(char){
                        unlist(lapply(object,
                                   function(df) {char %in% colnames(df@metadata)}
                                   ))

                        }
                      )
    present.sum <- colSums(present)

    for(l in names(group.by)){
      message("** ", l, " present in ", present.sum[[l]], " / ", length(object), " HiT objects.")
    }

  }

  if(!is.null(metadata.vars)){
    message("\n#### Metadata ####")
    if(suppressWarnings(!all(lapply(object, function(x){any(metadata.vars %in% names(x@metadata))})))) {
      message("Not all supplied HiT object contain ", paste(metadata.vars, collapse = ", "),
           "metadata elements in their metadata")
    }
      in.md <- sapply(metadata.vars,
                        function(char){
                          unlist(lapply(object,
                                        function(df) {char %in% colnames(df@metadata)}
                          ))

                        }
      )
      in.md.sum <- colSums(in.md)

      for(l in metadata.vars){
        message("** ", l, " present in ", in.md.sum[[l]], " / ", length(object), " HiT objects.")
      }

    }


# set paralelization parameters
  param <- set_paralel_params(ncores = ncores,
                              bparam = bparam,
                              progressbar = progressbar)

# prepare metadata
  md.all <- sapply(metadata.vars, function(x){
    BiocParallel::bplapply(X = names(object),
                           BPPARAM = param,
                           FUN =function(y){
                             unique(object[[y]]@metadata[[x]])
                           }
    ) %>% unlist()
  }) %>% as.data.frame() %>%
    dplyr::mutate(sample = names(object))
# Join data from all HiT objects

  # list to store joined data for each grouping variable
  join.list <- list()
  join.list[["aggregated"]] <- list()
  join.list[["composition"]] <- list()

  for(gb in names(group.by)){

    layer.present <- rownames(present)[present[,gb]]
    message("Computing metrics for composition of " , gb, "...")
      count  <-
        BiocParallel::bplapply(
          X = layer.present,
          BPPARAM = param,
          FUN = function(x){
              mat <- object[[x]]@composition[[gb]]$freq %>%
                      dplyr::mutate(sample = x)
              return(mat)
        }) %>%
        data.table::rbindlist(., use.names = T, fill = T) %>%
        # convert NA to 0
        dplyr::mutate_if(is.numeric, ~ifelse(is.na(.), 0, .)) %>%
        tibble::column_to_rownames("sample") %>%
        t() %>% as.matrix()

      join.list[["composition"]][[gb]] <- list("Composition" = count,
                                               "Metadata" = md.all)

      message("Computing metrics for aggregated profile of " , gb, "...")

      # list for each gene subset


    gf <- BiocParallel::bplapply(
      X = layer.present,
      BPPARAM = param,
      FUN = function(x){
        # remove cells with less than min.cells
        tab <- object[[x]]@composition[[gb]]$cell_counts
        names(tab) <- gsub("-", "_", names(tab))
        keep <- names(tab)[tab>=min.cells]

        if(length(keep) > 0){
        dat <- object[[x]]@aggregated_profile$Pseudobulk[[gb]]
        # in case _ and - are not match
        colnames(dat) <- gsub("-", "_", colnames(dat))

        #keep only cell type with more than min.cells
        dat <- dat[, keep, drop = F]
        celltype <- colnames(dat)[colnames(dat)!= "gene"]

        # accommodate colnames to merge then
        colnames(dat) <- paste(colnames(dat), x, sep = "_")
        dat <- dat %>%
                as.data.frame() %>%
                tibble::rownames_to_column("gene")
        md <- data.frame(rn = colnames(dat)[colnames(dat)!= "gene"],
                         sample = rep(x, ncol(dat)-1),
                         celltype = celltype)
        return(list("data" = dat,
                    "metadata" = md))
        } else {
          return(NULL)
        }
      })
      # remove NULL if generated
      gf <- gf[!sapply(gf, is.null)]

      # join all gene expression matrices using full_join
      data.all <- lapply(gf, function(x) x[["data"]]) %>%
                    reduce(full_join, by = "gene") %>%
                    # convert NA to 0
                    mutate_if(is.numeric, ~ifelse(is.na(.), 0, .)) %>%
                    tibble::column_to_rownames("gene") %>%
                    as.matrix()

      ## Summarize metadata
      md.agg <- lapply(gf, function(x) x[["metadata"]]) %>%
            data.table::rbindlist() %>%
            left_join(., md.all, by = "sample") %>%
            as.data.frame() %>%
            tibble::column_to_rownames("rn")

      join.list[["aggregated"]][[gb]] <- list("data" = data.all,
                                              "metadata" = md.agg)



      message("\nComputing clustering metrics of aggregated for ", gb)

      cc <- get.cluster.score(matrix = join.list[["aggregated"]][[gb]]$data,
                              metadata = join.list[["aggregated"]][[gb]]$metadata,
                              cluster.by = c("celltype", "sample"),
                              score = score,
                              ndim = ndim,
                              ntests = 100,
                              gene.filter = gene.filter,
                              GO_accession = GO_accession,
                              black.list = black.list,
                              nVarGenes = nVarGenes,
                              dist.method = dist.method)

      join.list[["aggregated"]][[gb]][["clustering"]] <- cc


  }





  return(join.list)

}




#' Render plots summarizing celltype proportions and distribution in samples
#'
#'
#' @param object List of Hit class object
#' @param group.by Grouping variable for
#' @param split.by If desired, indicate how to split plots by metadata variable.
#' @param by.x Determine which variable show on x-axis, if celltype it shows boxplots of the number of each cell type for each sample. If sample is indicated a barplot is shown with the relative proportions of each celltype within each sample.


#' @importFrom dplyr mutate group_by distinct
#' @importFrom ggplot2 aes geom_point facet_wrap geom_col labs geom_boxplot theme_bw
#' @importFrom data.table rbindlist

#' @return Plots showing the compositional cell type data across the different samples.
#' @export plot.celltype.freq
#'

plot.celltype.freq <- function(object = NULL,
                                group.by = list("layer1" = c("scGate_multi"),
                                                "layer2" = c("functional.cluster")
                                ),
                                split.by = NULL,
                                by.x = "celltype"
){


  if (is.null(object)) {
    stop("Please provide a single one or a list of HiT object")
  }

  if(!is.list(object)){
    object <- list(object)
  }

  if(suppressWarnings(!all(lapply(object, function(x){class(x) == "HiT"})))){
    stop("Not all components of the list are HiT objects.")
  }

  #give name to list of hit objects
  for(v in seq_along(object)){
    if(is.null(names(object)[[v]]) || is.na(names(object)[[v]])){
      names(object)[[v]] <- paste0("Sample", v)
    }
  }


  if(is.null(group.by)){
    stop("Please provide at least one grouping variable for cell type classification common in all elements of the list of Hit objects")
  } else {
    if(!is.list(group.by)){
      group.by <- as.list(group.by)
    }
    # give name to list of grouping.by variables
    for(v in seq_along(group.by)){
      if(is.null(names(group.by)[[v]]) || is.na(names(group.by)[[v]])){
        names(group.by)[[v]] <- paste0("layer", v)
      }
    }
  }

  if(suppressWarnings(!all(lapply(object, function(x){any(group.by %in% names(x@metadata))})))) {
    stop("Not all supplied HiT object contain ", paste(group.by, collapse = ", "),
         " group.by elements in their metadata")
}
  # remove group by if not present in any sample
    present <- sapply(group.by,
                      function(char){
                        any(unlist(lapply(object,
                                      function(df) {char %in% colnames(df@metadata)}
                        )))

                      }
    )
    if(length(group.by[present]) != length(group.by)){
      warning("Celltype classification for ", paste(as.vector(group.by[!present]), collapse = ", "),
              " not present in any object")

      group.by <- group.by[present]
    }


  # check splitting by variables
  if(!is.null(split.by)){
    if(suppressWarnings(!all(lapply(object, function(x){any(split.by %in% names(x@metadata))})))) {
      warning("Not all supplied HiT objects contain ", paste(split.by, collapse = ", "),
              " metadata elements in their metadata")
    }
    present <- sapply(split.by,
                    function(char){
                      any(unlist(lapply(object,
                                    function(df) {char %in% colnames(df@metadata)}
                      )))

                    }
    )
    if(length(split.by[present]) != length(split.by)){
      warning("Metadata variable ", paste(as.vector(split.by[!present]), collapse = ", "),
              " not present in any object")
      split.by <- split.by[present]
    }
  }

  # build all table
  df <- lapply(names(object), function(y){
    sel <- names(object[[y]]@metadata) %in% c(group.by, split.by)
    a <- object[[y]]@metadata[,sel, drop = F] %>%
          dplyr::mutate(sample = y)
    return(a)
  }) %>%
    data.table::rbindlist(fill = T) %>% as.data.frame()

  # compute counts
  c.list <- list()
  for(gr.by in as.vector(group.by)){
    #keep only needed columns
    rm <- group.by[group.by != gr.by] %>% as.character()
    kp <- which(!names(df) %in% rm)

    # count celltypes
    c.list[[gr.by]] <- df[ , kp] %>%
      dplyr::group_by(.data[[gr.by]], sample) %>%
      dplyr::mutate(count = dplyr::n()) %>%
      dplyr::distinct() %>%
      dplyr::group_by(sample) %>%
      dplyr::mutate(Freq = count/sum(count))


  }

  # list to store plots
  pl.list <- list()
  if(by.x == "celltype"){
    # count number of celltypes
    for(gr.by in as.vector(group.by)){
      pos <- ggplot2::position_dodge(width = .9)
      pl.list[[gr.by]] <-
        c.list[[gr.by]] %>%
        ggplot2::ggplot(ggplot2::aes(.data[[gr.by]], Freq,
                                     fill = .data[[gr.by]])) +
        ggplot2::geom_boxplot(outlier.colour = NA,
                              position = pos,
                              show.legend = F,
                              width = 0.6) +
        ggplot2::geom_point(position = pos,
                            show.legend = F) +
        ggplot2::labs(title = paste0(gr.by), " classification") +
        ggplot2::theme_bw()


    }


  } else if(by.x == "sample"){
    for(gr.by in as.vector(group.by)){
      pl.list[[gr.by]] <-
        c.list[[gr.by]] %>%
        ggplot2::ggplot(ggplot2::aes(sample, Freq, fill = .data[[gr.by]])) +
        ggplot2::geom_col() +
        ggplot2::labs(title = paste0(gr.by), " classification") +
        ggplot2::theme_bw()


    }
  }

  # rotat x-axis label
  pl.list <- lapply(pl.list, function(x){
    x + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                       hjust = 1,
                                                       vjust = 1))
  })

  for(spl in split.by){
    pl.list <- lapply(pl.list, function(x){
      x + ggplot2::facet_wrap(~.data[[spl]])
    })
  }



  return(pl.list)

}



#' Built confusion matrix-like plots between different cell type classification approaches
#'
#'
#' @param object List of Hit class object
#' @param var.1 1st of grouping variables on the x-axis. If using \code{relative = TRUE} proportions will be normalized to this variable.
#' @param var.2 2nd of grouping variables on the x-axis.
#' @param relative Wheter to show absolute number of cells, or relative number cells out of the 1st variable indicted. Default is FALSE.
#' @param useNA Whether to include not annotated cells or not (labelled as "NA"). Can be "no", "ifany", or "always". See \code{?table} for more information. Default is "ifany".
#' @param type Type of plot to render. Either "tile" (default) as a confusion matrixs or "barplot".

#' @importFrom dplyr mutate group_by distinct
#' @importFrom ggplot2 aes geom_point geom_tile geom_col scale_fill_gradient labs geom_label theme
#' @importFrom cowplot theme_cowplot
#' @importFrom data.table rbindlist

#' @return Plots to evaluate the correspondance between different classification methods.
#' @export plot.confusion.matrix
#'

plot.confusion.matrix <- function(object = NULL,
                                  var.1 = NULL,
                                  var.2 = NULL,
                                  relative = FALSE,
                                  useNA = "ifany",
                                  type = "tile"
){


  if (is.null(object)) {
    stop("Please provide a single one or a list of HiT object")
  }

  if(!is.list(object)){
    object <- list(object)
  }

  if(suppressWarnings(!all(lapply(object, function(x){class(x) == "HiT"})))){
    stop("Not all components of the list are HiT objects.")
  }

  #give name to list of hit objects
  for(v in seq_along(object)){
    if(is.null(names(object)[[v]]) || is.na(names(object)[[v]])){
      names(object)[[v]] <- paste0("Sample", v)
    }
  }

  # join vars to plot
  vars <- c(var.1, var.2)

  if(any(is.null(vars))){
    stop("Please provide 2 cell type classification labels common in all elements of the list of Hit objects: var.1 and var.2")
  }

  if(suppressWarnings(!all(lapply(object, function(x){any(vars %in% names(x@metadata))})))) {
    stop("Not all supplied HiT object contain ", paste(vars, collapse = ", "),
         " group.by elements in their metadata")
  }


  # build all table
  data <- lapply(names(object), function(y){
    sel <- names(object[[y]]@metadata) %in% vars
    a <- object[[y]]@metadata[,sel, drop = F] %>%
      dplyr::mutate(sample = y)
    return(a)
  }) %>%
    data.table::rbindlist(fill = T) %>% as.data.frame()

  # Get counts for every group
  pz <- table(data[[var.1]],
              data[[var.2]],
              useNA = useNA)
  if(relative){
    pz <- pz %>% prop.table(margin = 1)
    legend <- "Relative counts"
  } else {
    legend <- "Absolute counts"
  }


  if(tolower(type) == "tile"){
  plot <- pz %>% as.data.frame() %>%
    ggplot2::ggplot(ggplot2::aes(Var1, Var2,
               fill = Freq)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(name = legend,
                                  low = "#FFFFC8",
                                  high = "#7D0025") +
    ggplot2::labs(x = var.1,
                  y = var.2) +
    ggplot2::geom_label(ggplot2::aes(label = ifelse(Freq>0,
                                                    round(Freq,1),
                                                    NA)),
                        color = "white",
                        alpha = 0.6,
                        fill = "black") +
    cowplot::theme_cowplot()
  } else if(grepl("bar|col", type, ignore.case = T)){
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
                        split.by = NULL
){

  if (is.null(object)) {
    stop("Please provide a Seurat object or a list of them")
  }

  if(is.list(object) && !is.null(split.by)){
    stop("split.by only supported for a single Seurat object, not a list.\n
         Merge list before running HiTME")
  }

  if (is.null(scGate.model)) {
    stop("Please provide a scGate model or list of them")
  }

  if(!is.list(scGate.model)){
    scGate.model <- list("scGate_model" = scGate.model)
  }

  # split object into a list if indicated
  if (!is.null(split.by)) {
    if(is.list(object)){
      stop("Split.by argument not supported when providing a list of Seurat objects. Set split.by = NULL or merge list.")
    }
    if(!split.by %in% names(object@meta.data)){
      stop(paste("split.by argument: ", split.by, " is not a metadata column in this Seurat object"))
    }
    object <- Seurat::SplitObject(object, split.by = split.by)

  }

  if (!is.null(group.by)) {
    if(!group.by %in% names(object@meta.data)){
      stop(paste("group.by argument: ", group.by, " is not a metadata column in this Seurat object"))
    }
  } else {
    group.by <- "orig.ident"
  }

  # if object is unique turn into a list
  if(!is.list(object)){
    object <- list("object" = object)
  }

  plots <- list()
  # render plots
  model <- scGate:::table.to.model(scGate.model)

  suppressWarnings(
    {
  for(ob in names(object)){
    if(class(object[[ob]]) != "Seurat"){
      stop("Not Seurat object included, cannot be processed.\n")
    }

    # keep only genes expressed
    sc_names <- rownames(object[[ob]])[rowSums(object[[ob]]) > 0]
    # make plot list for each level
    pl.list <- list()

    for(a in names(model)){
      m <- model[[a]] %>% unlist(recursive = F)
      pl.sublist <- list()
      for(e in names(m)){
        # remove - sign on negative markers
        feat <- m[[e]] %>% gsub("-", "", .)

        feat <- intersect(feat, sc_names)

        if(length(feat)>0){
        # do not stack if only one gene is present
        stack <- ifelse(length(feat)>1, T, F)

        pl.sublist[[e]] <-
            Seurat::VlnPlot(object[[ob]],
                                     features = feat,
                                     group.by = group.by,
                                     stack = stack,
                                      pt.size = 0,
                                     flip = T) +
                            ggplot2::ggtitle(e) +
                            Seurat::NoLegend() +
                            ggplot2::xlab("") +
                            {if(length(feat) == 1){
                              ggplot2::ylab(feat)
                            }} +
                          ggplot2::theme(plot.title = element_text(hjust = 0.5))

        }
      }

      pl.list[[a]] <- pl.sublist
    }

    # max number of plots
    max <- lapply(pl.list, length) %>% unlist() %>% max()

    for(p in names(pl.list)){

      # add blank plots if needed
      while(length(pl.list[[p]])<max){
        void <- ggplot2::ggplot() + ggplot2::theme_void()
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
                      progressbar = T){
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
                      ncores = parallelly::availableCores() - 2,
                      progressbar = T){
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




