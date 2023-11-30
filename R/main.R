

#' Classify cells using scGate and ProjecTILs.
#'
#' @param object A seurat object or a list of seurat objects
#' @param scGate.model The scGate model to use (use get_scGateDB() to get a list of available models), by default fetch the HiTME models: \code{ scGate::get_scGateDB(branch = scGate.model.branch)[[species]][["HiTME"]]}.
#' @param scGate.model.branch From which branch Run.HiTME fetch the scGate models, by default models are retrieved from \code{master} branch.
#' @param additional.signatures UCell additional signatures to compute on each cell, by default \code{SignatuR} Programs are included.
#' @param ref.maps A named list of the ProjecTILs reference maps to use. They ought to be Seurat objects. It is recommended to add in reference object slost \code{misc} the identifier connecting to layer 1 classification (scGate): \code{ref.map@misc$layer1_link}
#' @param split.by A Seurat object metadata column to split by (e.g. sample names).
#' @param layer1_link Column of metadata linking layer1 prediction (e.g. scGate ~ scGate_multi) in order to perform subsetting for second layer classification.
#' @param return.Seurat Whether to return a Hit object or add the data into Seurat object metadata
#' @param group.by If return.Seurat = F, variables to be used to summarize HiTME classification data in HiT object.
#' @param useNA Whether to include not annotated cells or not (labelled as "NA" in the group.by.composition) when summarizing into HiT object.
#' @param remerge When setting split.by, if remerge = TRUE one object will be returned. If remerge = FALSE a list of objects will be returned.
#' @param ncores The number of cores to use
#' @param bparam A \code{BiocParallel::bpparam()} object that tells Run.HiTME how to parallelize. If provided, it overrides the `ncores` parameter.
#' @param progressbar Whether to show a progressbar or not
#'
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom parallelly availableCores
#' @importFrom dplyr mutate filter %>%
#' @importFrom tibble column_to_rownames
#' @importFrom scGate scGate get_scGateDB
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
                       progressbar = TRUE
                       ){

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
      stop(paste("split.by argument:", split.by, "is not a metadata column in this Seurat object"))
    }
    object <- Seurat::SplitObject(object, split.by = split.by)

  } else {
    # if not applying split.by not merge by default
    remerge <-  FALSE
  }

  if(!return.Seurat && is.null(group.by)){
    message("If setting return Seurat as FALSE, HiT summarized object will be returned. Need to indicate group.by variable indicating cell type classification\n
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

  if(is.null(ncores)){
    ncores <- 1
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


                    x@misc[["layer1_param"]] <- list()
                    x@misc[["layer1_param"]][["layer1_models"]] <- names(scGate.model)
                    x@misc[["layer1_param"]][["additional.signatures"]] <- names(additional.signatures)

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


      x@misc[["layer1_param"]] <- list()
      x@misc[["layer1_param"]][["additional.signatures"]] <- names(additional.signatures)

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
      message("Some or all reference maps are not a Seurat object, please prodive refenrece maps as Seurat objects.\nNot running Projectils.")
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


  if (is.list(object)) {
    if (remerge && length(object)>1) {
        object <- merge(object[[1]],object[2:length(object)])
        # add misc slot, removed when merging
        object@misc[["layer1_param"]] <- list()
        object@misc[["layer1_param"]][["layer1_models"]] <- names(scGate.model)
        object@misc[["layer1_param"]][["additional.signatures"]] <- names(additional.signatures)

        object <- list(object)
    }
  } else {
    object <- list(object)
  }


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
#' Get a HiT object summarizing cell type classification and aggregated profiles.
#'
#'
#'
#' @param object A Seurat object
#' @param group.by List with one or multiple Seurat object metadata columns with cell type predictions to group by (e.g. layer 1 cell type classification)
#' @param name.additional.signatures Names of additional signatures as found in object metadata to take into account.
#' @param clr_zero_impute_perc Parameter for internal \link{get.celltype.composition}.
#' @param gene.filter List of genes to subset for aggregated expression. Parameter for internal \link{get.aggregated.profile}.
#' @param assay Parameter for internal \link{get.aggregated.profile}.
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
                          gene.filter = NULL,
                          assay = "RNA"
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
    if(is.null(names(group.by[[v]])) || is.na(names(group.by[[v]]))){
      names(group.by[[v]]) <- paste0("layer", v)
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
  message("Computing cell type composition...\n")

  comp.prop <- get.celltype.composition(object,
                                        group.by.composition = group.by,
                                        useNA = useNA,
                                        clr_zero_impute_perc = clr_zero_impute_perc)
  # Compute avg expression
  message("Computing aggregated profile...\n")

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
#' @param object A seurat object or metadata dataframe.
#' @param group.by.composition The Seurat object metadata column(s) containing celltype annotations (provide as character vector, containing the metadata column name(s))
#' @param split.by A Seurat object metadata column to split by (e.g. sample names)
#' @param min.cells Set a minimum threshold for number of cells to calculate relative abundance (e.g. less than 10 cells -> no relative abundance will be calculated)

#' @param useNA Whether to include not annotated cells or not (labelled as "NA" in the group.by.composition). Can be defined separately for each group.by.composition (provide single boolean or vector of booleans)
#' @param clr_zero_impute_perc To calculate the clr-transformed relative abundance ("clr_freq"), zero values are not allowed and need to be imputed (e.g. by adding a pseudo cell count). Instead of adding a pseudo cell count of flat +1, here a pseudo cell count of +1% of the total cell count will be added to all cell types, to better take into consideration the relative abundance ratios (e.g. adding +1 cell to a total cell count of 10 cells would have a different, i.e. much larger effect, than adding +1 to 1000 cells).

#' @importFrom Hotelling clr
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
                                     clr_zero_impute_perc = 1) {

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
  } else if (class(object) == "data.frame"){
    meta.data <- object
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

  if(length(useNA) != 1 && length(useNA) == length(group.by.composition)){
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


  # Convert TRUE/FALSE to ifany or no
  useNA <- ifelse(useNA == TRUE, "ifany", "no")
  # Rep useNa parameters according to group.by.composition variable
  if (length(useNA) == 1) {
    useNA <- rep(useNA, length(group.by.composition))
  }

  # list to store compositions
  celltype.compositions <- list()


  for (i in seq_along(group.by.composition)) {
    if (is.null(split.by)) {
      comp_table <- table(meta.data[[group.by.composition[[i]]]],
                          useNA = useNA[i])
      comp_table_freq <- prop.table(comp_table) * 100 # To get percentage
    } else {
      comp_table <- table(meta.data[[group.by.composition[[i]]]],
                          meta.data[[split.by]],
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
    celltype.compositions[[names(group.by.composition)[i]]][["cell_counts"]] <- as.data.frame.matrix(t(comp_table))
    celltype.compositions[[names(group.by.composition)[i]]][["freq"]] <- as.data.frame.matrix(t(comp_table_freq))
    celltype.compositions[[names(group.by.composition)[i]]][["freq_clr"]] <- as.data.frame.matrix(comp_table_clr)
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
#' @param gene.filter Additional named list of genes to subset their aggregated expression, by default all genes,
#' ribosomal genes, transcription factor (GO:0003700), cytokines (GO:0005125),
#'   cytokine receptors (GO:0004896), chemokines (GO:0008009), and chemokines receptors (GO:0004950) are subsetted.
#' @param GO_accession Additional GO accessions to subset genes for aggregation, by default the list indicated in \code{gene.filter} are returned.
#' @param assay Assay to retrieve information. By default "RNA".
#' @param useNA logical whether to return aggregated profile for NA (undefined) cell types, default is FALSE.
#' @param ... Extra parameters for internal Seurat functions: AverageExpression, AggregateExpression, FindVariableFeatures

#' @importFrom Seurat AverageExpression AggregateExpression FindVariableFeatures
#' @return Average and aggregated expression as a list of matrices for all genes and indicated gene lists filtering.
#' @export get.aggregated.profile


get.aggregated.profile <- function(object,
                                   group.by.aggregated = NULL,
                                   gene.filter = NULL,
                                   GO_accession = NULL,
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

  if(!is.null(gene.filter) && !is.list(gene.filter)){
      gene.filter <- list(gene.filter)
  }


  gene.filter.list <- list()
  # make a list of default subsetting genes

  # Dorothea transcription factors
  # data("entire_database", package = "dorothea")
  # gene.filter.list[["Dorothea_Transcription_Factors"]] <- entire_database$tf %>% unique()

  # Ribosomal genes
  gene.filter.list[["Ribosomal"]] <- SignatuR::GetSignature(SignatuR::SignatuR$Hs$Compartments$Ribo)

  # Add list genes from GO accessions
  # check if indicted additional GO accessions
  if(!is.null(GO_accession)){
    GO_additional <- get.GOList(GO_accession, ...)
  } else {
    GO_additional <- NULL
  }

  #default gene list
  data("GO_accession_default")
  gene.filter.list <- c(gene.filter.list, GO_default, GO_additional)

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

  for(i in names(group.by.aggregated)){
    avg.exp[[i]] <- list()

    # if useNA = TRUE, transform NA to character
    if(useNA){
      object@meta.data[is.na(object@meta.data[[group.by.aggregated[[i]]]]),
                       group.by.aggregated[[i]]] <- "NA"
    }

  # compute pseudobulk
    suppressWarnings(
      {
        avg.exp[[i]][["All.genes"]] <-
          Seurat::AggregateExpression(object,
                                      group.by = group.by.aggregated[[i]],
                                      assays = assay,
                                      verbose = F,
                                      ...)[[assay]]
      })


    if(ncol(avg.exp[[i]][["All.genes"]]) == 1){
      for(av in names(avg.exp)){
        colnames(avg.exp[[av]][[i]][["All.genes"]]) <-
          unique(object@meta.data[!is.na(object@meta.data[[group.by.aggregated[[i]]]]),
                                  group.by.aggregated[[i]]])
      }
    }

      # add colnames if only one cell type is found
    if(ncol(avg.exp[[i]][["All.genes"]]) == 1){
        colnames(avg.exp[[i]][["All.genes"]]) <-
          unique(object@meta.data[!is.na(object@meta.data[[group.by.aggregated[[i]]]]),
                                  group.by.aggregated[[i]]])
      }

      # subset genes accroding to gene filter list
      for(e in names(gene.filter.list)){
        keep <- gene.filter.list[[e]][gene.filter.list[[e]] %in%
                                        rownames(avg.exp[[i]][["All.genes"]])]
        avg.exp[[i]][[e]] <- avg.exp[[i]][["All.genes"]][keep, , drop = FALSE]
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
#' @param ncores The number of cores to use
#' @param bparam A \code{BiocParallel::bpparam()} object that tells how to parallelize. If provided, it overrides the `ncores` parameter.
#' @param progressbar Whether to show a progressbar or not

#' @importFrom dplyr mutate filter %>% coalesce mutate_all full_join
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom parallelly availableCores
#' @importFrom tibble rownames_to_column
#' @importFrom caret nearZeroVar
#' @importFrom ggplot2 aes geom_point guides theme geom_col labs geom_hline guide_legend

#' @return List of genes for each GO accession requested. If named vector is provided, lists output are named as GO:accession.ID_named (e.g. "GO:0004950_cytokine_receptor_activity")
#' @export get.cluster.samples
#'

get.cluster.samples <- function(object,
                            group.by = list("layer1" = c("scGate_multi"),
                                            "layer2" = c("functional.cluster")
                            ),
                            min.cells = 10,
                            metadata.vars = NULL,
                            score = c("silhouette"),
                            dist.method = "euclidean",
                            ndim = 10,
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
    stop("Not all supplied HiT object contain", paste(group.by, collapse = ", "),
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
      message("Not all supplied HiT object contain", paste(metadata.vars, collapse = ", "),
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
  if(is.null(ncores)){
    ncores <- 1
  }

  if(ncores >= parallelly::availableCores()){
    ncores <- parallelly::availableCores() - 1
    message("Using all or more cores available in this computer, reducing number of cores to ", ncores)
  }

  # set paralelization parameters
  if (is.null(bparam)) {
    if (ncores>1) {
      param <- BiocParallel::MulticoreParam(workers =  ncores,
                                            progressbar = progressbar)
    } else {
      param <- SerialParam()
    }
  } else {
    param <- bparam
  }

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
        mutate_if(is.numeric, ~ifelse(is.na(.), 0, .)) %>%
        tibble::column_to_rownames("sample") %>%
        t() %>% as.matrix()

      join.list[["composition"]][[gb]] <- list("Composition" = count,
                                               "Metadata" = md.all %>%
                                                            tibble::column_to_rownames("sample"))

      message("Computing metrics for aggregated profile of " , gb, "...")

      # get gene subsets
      gene.filter <- unique(unlist(lapply(object, function(x){
                      names(x@aggregated_profile$Gene_expression$Average[[gb]])
                     })))
      # list for each gene subset

      join.list[["aggregated"]][[gb]] <-
        BiocParallel::bplapply(X = gene.filter,
                               BPPARAM = param,
                               FUN =function(y){
                gf <- lapply(
                  X = layer.present,
                  FUN = function(x){
                    # remove cells with less than min.cells
                    tab <- object[[x]]@composition[[gb]]$cell_counts
                    names(tab) <- gsub("-", "_", names(tab))
                    keep <- names(tab)[tab>=min.cells]
                    dat <- object[[x]]@aggregated_profile$Gene_expression$Aggregated[[gb]][[y]]
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
            })
            data.all <- lapply(gf, function(x) x[["data"]]) %>%
                    reduce(full_join, by = "gene") %>%
            # convert NA to 0
            mutate_if(is.numeric, ~ifelse(is.na(.), 0, .)) %>%
            tibble::column_to_rownames("gene") %>%
            as.matrix()

            ## Summarize metadata
            md.all <- lapply(gf, function(x) x[["metadata"]]) %>%
                  data.table::rbindlist() %>%
                  left_join(., md.all, by = "sample") %>%
                  as.data.frame() %>%
                  tibble::column_to_rownames("rn")

            # compute clustering metrics
            # return list
            return(list("data" = data.all,
                        "metadata" = md.all))
        })
      # give names to list
      names(join.list[["aggregated"]][[gb]]) <- gene.filter

      # there seems to be problems in running this in pararlel...
      message("\nComputing clustering metrics of aggregated for ", gb)
      for(c in names(join.list[["aggregated"]][[gb]])){
        message("Computing clustering metrics of aggregated for ", gb, " ", c)
        cc <- get.cluster.score(matrix = join.list[["aggregated"]][[gb]][[c]]$data,
                                metadata = join.list[["aggregated"]][[gb]][[c]]$metadata,
                                cluster.by = c("sample", "celltype"),
                                score = score,
                                dist.method = dist.method)
        join.list[["aggregated"]][[gb]][[c]][["clustering"]] <- cc
      }

      md <- join.list[["aggregated"]][[gb]][[1]]$metadata %>% data.frame()
      data <- lapply(join.list[["aggregated"]][[gb]], function(x) {x[-which(names(x) == "metadata") ]})
      join.list[["aggregated"]][[gb]] <- c(data,list("metadata" = md))


  }





  return(join.list)

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




