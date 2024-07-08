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
