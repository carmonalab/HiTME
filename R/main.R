#' Annotate cell types at multiple levels of granularity in a scRNA-seq dataset.
#'
#'The function takes as input Seurat objects (or list of them).
#'These should be split by sample to avoid batch effects, or split internally in by indicating the parameter \code{split.by}.
#'<br>
#'This function firstly runs \link[scGate]{scGate} (easily customizable) marker-based classification, resulting in a coarse-grained cell type classification (CD4T, B cell, Dendritic cell...).
#'Next, it runs for each broad cell type \link[ProjecTILs]{ProjecTILs} for a finer cell type classification (CD4+ TFH, Tex CD8+, cDC1...) based on cell mapping onto expert-curated single-cell reference maps.
#'
#' @param object A seurat object or a list of seurat objects
#' @param scGate.model The \link[scGate]{scGate} model to use. Use \link[scGate]{get_scGateDB} to get a list of available models. By default it fetchs the HiTME models: \code{ scGate::get_scGateDB(branch = scGate.model.branch)[[species]][["HiTME"]]}
#' @param scGate.model.branch From which branch Run.HiTME fetch the scGate models, by default models are retrieved from \code{master} branch.
#' @param multi.asNA How to label cells that are "Pure" for multiple annotations: "Multi" (FALSE) or NA (TRUE)
#' @param additional.signatures Additional signature(s) to compute on each cell using \link[UCell]{UCell}. Multiple signatures should be provided as a list.
#' @param ref.maps A named list of the \link[ProjecTILs]{ProjecTILs} reference maps to use. Download default reference maps using \link[ProjecTILs](get.reference.maps). For custom reference maps, they must be Seurat objects turn into reference with \link[ProjecTILs](make.reference). It is recommended to add in reference object slost \code{misc} the identifier connecting to layer 1 classification (scGate): \code{ref.map@misc$layer1_link}
#' @param layer3 Gene signature(s) to define gene programs or cell states on top of layer2 (\link[ProjecTILs]{ProjecTILs} classification). By default \link[SignatuR]{SignatuR} programs for cell cycling, IFN response and HeatShock response are run.
#' @param layer3.threshold Minimum value threshold to apply on \link[UCell]{UCell} scores of layer3 signatures to classify cell types.
#' @param split.by A Seurat object metadata column to split by (e.g. sample names).
#' @param layer1_link Column of metadata linking layer1 prediction (e.g. scGate ~ CellOntology_ID) in order to perform subsetting for second layer classification.
#' @param remerge When setting split.by or providing a list of objects, if \code{remerge = TRUE} one object will be returned (default). If \code{remerge = FALSE} a list of objects will be returned.
#' @param species Define species to get the default \link[scGate]{scGate} models and \link[SignatuR]{SignatuR} signatures. Currently only human and mouse are supported, if other species are studied, set to \code{NULL}
#' @param ncores The number of cores to use, by default all available cores minus 2 are used.
#' @param bparam A \link[BiocParallel]{bpparam} object that tells Run.HiTME how to parallelize. If provided, it overrides the `ncores` parameter.
#' @param progressbar Whether to show a progress bar or not
#' @param verbose Verbose output, by default output messsage are returned
#'
#' @importFrom BiocParallel MulticoreParam bplapply
#' @importFrom parallel detectCores
#' @importFrom dplyr mutate filter rowwise all_of ungroup across c_across starts_with %>%
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom scGate scGate get_scGateDB
#' @import SignatuR
#' @import scGate
#' @import ProjecTILs
#' @importFrom Seurat SplitObject
#' @importFrom data.table rbindlist setDT
#'
#' @return Seurat object with additional metadata showing cell type classification.
#' @export Run.HiTME
#'
#' @examples
#' library(ProjecTILs)
#' library(HiTME)
#'
#' # get reference maps for human
#' ref.maps <- get.reference.maps(collection  = "human")
#'
#' # Run HiTME with default settings, specifying the human ref.maps
#' # Internally default scGate models and gene programs signatures will be fetched
#'
#'  query <- Run.HiTME(object = query,
#'                     ref.maps = ref.maps[["human"]]
#'                     )
#'

Run.HiTME <- function(object = NULL,
                      scGate.model = "default",
                      scGate.model.branch = c("master", "dev"),
                      multi.asNA = TRUE,
                      additional.signatures = NULL,
                      ref.maps = NULL,
                      layer3 = "default",
                      layer3.threshold = 0.2,
                      split.by = NULL,
                      layer1_link = "CellOntology_ID",
                      remerge = TRUE,
                      species = "human",
                      bparam = NULL,
                      ncores = parallel::detectCores() - 2,
                      progressbar = TRUE,
                      verbose = TRUE) {

  if (is.null(object)) {
    stop("Please provide a Seurat object or a list of them")
  }

  # check if objects are Seurat objects
  suppressWarnings(
    {
      if (!(inherits(object, "Seurat") ||
            all(lapply(object, inherits, "Seurat")))) {
        stop("All or some objects are not Seurat objects, cannot be processed\n")
      }
    })

  # split object into a list if indicated
  if (!is.null(split.by)) {
    if (is.list(object)) {
      stop("split.by only supported for a single Seurat object, not a list.\n
         Merge list before running HiTME or set split.by = NULL")
    }
    if (!split.by %in%
        names(object@meta.data)) {
      stop(paste("split.by argument: ",
                 split.by, " is not a metadata column in this Seurat object"))
    }
    object <- Seurat::SplitObject(object,
                                  split.by = split.by)
  }

  # adapt species
  if (is.null(species)) {
    warning("Not using default scGate models and signatures for supported species (human or mouse).\n
            Make sure to provide models for scGate.model and additional.signatures for cell type annotation.")
  } else {

    species <- tolower(species)
    if (grepl("homo|sap|huma", species)) {
      species <- "human"
      sig.species <- "Hs"
    }
    else if (grepl("mice|mus", species)) {
      species <- "mouse"
      sig.species <- "Mm"
    } else {
      stop("Only supported species for default scGate models and SignatuR signatures are human and mouse. \n
         If other species are studied, set species = NULL and provide models for scGate.model and additional.signatures")
    }
  }

  # if object is unique turn into a list
  if (!is.list(object)) {
    # not remerge because it is a single object
    remerge <-  FALSE
    # convert into a length=1 list to run
    object <- list(object)
  }


  # Warning to reduce number of cores if file is huge
  if (suppressWarnings(any(lapply(object, ncol)>=15000))) {
    warning("Huge Seurat object, consider reducing number of cores to avoid memory issues")
  }


  # set parallelization parameters
  param <- set_parallel_params(ncores = ncores,
                               bparam = bparam,
                               progressbar = progressbar)


  # Adapt additional signatures
  additional.signatures <- adapt_vector(additional.signatures,
                                        prefix = "Additional_signature_")


  # Set up layer3 params
  if (!is.null(layer3)) {
    if (length(layer3) == 1 &&
        tolower(layer3) == "default") {

      if (!is.null(species)) {
        layer3 <- SignatuR::GetSignature(SignatuR::SignatuR[[sig.species]][["Programs"]])
        selected.SignatuR.programs <- c("IFN",
                                        "HeatShock",
                                        "cellCycle.G1S",
                                        "cellCycle.G2M")
        layer3 <- layer3[selected.SignatuR.programs]
      } else {
        layer3 <- NULL
      }
    } else {
      layer3 <- adapt_vector(layer3,
                             prefix = "Sig3_")
    }
  }



  # Run scGate, if model is provided
  if (!is.null(scGate.model)) {

    # Retrieve default scGate models if default
    if (length(scGate.model) == 1 &&
        tolower(scGate.model) == "default" &&
        !is.null(species)) {

      if (verbose) {message("- Retrieving default HiTME scGate models for ", species, "\n")}

      scGate.model.branch <- scGate.model.branch[1]
      scGate.model <- scGate::get_scGateDB(branch = scGate.model.branch,
                                           verbose = FALSE)[[species]][["HiTME"]]
    }
    if (verbose) {message("### Running scGate\n")}

    ## Run scGate
    object <- lapply(
      X = object,
      function(x) {

        x <- scGate::scGate(x,
                            model=scGate.model,
                            additional.signatures = c(additional.signatures, layer3),
                            BPPARAM = param,
                            multi.asNA = multi.asNA,
                            verbose = verbose)
        return(x)
      }
    )

    message("Finished scGate\n#########################\n")
  } else {
    message("Not running coarse cell type classification as no scGate model was indicated.\n")
    # Add anyway a column with only NA
    object <- lapply(
      X = object,
      function(x) {
        x@meta.data[["scGate_multi"]] <- NA
        return(x)
      }
    )
  }

  # Instance if we want to run additional signatures but not scGate
  if (is.null(scGate.model) &&
      !is.null(c(additional.signatures, layer3))) {

    if (verbose) {message("### Running Gene Signatures but not Running scGate classification\n")}

    object <- lapply(
      X = object,
      function(x) {

        x <- UCell::AddModuleScore_UCell(x,
                                         features = c(additional.signatures, layer3),
                                         BPPARAM = param)
        return(x)
      }
    )
  }

  # adapt reference maps if needed
  ref.maps <- adapt_vector(ref.maps,
                           prefix = "ReferenceMap_")

  # Run ProjecTILs if ref.maps is provided
  if (!is.null(ref.maps)) {

    # check that all ref maps are Seurat objects
    if (suppressWarnings(!all(lapply(ref.maps, function(x) {inherits(x, "Seurat")})))) {
      warning("Some or all reference maps are not a Seurat object, please provide reference maps as Seurat objects.\nNot running Projectils.")
    } else {

      if (verbose) {message("### Running Projectils\n")}

      object <- lapply(
        X = object,
        function(x) {

          x <- ProjecTILs.classifier.multi(x,
                                           ref.maps = ref.maps,
                                           bparam = param,
                                           layer1_link = layer1_link,
                                           verbose = verbose)

          return(x)
        }
      )
      if (verbose) {message("Finished Projectils\n###########################\n")}
    }

  } else {
    if (verbose) {message("Not running reference mapping as no reference maps were indicated.\n")}
  }


  # add columns for functional.cluster if not present

  object <- lapply(
    X = object,
    function(x) {
      # Check if ProjecTIL columns were added
      # If not, add NA columns
      ProjecTILs_cols <- c("functional.cluster",
                           "functional.cluster.conf")
      if (!any(ProjecTILs_cols %in% names(x@meta.data))) {
        x@meta.data[ProjecTILs_cols[!ProjecTILs_cols %in% names(x@meta.data)]] <- NA
      }
      return(x)
    }
  )


  # Join layer1(scGate) and layer2(ProjecTILs)
  object <- lapply(
    X = object,
    function(x) {

      x@meta.data <- x@meta.data %>%
        mutate(layer2 = ifelse(is.na(functional.cluster),
                               scGate_multi, functional.cluster))

      return(x)
    }
  )

  # Add layer 3
  if (!is.null(layer3)) {



    if (verbose) {message("### Running Layer3 cell types\n")}

    # get signatures for indicated layer3
    sigs.cols <- paste0(names(layer3), "_UCell")

    object <- BiocParallel::bplapply(X = object,
                                     BPPARAM = param,
                                     function(x) {

                                       get.layer3(x,
                                                  ann.col = "layer2",
                                                  sigs.cols = sigs.cols,
                                                  layer3.threshold = layer3.threshold)

                                     })
  } else {
    if (verbose) {message("Not running Layer3 cell types as no signatures were indicated\n")}
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
                     x@misc[["additional.signatures"]] <- names(additional.signatures)
                     x@misc[["layer3_param"]] <- names(layer3)

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

  # if list is of 1, return object not list
  if (length(object)==1) {
    object <- object[[1]]
  }

  return(object)

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



#' Build confusion matrix-like or bar plots between different cell type classification approaches
#'
#'
#' @param object Seurat object or its metadata with cell type annotation
#' @param var.1 1st of grouping variables on the x-axis. If using \code{relative = TRUE} proportions will be normalized to this variable.
#' @param var.2 2nd of grouping variables on the x-axis.
#' @param relative Whether to show absolute number of cells, or relative number cells out of the 1st variable indicted. Default is FALSE.
#' @param useNA Whether to include not annotated cells or not (labelled as "NA"). Can be "no", "ifany", or "always". See \code{?table} for more information. Default is "ifany".
#' @param type Type of plot to render. Either "tile" (default) as a confusion matrices or "barplot".
#' @param xlab Label for x-axis
#' @param ylab Label for y-axis
#' @param plot.title Plot title

#' @importFrom ggplot2 aes geom_point geom_tile geom_col scale_fill_gradient labs geom_label theme theme_bw element_text

#' @return Plots to evaluate the correspondence between different classification methods.
#' @export plot.confusion
#'

plot.confusion <- function(object = NULL,
                           var.1 = NULL,
                           var.2 = NULL,
                           relative = FALSE,
                           useNA = "ifany",
                           type = "tile",
                           xlab = NULL,
                           ylab = NULL,
                           plot.title = "") {

  if (is.null(object)) {
    stop("Please provide a Seurat object or its metadata")
  }

  if (inherits(object, "Seurat")) {
    data <- object@meta.data
  } else if (is.data.frame(object)) {
    data <- object
  } else {
    stop("Please provide a Seurat object or its metadata")
  }

  legend.title <- "Number of cells"
  data <- as.data.frame(data)


  if (is.null(xlab)) {
    xlab <-  var.1
  }
  if (is.null(ylab)) {
    ylab <-  var.2
  }

  vars <- c(var.1, var.2)
  if (any(is.null(vars))) {
    stop("Please provide 2 cell type classification columns in metadata as var.1 and var.2")
  }

  if (any(!vars %in% names(data))) {
    stop("Not all classification variables: ",
         paste(vars, collapse = ", "),
         " were found in metadata")
  }

  # Get counts for every group
  pz <- table(data[[var.1]],
              data[[var.2]],
              useNA = useNA)

  if (relative) {
    pz <- prop.table(pz, margin = 1) %>%
      round(digits = 2)
    legend.title <- "Proportion of cells"
  }


  if (tolower(type) == "tile") {
    plot <- pz %>%
      as.data.frame() %>%
      ggplot2::ggplot(ggplot2::aes(Var1, Var2,
                                   fill = Freq)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient(name = legend.title,
                                   low = "#FFFFC8",
                                   high = "#7D0025") +
      ggplot2::labs(x = xlab,
                    y = ylab,
                    title = plot.title) +
      ggplot2::geom_label(ggplot2::aes(label = ifelse(Freq > 0,
                                                      round(Freq, 1),
                                                      NA)),
                          color = "white",
                          alpha = 0.6,
                          fill = "black")

  } else if (grepl("bar|col", type, ignore.case = TRUE)) {
    plot <- pz %>%
      as.data.frame() %>%
      ggplot2::ggplot(ggplot2::aes(Var1, Freq,
                                   fill = Var2)) +
      ggplot2::labs(x = xlab,
                    y = legend.title,
                    fill = ylab,
                    title = plot.title) +
      ggplot2::geom_col()
  }

  plot <- plot +
    theme_bw() +
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
#' @export plot.geneGating
#'

plot.geneGating <- function(object = NULL,
                            scGate.model = NULL,
                            group.by = NULL,
                            split.by = NULL) {

  if (is.null(object)) {
    stop("Please provide a Seurat object or a list of them")
  }

  if (is.list(object) && !is.null(split.by)) {
    stop("split.by only supported for a single Seurat object, not a list.\n")
  }

  if (is.null(scGate.model)) {
    stop("Please provide a scGate model or list of them")
  }

  if (!is.list(scGate.model)) {
    scGate.model <- list("scGate_model" = scGate.model)
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




#' Infer the sex of single-cell transcriptomic data based on the expression of X and Y chromosome genes
#'
#' This function takes single-cell (Seurat object or count matrix) transcriptomics data and returns the inferred sex of the cells and individuals samples provided.
#' Given the sparsity of single-cell RNA-seq data inferring of the sex in many cells is not possible (NA are returned), however, sex inference per sample is highly confident based on pseudobulk data.
#'<br>
#'Currently only suppoted for human data.
#'
#'
#' @param object Either one or a list of Seurat objects or count matrix (matrix or dcGMatrix)
#' @param infer.level Whether to infer sex for cells, samples or both.
#' @param split.by Split by sample based on a variable of metadata, either a metadata columns for Seurat object (category) or a vector indicating the sample procedence for each cell
#' @param return.Seurat return Seurat object or dataframe summarizing sex imputation. If Seurat object is provided by default Seurat object with additional metadata with infered sex on cell- and sample-wise is returned. If not Seurat object are provided, but matrices, by default dataframe with results are returned
#' @param ncores The number of cores to use, by default all available cores minus 2 are used.
#' @param bparam A \link[BiocParallel]{bpparam} object that tells Run.HiTME how to parallelize. If provided, it overrides the `ncores` parameter.
#' @param progressbar Whether to show a progress bar or not
#' @param verbose Verbose output, by default output messsage are returned

#' @import scGate
#' @importFrom Seurat CreateSeuratObject AggregateExpression SplitObject
#' @importFrom dplyr %>% left_join mutate select everything
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom scGate scGate get_scGateDB

#' @return List of genes for each GO accession requested. If named vector is provided, lists output are named as GO:accession.ID_named (e.g. "GO:0004950_cytokine_receptor_activity")
#' @export infer.Sex


infer.Sex <- function(object = NULL,
                      infer.level = c("sample", "cell", "both"),
                      split.by = NULL,
                      return.Seurat = TRUE,
                      ncores = parallel::detectCores() - 2,
                      bparam = NULL,
                      progressbar = TRUE,
                      verbose = TRUE) {

  if (is.null(object)) {
    stop("Please provide a Seurat object, count matrix or a list of them")
  }

  if (!is.character(infer.level) ||
      !tolower(infer.level[1]) %in% c("sample", "cell", "both")) {
    stop("Please check infer.level parameter.")
  }
  infer.level <- tolower(infer.level[1])

  # check if objects are Seurat or count matrices

  is_valid_object <- function(obj) {
    is_valid <- function(x) {
      inherits(x, "Seurat") || is.matrix(x) || inherits(x, "dgCMatrix")
    }

    if (is.list(obj)) {
      return(all(sapply(obj, is_valid)))
    } else {
      return(is_valid(obj))
    }
  }

  if (suppressWarnings(!is_valid_object(object))) {
    stop("Please provide a Seurat object, count matrix or a list of them")
  }


  if (!is.null(split.by) && is.list(object)) {
    stop("split.by only supported for a single object, not a list.\n
         Merge list of objects before running infer.Sex or set split.by = NULL")
  }

  if (!is.null(split.by)) {
    if (!inherits(object, "Seurat") && ncol(object) != length(split.by)) {
      stop("When providing matrix or dgCMatrix, if willing to split by sample, a vector with the same length of the number of cells (number columns of the matrix) should be provided")
    }

    if (inherits(object, "Seurat")) {
      if (!split.by %in% names(object@meta.data)) {
        stop(paste("split.by argument: ",
                   split.by, " is not a metadata column in this Seurat object"))
      }
    }
  }


  # set parallelization parameters
  param <- set_parallel_params(ncores = ncores,
                               bparam = bparam,
                               progressbar = progressbar)


  if (!inherits(object, "Seurat") && !is.list(object)) {
    object <- Seurat::CreateSeuratObject(counts = object)
    object <- Seurat::NormalizeData(object)

    # add metadata if provided
    if (!is.null(split.by)) {
      object@meta.data[["sex.splitby"]] <- split.by
      split.by <- "sex.splitby"
    }
  }

  if (!is.list(object) && !is.null(split.by)) {
    object <- Seurat::SplitObject(object,
                                  split.by = split.by)
  }



  # if object is unique, turn into Seurat if matrix, and turn into a list
  if (!is.list(object)) {
    # convert into a length=1 list to run
    object <- list(object)
  }

  # name objects
  object <- adapt_vector(object,
                         prefix = "sample_")

  # Load male and female scGate models
  suppressMessages({
    models <- scGate::get_scGateDB()
  })
  sigs <- list(Male = models$human$generic$Male,
               Female = models$human$generic$Female)

  if (infer.level == "cell" | infer.level == "both") {
    if (return.Seurat) {
      if (verbose) {message("Computing sex inference per cell...\n")}

      # Run scGate to get inference per cell
      object <- lapply(object,
                       function(s) {
                         # avoid overwritting scGate annotations
                         if ("scGate_multi" %in% names(s@meta.data)) {
                           scgate.col <- s@meta.data[["scGate_multi"]]
                         } else {
                           scgate.col <- NULL
                         }

                         s <- scGate::scGate(s,
                                             model = sigs,
                                             smooth.decay = 1,
                                             pos.thr = 0.2,
                                             maxRank = 5000,
                                             BPPARAM = param,
                                             multi.asNA = T,
                                             verbose = verbose,
                                             return.CellOntology = F)

                         s@meta.data[["sex.cell"]] <- s@meta.data[["scGate_multi"]]
                         s@meta.data[["scGate_multi"]] <- scgate.col

                         return(s)
                       })
    }
  }


  if (infer.level == "sample" | infer.level == "both") {
    if (verbose) {message("Computing sex inference per sample...\n")}

    # AggregateExpression requires some group.by
    # let's add it artificially
    if (is.null(split.by)) {
      for(a in names(object)) {
        object[[a]]@meta.data[["sex.splitby"]] <- a
      }
      split.by <- "sex.splitby"
    }


    # Compute pseudobulk
    suppressMessages({
      pseudobulk <- lapply(names(object),
                           function(s) {

                             # Compute pseudobulk and get counts
                             AggregateExpression(object[[s]],
                                                 group.by = split.by)[["RNA"]]
                           }
      )
    })

    names(pseudobulk) <- names(object)

    # infer sex on pseudobulk
    ps.sex <- lapply(pseudobulk,
                     function(s) {
                       pseudobulk_infer.Sex(s)
                     })
  }


  if (return.Seurat) {
    sex.res <- lapply(names(object),
                      function(s) {
                        seu <- object[[s]]
                        sex <- as.data.frame(ps.sex[[s]])
                        sex[[split.by]] <- s

                        seu@meta.data <- seu@meta.data %>%
                          tibble::rownames_to_column("rn") %>%
                          left_join(., sex, by = split.by) %>%
                          tibble::column_to_rownames("rn")

                        # remove extra column added for linking
                        seu@meta.data[["sex.splitby"]] <- NULL

                        return(seu)

                      })
    names(sex.res) <- names(object)

    #if list of length 1, do not return list
    if (length(sex.res) == 1) {
      sex.res <- sex.res[[1]]
    }


  } else {
    sex.res <- lapply(names(ps.sex),
                      function(s) {
                        sex <-  as.data.frame(ps.sex[[s]]) %>%
                          dplyr::mutate(Sample = s) %>%
                          dplyr::select(Sample, dplyr::everything())
                      }
    )  %>%
      data.table::rbindlist() %>%
      as.data.frame()
  }

  return(sex.res)
}


