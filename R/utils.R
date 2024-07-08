# ProjecTILs on multiple reference maps #############################################################

ProjecTILs.classifier.multi <- function(object,
                                        ref.maps,
                                        bparam = NULL,
                                        layer1_link = "CellOntology_ID") {

  # count how many cells are found by scGate for each map reference
  if (layer1_link %in% names(object@meta.data)) {
    filter.cells <- F
    map.celltypes <- lapply(ref.maps,
                            function(x) {
                              ct <- x@misc$layer1_link
                              nrow(object@meta.data[object@meta.data[[layer1_link]] %in% ct,])
                            })

    present <- names(ref.maps[map.celltypes > 1])
    no.present <- names(ref.maps)[!names(ref.maps) %in% present]

    if (length(present) == 0) {
      message(paste("No cells linked to reference maps found by", layer1_link,
                    ".\nNot running Projectils"))
      run <- FALSE
    } else {
      run <- TRUE
      message("Performing Projectils classification for ",
              paste(present, collapse = ", "))

      if (length(ref.maps) != length(present)) {
        message("Not doing mapping for ",
                paste(no.present, collapse = ", ") )
        ref.maps <- ref.maps[present]
      }
    }

  } else {
    filter.cells <- TRUE

    message("Not scGate classification found in this object\n",
            "Running default scGate (or layer1 classification) filtering within Projectils)")
    run <- TRUE
  }

  if (run) {
    suppressWarnings(require(ProjecTILs))

    functional.clusters <-
      BiocParallel::bplapply(
        X = names(ref.maps),
        BPPARAM = bparam,
        function(m) {
          if (filter.cells) {
            # don't subset if no scGate gating is found in metadata, instead run scGate within projectils
            subset.object <- object
          } else {
            map.celltype <- ref.maps[[m]]@misc$layer1_link
            subset.object <- object[,object@meta.data[[layer1_link]] %in% map.celltype]

          }

          if (ncol(subset.object)>0) {
            message("\nRunning Projectils for ", m, " reference")
            # it is mandatory to make run in serial (ncores = 1 and BPPARAM SerialParam)
            subset.object.pred <- ProjecTILs:::classifier.singleobject(subset.object,
                                                                       ref = ref.maps[[m]],
                                                                       filter.cells = filter.cells,
                                                                       ncores = 1)

            return(subset.object.pred)

          } else {
            message("Not Running Projectils classifier for reference map",
                    paste(map.celltype, collapse = ", "), ". No cells found")
          }
        }
      )


    if (!all(unlist(lapply(functional.clusters, function(x) {is.null(x)})))) {
      functional.clusters <- data.table::rbindlist(lapply(functional.clusters,
                                                          data.table::setDT,
                                                          keep.rownames = TRUE)) %>%
        #remove duplicated row.names (classified by different ref.maps)
        dplyr::filter(!duplicated(rn)) %>%
        tibble::column_to_rownames("rn")

      # remove already present functional.cluster columns
      rm <- grep("^functional.cluster",
                 names(object@meta.data))
      if (length(rm) > 0) {
        object@meta.data <- object@meta.data[, -rm, drop = F]
      }
      object@meta.data <- merge(object@meta.data,
                                functional.clusters,
                                by = 0,
                                all.x = T) %>%
        tibble::column_to_rownames("Row.names")

      # save names of reference maps run (will be lost if objects remerged)
      object@misc[["layer2_param"]][["functional.cluster"]][["References_executed"]] <- names(ref.maps)

    } else {
      object@meta.data <- object@meta.data %>%
        dplyr::mutate(functional.cluster = NA,
                      functional.cluster.conf = NA)
    }
  }
  return(object)
}




match_dictionary <- function(cell_type,
                             dictionary = NULL) {

  if (is.null(dictionary)) {
    dict <- c("^T|CTL" = "T cell",
              "^B|B$" = "B cell",
              "^NK" = "NK",
              "Mono" = "Monocyte-like",
              "^Fibro" = "Fibroblast",
              "^Malig|cancer|tumo[u]r" = "Malignant",
              "Myeloid" = "Myeloid",
              "HSC" = "HSC",
              "DC|Dendri" = "DC",
              "^Epi|^EC" = "Epithelial",
              "Mast" = "Mastocyte",
              "Macroph" = "Macrophage",
              "^Endo" = "Endothelial",
              "Neutro" = "Neutrophil",
              "Strom" = "Stromal")
  }

  for (keyword in names(dict)) {
    if (grepl(keyword, cell_type, ignore.case = TRUE)) {
      cell_type <- dict[[keyword]]
      break
    }
  }
  return(cell_type)
}



StandardizeCellNames <- function(cell.names, dictionary = NULL) {
  standarized_cellnames <- sapply(cell.names, match_dictionary, dictionary)
  standarized_cellnames <- unname(standarized_cellnames)

  return(standarized_cellnames)
}



set_parallel_params <- function(ncores,
                                bparam,
                                progressbar)
{
  if (is.null(ncores)) {
    ncores <- 1
  }

  if (ncores > parallelly::availableCores()) {
    ncores <- parallelly::availableCores()
    message("Using more cores available in this computer, reducing number of cores to ", ncores)
  }

  # set parallelization parameters
  if (is.null(bparam)) {
    if (ncores > 1) {
      param <- BiocParallel::MulticoreParam(workers =  ncores,
                                            progressbar = progressbar)
    } else {
      param <- BiocParallel::SerialParam()
    }
  } else {
    param <- bparam
  }
  return(param)
}
