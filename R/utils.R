# ProjecTILs on multiple reference maps #############################################################

ProjecTILs.classifier.multi <- function(object,
                                        ref.maps,
                                        bparam = NULL,
                                        layer1_link = "CellOntology_ID",
                                        verbose = TRUE) {

  # count how many cells are found by scGate for each map reference
  if (layer1_link %in% names(object@meta.data)) {
    filter.cells <- F
    map.celltypes.count <- lapply(ref.maps,
                                  function(x) {
                                    ct <- x@misc$layer1_link
                                    nrow(object@meta.data[object@meta.data[[layer1_link]] %in% ct,])
                                  })

    present <- names(ref.maps[map.celltypes.count > 1])
    no.present <- names(ref.maps)[!names(ref.maps) %in% present]

    if (length(present) == 0) {
      if(verbose){message(paste("No cells linked to reference maps found by",
                                layer1_link,
                                ".\nNot running Projectils"))}
      run <- FALSE
    } else {
      run <- TRUE
      if(verbose){message("Performing Projectils classification for ",
                          paste(present, collapse = ", "))}

      if (length(ref.maps) != length(present)) {
        if(verbose){message("Not doing mapping for ",
                            paste(no.present, collapse = ", ") )}
        ref.maps <- ref.maps[present]
      }
    }

  } else {
    filter.cells <- TRUE

    if(verbose){message("Not scGate classification found in this object\n",
                        "Running default scGate (or layer1 classification) filtering within Projectils for all provided reference maps)")}
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
            if(verbose){message("\nRunning Projectils for ", m, " reference")}
            # it is mandatory to make run in serial (ncores = 1 and BPPARAM SerialParam)
            subset.object.pred <- ProjecTILs:::classifier.singleobject(subset.object,
                                                                       ref = ref.maps[[m]],
                                                                       filter.cells = filter.cells,
                                                                       ncores = 1)

            return(subset.object.pred)

          } else {
            if(verbose){message("Not Running Projectils classifier for reference map",
                                paste(map.celltype, collapse = ", "), ". No cells found")}
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



set_parallel_params <- function(ncores = NULL,
                                bparam = NULL,
                                progressbar = TRUE)
{
  if (is.null(ncores)) {
    ncores <- 1
  }

  if (ncores > parallel::detectCores()) {
    ncores <- parallel::detectCores()
    message("Using more cores available in this computer, reducing number of cores to ",
            ncores)
  }

  # set parallelization parameters
  if (is.null(bparam)) {
    if (ncores > 1) {
      if(.Platform$OS.type == "windows"){
        param <- BiocParallel::SnowParam(workers=ncores,
                                         progressbar = progressbar)
      } else {
        param <- BiocParallel::MulticoreParam(workers =  ncores,
                                              progressbar = progressbar)
      }
    } else {
      param <- BiocParallel::SerialParam(progressbar = progressbar)
    }
  } else {
    param <- bparam
  }
  return(param)
}


# adapt vectors or list for signatures
adapt_vector <- function(vector,
                         prefix){

  if (!is.null(vector)) {
    # convert to list if not
    if (!is.list(vector)) {
      vector <- list(vector)
    }

    # if some elements of list is not named, name it
    for (v in seq_along(vector)) {
      if (is.null(names(vector)[[v]]) ||
          is.na(names(vector)[[v]]) ||
          names(vector)[[v]] == "") {
        names(vector)[[v]] <- paste0(prefix, v)
      }
    }
  }
  return(vector)
}


# combine different tags in different columns for each row
combine_tags <- function(row,
                         default = "resting",
                         collapse = NULL){
  r <- unique(row)
  # remove default value, only return if all are default
  if(length(r)>1){
    r <- r[r != default]
  }

  # option to collapse
  if(!is.null(collapse)){
    r <- paste(r,
               collapse = collapse)
  }

  return(r)
}


# accommodate layer3 classification
get.layer3 <- function(s,
                       ann.col = "layer2",
                       sigs.cols,
                       layer3.threshold = 0.2,
                       default.label = "resting"){

  s@meta.data <- s@meta.data %>%
    tibble::rownames_to_column("cell.id") %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(sigs.cols),
                         ~ ifelse(. > layer3.threshold,
                                  gsub("_UCell", "", cur_column()),
                                  default.label),
                         .names = "is.pureL3_{col}")
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      layer3_annotation = combine_tags(dplyr::c_across(dplyr::starts_with("is.pureL3_")),
                                       collapse = "_",
                                       default = default.label),
      layer3 = ifelse(!is.na(.data[[ann.col]]),
                      paste(c(as.character(.data[[ann.col]]),
                              layer3_annotation),
                            collapse = "_"),
                      NA)
    ) %>%
    dplyr::ungroup() %>%
    tibble::column_to_rownames("cell.id")

  return(s)
}


# takes pseudobulk or bulk matrix and make ratios of sex genes
pseudobulk_infer.Sex <- function(matrix,
                                 male.gene = "RPS4Y1",
                                 female.gene = "XIST",
                                 log2.threshold = 1){
  female.counts <- 1
  male.counts <- 1

  # add 1 to mantain ratio in case not found
  try({
    female.counts <- matrix[female.gene,] + 1
  })
  try({
    male.counts <- matrix[male.gene,] + 1
  })

  marker_ratio <- log2(female.counts/male.counts)

  sample.sex <- ifelse(abs(marker_ratio) < log2.threshold,
                       NA,
                       ifelse(marker_ratio > 0,
                              "female",
                              "male"
                       )
  )

  ret <- list("A" = male.counts,
              "B" = female.counts,
              log2_female.male_ratio = marker_ratio,
              sample.sex = sample.sex)

  names(ret)[1:2] <- c(paste(male.gene, "pseudocounts"),
                  paste(female.gene, "pseudocounts"))

  return(ret)

}


