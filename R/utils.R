#############################################################################################################
# HiT class definition
## Build S4 objects to store data
methods::setClass(Class = "HiT",
                  slots = list(
                    metadata = "data.frame",
                    predictions = "list",
                    composition = "list",
                    aggregated_profile = "list",
                    version = "list"
                  )
)


#############################################################################################################
# Projectils on multiple reference maps
ProjecTILs.classifier.multi <- function(object,
                                        ref.maps,
                                        bparam = NULL,
                                        layer1_link = "CellOntology_ID"){

  # count how many cells are found by scGate for each map reference
  if(layer1_link %in% names(object@meta.data)){
    filter.cells <- F
    map.celltypes <- lapply(ref.maps,
                            function(x){
                              ct <- x@misc$layer1_link
                              # collpase if multiple cell types, such as MoMac
                              if(length(ct)>1){ct <- paste(ct, collapse = "_")}

                              nrow(object@meta.data[object@meta.data[[layer1_link]] %in% ct,])
                            })
    present <- names(ref.maps[map.celltypes > 0])
    no.present <- names(ref.maps)[!names(ref.maps) %in% present]

    if(length(present) == 0){
      message(paste("No cells linked to reference maps found by", layer1_link, ".\nNot running Projectils"))
      run <- FALSE
    } else {
      run <- TRUE
      message("Performing Projectils classification for ", paste(present, collapse = ", "))

      if(length(ref.maps) != length(present)){
        message("Not doing mapping for ", paste(no.present, collapse = ", ") )
        ref.maps <- ref.maps[present]
      }
    }

  } else {
    filter.cells <- T

    message("Not scGate classification found in this object\n",
            "Running default scGate (or layer1 classification) filtering within Projectils)")
    run <- TRUE
  }

  if(run){
  functional.clusters <-
    BiocParallel::bplapply(
      X = names(ref.maps),
      BPPARAM = bparam,
      FUN = function(m){
        if(filter.cells){
          # don't subset if no scGate gating is found in metadata, instead run scGate within projectils
          subset.object <- object
        } else {
        map.celltype <- ref.maps[[m]]@misc$layer1_link
        subset.object <- object[,object@meta.data[[layer1_link]] %in% map.celltype]

        }

        if(ncol(subset.object)>0){
          message("\nRunning Projectils for ", m, " reference")
          data("Hs2Mm.convert.table")
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


  if(!all(unlist(lapply(functional.clusters, function(x){is.null(x)})))){
    functional.clusters <- data.table::rbindlist(lapply(functional.clusters,
                                                        data.table::setDT,
                                                        keep.rownames = TRUE)) %>%
                            #remove duplicated rownames (classified by different ref.maps)
                            dplyr::filter(!duplicated(rn)) %>%
                            tibble::column_to_rownames("rn")
    object@meta.data <- merge(object@meta.data, functional.clusters, by = 0, all.x = T) %>%
                        tibble::column_to_rownames("Row.names")

    # save names of reference maps run
    object@misc[["Projectils_param"]] <- names(ref.maps)
  } else {
    object@meta.data <- object@meta.data %>%
                        dplyr::mutate(functional.cluster = NA,
                                      functional.cluster.conf = NA)
  }
  }


  return(object)

}


#############################################################################################################
# function to match words with grep
match_dictionary <- function(cell_type,
                             dictionary = NULL) {

  if(is.null(dictionary)){
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
              "Strom" = "Stromal"
    )
  }

  for (keyword in names(dict)) {
    if (grepl(keyword, cell_type, ignore.case = TRUE)) {
      cell_type <- dict[[keyword]]
      break
    }
  }
  return(cell_type)
}

#############################################################################################################
# complete function applying previous function with sapply to all vector
StandardizeCellNames <- function(cell.names, dictionary = NULL){

  standarized_cellnames <- sapply(cell.names,
                                  match_dictionary,
                                  dictionary)
  standarized_cellnames <- unname(standarized_cellnames)
  return(standarized_cellnames)

}



