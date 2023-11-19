# projectils on multiple reference maps
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

    message("Performing Projectils classification for ", paste(present, collapse = ", "))

    if(length(ref.maps) != length(present)){
      message("Not doing mapping for ", paste(no.present, collapse = ", ") )
      ref.maps <- ref.maps[present]
    }
  } else {
    filter.cells <- T

    message("Not scGate classification found in this object\n",
            "Running default scGate (or layer1 classification) filtering within Projectils)")
  }

  if(length(present) == 0){
    warning(paste("No cells linked to reference maps found by", layer1_link, ". Not running Projectils."))
  } else {


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
          message("\nRunning Projectils for", paste(map.celltype, collapse = ", "), " celltypes")
          }

          if(ncol(subset.object)>0){

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


    functional.clusters <- data.table::rbindlist(lapply(functional.clusters,
                                                        setDT, keep.rownames = TRUE)) %>%
                            tibble::column_to_rownames("rn")
    object@meta.data <- merge(object@meta.data, functional.clusters, by = 0, all.x = T) %>%
                        tibble::column_to_rownames("Row.names")

    # save names of reference maps run
    object@misc[["Projectils_param"]] <- names(ref.maps)
  }


  return(object)

}



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

# complete function applying previous function with sapply to all vector
StandardizeCellNames <- function(cell.names, dictionary = NULL){

  standarized_cellnames <- sapply(cell.names,
                                  match_dictionary,
                                  dictionary)
  standarized_cellnames <- unname(standarized_cellnames)
  return(standarized_cellnames)

}


### Function to get the list of certain genes
GetGOList <- function(GO_accession = NULL,
                      species = "Human",
                      host = "https://dec2021.archive.ensembl.org/"
){

  GO_acc <- c(
    "DNA-binding transcription factor activity" = "GO:0003700",
    "Cytokine activity" = "GO:0005125",
    "cytokine receptor activity" = "GO:0004896",
    "chemokine activity" = "GO:0008009",
    "chemokine receptor activity" = "GO:0004950"
  )

  if(!is.null(GO_accession)){
    if(length(names(GO_accession)) == 0){
      names(GO_accession) <- GO_accession
    }
    GO_acc <- c(GO_acc, GO_accession)
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

  # Retrieve genes for each GO
  ensembl = biomaRt::useMart("ensembl",
                             dataset = dataset)

  gene.data.bm <- biomaRt::getBM(attributes=c('hgnc_symbol',
                                              'go_id'),
                                 filters = 'go',
                                 values = GO_acc,
                                 mart = ensembl)
  gene.data <- gene.data.bm %>%
    filter(go_id %in% GO_acc) %>%
    filter(hgnc_symbol != "") %>%
    split(., as.factor(.$go_id)) %>%
    lapply(., function(x) x[,1])

  if(length(names(gene.data)) != length(GO_acc)){
    warning("Additional GO accession provided not found in GO database")
  }
  names(gene.data) <- paste0(names(gene.data), "_", names(GO_acc)[GO_acc %in% names(gene.data)])


  return(gene.data)

}
