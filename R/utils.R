# projectils on multiple reference maps
ProjecTILs.classifier.multi <- function(object,
                                        ref.maps,
                                        param = param,
                                        layer1 = "scGate_multi"){

  # count how many cells are found by scGate for each map reference
  if(layer1 %in% names(object@meta.data)){
    filter.cells <- F
    map.celltypes <- lapply(ref.maps,
                            function(x){
                              ct <- x@misc$scGate_link
                              nrow(object@meta.data[object@meta.data[[layer1]] %in% ct,])
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

  if(length(ref.maps) == 0){
    stop(paste("No cells linked to reference maps found by", layer1, ". Not running Projectils."))
  }

  functional.clusters <-
    BiocParallel::bplapply(
      X = names(ref.maps),
      BPPARAM = param,
      FUN = function(m){
        map.celltype <- ref.maps[[m]]@misc$scGate_link
        subset.object <- object[,object@meta.data[[layer1]] %in% map.celltype]

        if(ncol(subset.object)>0){
          message("\nRunning Projectils for", paste(map.celltype, collapse = ", "), " celltypes")

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


  return(object)

}






# function to standarize the cell type names across studies
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

# function to match words with grep
match_dictionary <- function(cell_type,
                             dict) {
  for (keyword in names(dict)) {
    if (grepl(keyword, cell_type, ignore.case = TRUE)) {
      cell_type <- dict[[keyword]]
      break
    }
  }
  return(cell_type)
}

# complete function applying previous function with sapply to all vector
StandardizeCellNames <- function(cell.names,
                                 dictionary = dict){

  standarized_cellnames <- sapply(cell.names,
                                  match_dictionary,
                                  dict)
  standarized_cellnames <- unname(standarized_cellnames)
  return(standarized_cellnames)

}



###############################################################



get.HiTObject <- function(object,
                            group.by = c("layer1",
                                         "layer2"),
                            group.by.aggregated = NULL,
                            group.by.composition = NULL){


  if (is.null(object)) {
    stop("Please provide either a Seurat object")
  } else if(class(object) != "Seurat"){
    stop("Not Seurat object included, cannot be processed.\n")
  }

  if(is.null(group.by.aggregated)){
    group.by.aggregated <- group.by
  }

  if(is.null(group.by.composition)){
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


  setClass(Class = "layer1",
           slots = list(
             UCell_score = "data.frame",
             scGate_is.pure = "data.frame",
             scGate_multi = "data.frame",
             additional.signatures = "ANY"
           )
  )


  setClass(Class = "layer2",
           slots = list(
             functional.cluster = "data.frame"
           )
           )

  # extract values from misc slot from object
  additional.signatures <- object@misc$scGate_param$scGate_additional.signatures
  scgate.models <- object@misc$scGate_param$scGate_models

    sig <- grep(paste(additional.signatures, collapse = "|"),
                names(object@meta.data), value = T)
    sig.df <- object@meta.data[, sig]
    scgate <- grep(paste(scgate.models, collapse = "|"),
                  names(object@meta.data), value = T)
    scgate.df <- object@meta.data[, scgate]

    ucell <- names(object@meta.data)[!names(object@meta.data) %in% c(sig, scgate)] %>%
                grep("_UCell$", ., value = T)

    ucell.df <- object@meta.data[, ucell]




  # build scgate S4 object
  layer1.s4 <- new(Class = "layer1",
                   UCell_score = ucell.df,
                   scGate_is.pure = scgate.df,
                   scGate_multi = object@meta.data[,grep("scGate_multi", names(object@meta.data)), drop = F],
                   additional.signatures = sig.df
                   )

  # Create empty s4 object to fill next
  layer2.s4 <- new(Class = "layer2",
                       functional.cluster = object@meta.data[,grep("functional.cluster", names(object@meta.data)), drop = F])


  # Compute proportions
  comp.prop <- get.celltype.composition(object,
                                       group.by.composition = group.by.composition
                                      )
  # Compute avg expressoin
  avg.expr <- get.aggregated.profile(object,
                                group.by.aggregated = group.by.aggregated)


  hit <- new("HiT",
             metadata = object@meta.data,
             predictions = list("layer1" = layer1.s4,
                                "layer2" = layer2.s4),
             aggregated_profile = avg.expr,
             composition = comp.prop
             )


  return(hit)


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

    if(species == "human") {
      dataset <- "hsapiens_gene_ensembl"
    } else if(species == "mouse"){
      dataset <- "mmusculus_gene_ensembl"
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


