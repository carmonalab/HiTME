

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



CreateHiTObject <- function(object,
                            ref.maps = ref.maps,
                            scGate.model = scGate.model,
                            group.by = c("scGate_multi",
                                         "Consensus_subtypes"),
                            annot.cols.freq = c("scGate_multi",
                                                "Consensus_subtypes"),
                            additional.signatures = NULL){

  ## Build S4 objects to store data
  setClass(Class = "HiT",
           slots = list(
             metadata = "data.frame",
             Predictions = "list",
             Proportions = "list",
             Avg.Expression = "list",
             version = "list"
           )
  )


  setClass(Class = "scGate",
           slots = list(
             UCell_score = "data.frame",
             scGate_is.pure = "data.frame",
             scGate_multi = "data.frame",
             additional.signatures = "ANY"
           )
  )

  # add slot for consensus
  ref.maps[["Consensus"]] <- NA
  setClass(Class = "Projectils",
           slots = setNames(rep("list", length(ref.maps)), names(ref.maps))
           )

  if(is.null(additional.signatures)){
    ucell_score <- object@meta.data[,grep("_UCell$", names(object@meta.data))]
    add.sig = NA
  } else{
    all.ucell <- object@meta.data[,grep("_UCell$", names(object@meta.data))]
    sig <- grep(paste(names(additional.signatures), collapse = "|"),
                names(all.ucell), value = T)
    ucell <- grep(paste(names(additional.signatures), collapse = "|"),
                  names(all.ucell), value = T, invert = T)

    ucell_score <- all.ucell[,ucell]
    add.sig <- all.ucell[,sig]
  }


  # build scgate S4 object
  scgate.s4 <- new(Class = "scGate",
                   UCell_score = ucell_score,
                   scGate_is.pure = object@meta.data[,grep("^is.pure_", names(object@meta.data))],
                   scGate_multi = object@meta.data[,grep("scGate_multi", names(object@meta.data)), drop = F],
                   additional.signatures = add.sig
                   )

  # Create empty s4 object to fill next
  projectils.s4 <- new(Class = "Projectils")

  slots <- names(ref.maps)

  for(i in slots){
    if(!is.na(ref.maps[[i]])){
    version <- ref.maps[[i]]@misc$projecTILs
  } else {
    version <- NA
  }
    li <- list(Classification = object@meta.data[,grep(paste0(i,".*(_subtypes$|_confidence$])"),
                                                   names(object@meta.data)), drop = F],
           version = version)
    slot(projectils.s4, i) <- li
  }

  # Compute proportions
  comp.prop <- calc_CTcomp(object,
                           annot.cols = annot.cols.freq
                            )
  # Compute avg expressoin
  avg.expr <- CellAvgExpression(object,
                                group.by = "annotation")


  hit <- new("HiT",
             metadata = object@meta.data,
             Predictions = list("scGate" = scgate.s4,
                                "Projectils" = projectils.s4),
             Avg.Expression = avg.expr,
             Proportions = comp.prop
             )


  return(hit)


}




### Change names of added projectils classification
addProjectilsClassification <- function(object,
                                        ref.map.name){
  object@meta.data[[paste0(ref.map.name, "_subtypes")]] <- object@meta.data[["functional.cluster"]]
  object@meta.data[["functional.cluster"]] <- NULL
  object@meta.data[[paste0(ref.map.name, "_confidence")]] <- object@meta.data[["functional.cluster.conf"]]
  object@meta.data[["functional.cluster.conf"]] <- NULL

  return(object)
}


### Create consensus from the various projectils classification maps

# Define the custom logic function
consensus <- function(...) {
  if (all(is.na(c(...)))) {
    return(NA)
  } else if (sum(!is.na(c(...))) == 1) {
    return(c(...)[!is.na(c(...))])
  } else {
    return("Multiple")
  }
}

classificationConsensus <- function(object,
                                    pattern = "_subtypes") {
  subtype.cols <- grep(pattern, names(object@meta.data), value = T)
  object@meta.data <- object@meta.data %>%
    mutate(Consensus_subtypes = pmap(select(., all_of(subtype.cols)),
                                    consensus),
           Consensus_subtypes = sapply(Consensus_subtypes, function(x){
                                        unname(unlist(x))
           }))
  return(object)
}


### Function to compute average expression

CellAvgExpression <- function(object,
                              group.by = group.by,
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

  if(is.null(group.by)){
    message("No grouping provided, grouping by consensus of Projectils\n")
    group.by <- "Consensus_subtypes"
  } else if (!is.vector(group.by)) {
    stop("Please provide one or various grouping variables as a vector")
  } else {
    group.by <- unique(c(group.by, "Consensus_subtypes"))
  }

  if(any(!group.by %in% names(object@meta.data))){
    g.missing <- group.by[!group.by %in% names(object@meta.data)]
    warning(paste(g.missing, "not in object metadata, it is not computed for average expression.\n"))
    group.by <- group.by[group.by %in% names(object@meta.data)]
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

  for(i in group.by){
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


