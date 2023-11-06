

# function to standarize the cell type names across studies
dict <- c("^T|CTL" = "T cell",
          "^B|B$" = "B cell",
          "^NK" = "NK",
          "Mono" = "Monocyte-like",
          "^Fibro" = "Fibroblast",
          "^Malig|cancer|tumo[u]r" = "Malignant",
          "Myeloid" = "Myeloid",
          "HSC" = "HSC",
          "DC" = "DC",
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
                            scGate.model = scGate.model){

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
             additional.signatures = "data.frame"
           )
  )

  setClass(Class = "Projectils",
           slots = setNames(rep("list", length(ref.maps)), names(ref.maps))
           )

  # build scgate S4 object
  scgate.s4 <- new(Class = "scGate",
                   UCell_score = object@meta.data[,grep("_UCell$", names(object@meta.data))],
                   scGate_is.pure = object@meta.data[,grep("^is.pure_", names(object@meta.data))],
                   scGate_multi = object@meta.data[,grep("scGate_multi", names(met2@meta.data)), drop = F]
                   )

  # Create empty s4 object to fill next
  projectils.s4 <- new(Class = "Projectils")

  slots <- c(names(ref.maps), "Consensus")

  for(i in slots){

    li <- list(Classification = object@meta.data[,grep(paste0(i,".*(_subtypes$|_confidence$])"),
                                                   names(object@meta.data)), drop = F],
           version = ref.maps[[i]]@misc$projecTILs)
    slot(projectils.s4, i) <- li
  }

  # Compute proportions




  # Compute average expression



  hit <- new("HiT",
             metadata = object@meta.data,
             Predictions = list("scGate" = scgate.s4,
                                "Projectils" = projectils.s4))


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
                                    pattern = "_subtype") {
  subtype.cols <- grep(pattern, names(object@meta.data), value = T)
  object@meta.data <- object@meta.data %>%
    mutate(Consensus_subtype = pmap(select(., all_of(subtype.cols)),
                                    consensus))
  return(object)
}


