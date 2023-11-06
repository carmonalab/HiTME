

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
