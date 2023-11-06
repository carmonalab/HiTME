

# function to standarize the cell type names across studies
dict <- c("^T" = "T cell",
          "^B|B$" = "B cell",
          "^NK" = "NK",
          "Mono" = "Monocyte-like",
          "^Fibro" = "Fibroblast",
          "^Malig|cancer" = "Malignant",
          "Myeloid" = "Myeloid",
          "HSC" = "HSC",
          "DC" = "DC",
          "^Epithe|^EC" = "Epithelial",
          "Mast" = "Mastocyte",
          "Macroph" = "Macrophage",




          )



Standarize_celltypes <- function(cell.names){

}
