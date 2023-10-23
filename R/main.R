# Annotate cell types in parallel
annotate_cells <- function(dir, scGate.model, ref.maps,
                           ncores = 1, progressbar = T){
  files <- list.files(dir)
  if (!all(endsWith(files, '.rds'))) {
    stop(paste("There are some files not ending on '.rds'"))
  }
  
  param <- BiocParallel::MulticoreParam(workers = ncores, progressbar = progressbar)
  
  print("Running scGate")
  BiocParallel::bplapply(
    X = files,
    BPPARAM =  param,
    FUN = function(file) {
         path <- file.path(dir, file)
         x <- readRDS(path)
         x <- scGate(x, model=models.TME)
         saveRDS(x, path)
      }
  )
  print("Finished scGate")
  
  for (i in 1:length(ref.maps)) {
    ref.map.name <- names(ref.maps)[i]
    print(paste("Running ProjecTILs with map:  ", ref.map.name))
    BiocParallel::bplapply(
      X = files, 
      BPPARAM =  param,
      FUN = function(file) {
        path <- file.path(dir, file)
        x <- readRDS(path)
        x <- ProjecTILs.classifier(x, ref.maps[[i]])
        x@meta.data[[paste0(ref.map.name, ".subtypes")]] <- x@meta.data[["functional.cluster"]]
        x@meta.data[["functional.cluster"]] <- NULL
        saveRDS(x, path)
      }
    )
    print(paste("Finished ProjecTILs with map:  ", ref.map.name))
  }
}


# Calculate cell type composition
# - Cell counts
# - relative abundance (freq, sums to 100% for each metadata column provided)
# - clr-transformed relative abundance (freq_clr)

calc_CTcomp <- function(object, sample.col = NULL, annot.cols = "scGate_multi",
                        min.cells = 10, useNA = TRUE, rename.Multi = F){
  useNA <- ifelse(useNA == TRUE, "always", "no")
  
  celltype.compositions <- list()
  
  for (annot.col in annot.cols) {
    if (length(object@meta.data[[annot.cols]]) == 0) {
      stop(paste("Annotation metadata column", annot.col,"could not be found"))
    }
    
    if (rename.Multi) {
      object@meta.data[[annot.col]][object@meta.data[[annot.col]] == "Multi"] <- NA
    }
    
    if (is.null(sample.col)) {
      comp_table <- table(object@meta.data[[annot.col]],
                          useNA = useNA)
      comp_table_freq <- prop.table(comp_table+1) * 100 # To get percentage
    } else {
      if (length(object@meta.data[[sample.col]]) == 0) {
        stop(paste("Sample column could not be found"))
      }
      comp_table <- table(object@meta.data[[annot.col]],
                          object@meta.data[[sample.col]],
                          useNA = useNA)
      comp_table_freq <- prop.table(comp_table+1, margin = 2) * 100 # To get percentage
    }

    if (sum(comp_table) < min.cells) {
      warning(paste("There are less than", min.cells,
                    "cells detected. This is too few to calculate a reasonable celltype composition, please check or adjust parameter min.cells."))
      next
    }
    
    ## clr-transform
    comp_table_clr <- Hotelling::clr(t(comp_table_freq))
    
    ## Append
    celltype.compositions[[annot.col]][["cell_counts"]] <- as.data.frame.matrix(t(comp_table))
    celltype.compositions[[annot.col]][["freq"]] <- as.data.frame.matrix(t(comp_table_freq))
    celltype.compositions[[annot.col]][["freq_clr"]] <- as.data.frame.matrix(comp_table_clr)
  }
  return(celltype.compositions)
}
