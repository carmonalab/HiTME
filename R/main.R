# Calculates:
# - Cell counts
# - relative abundance (freq, sums to 100% for each metadata column provided)
# - clr-transformed relative abundance (freq_clr)

calc_CTcomp <- function(object, sample.col = NULL, annot.cols = "scGate_multi",
                        min.cells = 10, useNA = TRUE, rename.Multi = F){
  useNA <- ifelse(useNA == TRUE, "always", "no")
  
  celltype.compositions <- list()
  
  for (annot.col in annot.cols) {
    if (length(pbmc3k@meta.data[[annot.cols]]) == 0) {
      stop(paste("Annotation metadata column", annot.col,"could not be found"))
    }
    
    if (rename.Multi) {
      object@meta.data[[annot.col]][object@meta.data[[annot.col]] == "Multi"] <- NA
    }
    
    if (is.null(sample.col)) {
      comp_table <- table(object@meta.data[[annot.col]],
                          useNA = useNA)
    } else {
      tryCatch(invisible(length(object@meta.data[[sample.col]])), error = function(e){ stop(paste("Sample column could not be found"))})
      comp_table <- table(object@meta.data[[annot.col]],
                          object@meta.data[[sample.col]],
                          useNA = useNA)
    }

    comp_table_freq <- prop.table(comp_table+1)
    
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
