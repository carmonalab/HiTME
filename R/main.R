# Extracts the broad cell type and subtype composition from scGate and ProjecTIL annotated datasets
# Calculates relative abundance (freq) and clr-transformed freq

aux <- file.path(getwd(), "aux")

## Use this to add new maps with new subtypes
# CD8T_human_ref_v1 <- load.reference.map("~/Dropbox/CSI/reference_atlases/CD8T_human_ref_v1.rds")
# CD4T_human_ref_v2 <- load.reference.map("~/Dropbox/CSI/reference_atlases/CD4T_human_ref_v2.rds")
# DC_human_ref_v1 <- load.reference.map("~/Dropbox/CSI/reference_atlases/DC_human_ref_v1.rds")
# MoMac_human_v1 <- load.reference.map("~/Dropbox/CSI/reference_atlases/MoMac_human_v1.rds")
# map.subtypes <- c()
# for(map in maps){
#   map.subtypes <- c(levels(get(map)$functional.cluster))
#   # Need to add "\n" (line break) at the end, otherwise read.csv will complain
#   cat(c(paste0(map.subtypes, collapse = ","), "\n"), file = file.path(aux, sprintf("%s.txt", map)), sep = "")
# }

# Default maps
m <- c("CD8T_human_ref_v1","CD4T_human_ref_v2", "DC_human_ref_v1", "MoMac_human_v1")

calc_CTcomp <- function(object, sample.col = "Sample", meta.df, maps = m,
                        min.cells = 10, useNA = TRUE){
  # Check maps input
  maps.available <- stringr::str_remove(list.files(aux), ".txt")
  if (any(!maps %in% maps.available)) {
    stop(paste0(maps[!maps %in% maps.available], " not found in available reference maps. Available reference maps:\n", paste(maps.available, collapse = "\n")))
  }
  
  # Broad cell types
  celltype.compositions <- list()
  celltype.compositions$sample_metadata <- meta.df
  
  useNA <- ifelse(useNA == TRUE, "always", "no")
  object$scGate_multi[object$scGate_multi == "Multi"] <- NA
  broadtype_comp_table <- table(object$scGate_multi,
                                object@meta.data[[sample.col]],
                                useNA = useNA)
  broadtype_comp_table_freq <- prop.table(broadtype_comp_table+1, margin = 2)
  
  ## clr-transform
  broadtype_comp_table_clr <- Hotelling::clr(t(broadtype_comp_table_freq))
  
  ## Append
  celltype.compositions$broad_cell_types[["cell_counts"]] <- as.data.frame.matrix(t(broadtype_comp_table))
  celltype.compositions$broad_cell_types[["freq"]] <- as.data.frame.matrix(t(broadtype_comp_table_freq))
  celltype.compositions$broad_cell_types[["freq_clr"]] <- as.data.frame.matrix(broadtype_comp_table_clr)
  
  
  # High-resolution subtypes
  
  ## Get subtypes for specific reference map version
  map.subtypes <- list()
  for(map in maps){
    map.subtypes[[map]] <- unlist(read.csv(file.path(aux, sprintf("%s.txt", map)), header = F))
  }
  for(map in maps){
    print(paste("Used subtypes from reference map:  ", map))
    subtypes <- map.subtypes[[map]]
    subtype_comp_table <- table(object$functional.cluster[which(object$functional.cluster %in% subtypes)],
                                object$Sample[which(object$functional.cluster %in% subtypes)],
                                useNA = useNA)
    subtype_comp_table <- subtype_comp_table[-nrow(subtype_comp_table),-ncol(subtype_comp_table)] 
    
    if (sum(subtype_comp_table) < min.cells) {
      warning(paste("Using the", map, "there are less than", min.cells, "cells detected, too few to calculate a reasonable celltype composition, please check."))
      next
    }
    
    subtype_comp_table_freq <- prop.table(subtype_comp_table+1, margin = 2)
    
    ## clr-transform
    subtype_comp_table_clr <- Hotelling::clr(t(subtype_comp_table_freq))
    
    ## Append
    map.name <- strsplit(map, "_")[[1]][1]
    celltype.compositions$subtypes[[map.name]][["cell_counts"]] <- as.data.frame.matrix(t(subtype_comp_table))
    celltype.compositions$subtypes[[map.name]][["freq"]] <- as.data.frame.matrix(t(subtype_comp_table_freq))
    celltype.compositions$subtypes[[map.name]][["freq_clr"]] <- as.data.frame.matrix(subtype_comp_table_clr)
  }
  return(celltype.compositions)
}
