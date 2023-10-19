ds_name <- "ZhangY_2022_34653365"
path_root <- file.path("~", "Dropbox", "CSI", "Standardized_SingleCell_Datasets", ds_name)
path_output <- file.path(path_root, "output")
cache_filename <- file.path(path_output, sprintf("%s.processed.rds", ds_name))
obj <- readRDS(cache_filename)
meta.df <- data.frame()



ref.CD8 <- load.reference.map("~/Dropbox/CSI/reference_atlases/CD8T_human_ref_v1.rds")
ref.CD4 <- load.reference.map("~/Dropbox/CSI/reference_atlases/CD4T_human_ref_v2.rds")
ref.DC <- load.reference.map("~/Dropbox/CSI/reference_atlases/DC_human_ref_v1.rds")
ref.MoMac <- load.reference.map("~/Dropbox/CSI/reference_atlases/MoMac_human_v1.rds")
maps <- c("ref.CD4","ref.CD8", "ref.DC", "ref.MoMac")

calc_CTcomp <- function(object, sample.col = "Sample", meta.df, maps = NULL){
  
  ## Broad cell types
  celltype.compositions <- list()
  celltype.compositions$sample_metadata <- meta.df
  
  broadtype_comp_table <- table(object$scGate_multi, object@meta.data[[sample.col]], useNA = "ifany")
  rownames(broadtype_comp_table)[nrow(broadtype_comp_table)] <- "Other"
  head(t(broadtype_comp_table))
  broadtype_comp_table_freq <- prop.table(broadtype_comp_table+1, margin = 2)
  
  head(t(broadtype_comp_table_freq))
  
  # clr-transform
  broadtype_comp_table_clr <- Hotelling::clr(t(broadtype_comp_table_freq))
  
  head(broadtype_comp_table_clr)
  
  celltype.compositions$broad_cell_types[["cell_counts"]] <- as.data.frame.matrix(t(broadtype_comp_table))
  celltype.compositions$broad_cell_types[["freq"]] <- as.data.frame.matrix(t(broadtype_comp_table_freq))
  celltype.compositions$broad_cell_types[["freq_clr"]] <- as.data.frame.matrix(broadtype_comp_table_clr)
  
  
  ## High-resolution subtypes

  for(map in maps){
    map.subtypes <- levels(get(map)$functional.cluster)
    subtype_comp_table <- table(object$functional.cluster[which(object$functional.cluster %in% map.subtypes)],
                                object$Sample[which(object$functional.cluster %in% map.subtypes)], useNA="always")
    if (sum(subtype_comp_table) < 10) {
      warning(paste("Using the", map, "there are less than 10 cells detected, too few to calculate a reasonable celltype composition, please check."))
      next
    }
    
    subtype_comp_table <- subtype_comp_table[-nrow(subtype_comp_table),-ncol(subtype_comp_table)] 
    
    head(t(subtype_comp_table))
    
    subtype_comp_table_freq <- prop.table(subtype_comp_table+1, margin = 2)
    
    head(t(subtype_comp_table_freq))
    
    subtype_comp_table_clr <- Hotelling::clr(t(subtype_comp_table_freq))
    
    head(subtype_comp_table_clr)
    
    map.name <- str_remove(map, "ref\\.")
    celltype.compositions$subtypes[[map.name]][["cell_counts"]] <- as.data.frame.matrix(t(subtype_comp_table))
    celltype.compositions$subtypes[[map.name]][["freq"]] <- as.data.frame.matrix(t(subtype_comp_table_freq))
    celltype.compositions$subtypes[[map.name]][["freq_clr"]] <- as.data.frame.matrix(subtype_comp_table_clr)
  }
  return(celltype.compositions)
}


# Test with large object, all celltypes present
test <- calc_CTcomp(object = obj, meta.df = meta.df, maps = maps)


# Test with object subsetted to contain only CD4T cells, most celltypes not present
# However, there seems to still be some CD8 cells around
object.CD4T <- obj.CD4T <- subset(obj, subset = scGate_multi == "CD4T")
table(object.CD4T$scGate_multi, object.CD4T@meta.data[["Sample"]], useNA = "ifany")
test <- calc_CTcomp(object = object.CD4T, meta.df = meta.df, maps = maps)