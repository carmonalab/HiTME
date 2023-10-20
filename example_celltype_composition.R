ds_name <- "ZhangY_2022_34653365"
path_root <- file.path("~", "Dropbox", "CSI", "Standardized_SingleCell_Datasets", ds_name)
path_output <- file.path(path_root, "output")
cache_filename <- file.path(path_output, sprintf("%s.processed.rds", ds_name))
obj <- readRDS(cache_filename)
df <- data.frame()

source("R/main.R")

# Test with large object, all celltypes present
test <- calc_CTcomp(object = obj, meta.df = df)
no_NA_test <- calc_CTcomp(object = obj, meta.df = meta.df, useNA = FALSE)
no_maps_test <- calc_CTcomp(object = obj, meta.df = meta.df, maps = NULL)
bad_maps_test <- calc_CTcomp(object = obj, meta.df = df, maps = c(m, "abc"))

# Test with object subsetted to contain only CD4T cells, most celltypes not present
# However, there seems to still be some CD8 cells around
object.CD4T <- obj.CD4T <- subset(obj, subset = scGate_multi == "CD4T")
CD4_test <- calc_CTcomp(object = object.CD4T, meta.df = meta.df)