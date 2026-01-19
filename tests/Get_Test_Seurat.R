require(Seurat)
require(dplyr)

# Get as light as possible Seurat object for testing

# we will use Bassez object found in https://carmonalab.github.io/HiTME_CaseStudies/HiTME_GetStarted.html
# Can be downloaded from (https://figshare.com/ndownloader/files/43848153). This includes 3 samples in a list.

#### PARAMS #####
path <- "~/Documents/Projects/HiTME/tests/BassezA_2021_33958794_3patients.rds" # path to Seurat object before filtering
sample_id <- "sample"
ncells <- 100
myseed <- 22

##########################

# load object
seu <- readRDS(path)

# split into samples
seu <- SplitObject(seu,
                   split.by = sample_id)
lapply(seu, dim)

# Use Seurat's sketch to reduce number of cells
slim <- lapply(seu,
               function(x){
                 # process
                 x <- NormalizeData(x,
                                    verbose = F)
                 x <- FindVariableFeatures(x,
                                           selection.method = "vst",
                                           nfeatures = 2000,
                                           verbose = F)
                 x<- ScaleData(x,
                                 vars.to.regress = NULL,
                                 assay = "RNA",
                                 verbose = F)
                 x <- RunPCA(x,
                               npcs=30,
                               verbose = F)
                 # sketching
                 ds <- SketchData(x,
                                  assay = "RNA",
                                  ncells = ncells,
                                  seed = myseed)

                 # Seurat create another assay with few cells,
                 # let's keep only the cells resulting from sketching
                 cells <- colnames(ds[["sketch"]])
                 ds <- subset(ds, cells = cells)

                 # remove sketch assay
                 DefaultAssay(ds) <- "RNA"
                 ds[["sketch"]] <- NULL
                 return(ds)
               })

lapply(slim, dim)

# keep only hvg genes and genes found in scGate models
models <- scGate::get_scGateDB()
gate_genes <- lapply(models$human$HiTME, scGate:::table.to.model) %>%
  unlist() %>%
  unname() %>%
  gsub("-", "", .) %>%
  unique()
models_sex <- list(Male = models$human$generic$Male,
                   Female = models$human$generic$Female)
sex_genes <- lapply(models_sex, scGate:::table.to.model) %>%
  unlist() %>%
  unname() %>%
  gsub("-", "", .) %>%
  unique()

# get genes from references

ref.maps <- ProjecTILs::get.reference.maps(collection = "human",
                                           as.list = F)

ref_maps_genes <- lapply(ref.maps, function(x){
  rownames(x@assays$integrated@data)
}) %>%
  unlist() %>%
  unique()

slim.genes <- lapply(slim, function(x){

  # Ensure we work with RNA assay
  DefaultAssay(x) <- "RNA"

  # 1. Find HVGs using only counts
  x <- FindVariableFeatures(
    x,
    selection.method = "vst",
    nfeatures = 300,
    verbose = FALSE
  )

  hvg <- VariableFeatures(x)
  genes <- c(hvg, gate_genes, sex_genes, ref_maps_genes) %>%
    unique()

  # 2. Subset to HVGs only
  x <- subset(x, features = genes)

  # 3. Keep only raw counts assay
  counts <- GetAssayData(x, layer = "counts")

  new <- CreateSeuratObject(counts = counts,
                            meta.data = x@meta.data)

  # 4. Keep only sample_id metadata
  if(sample_id %in% colnames(x@meta.data)){
    # keep cellontology to test layer2
    new@meta.data <- new@meta.data %>%
      select(sample_id = sample)
  } else {
    stop("sample_id column not found in metadata")
  }

  return(new)
})

lapply(slim.genes, dim)

# merge object
seu.final <- merge(slim.genes[[1]], slim.genes[-1])
seu.final <- JoinLayers(seu.final)
dim(seu.final)
format(object.size(seu.final), units = "MB")

# save
saveRDS(seu.final, "~/Documents/Projects/HiTME/tests/testthat/lite_Seurat_test.rds")
