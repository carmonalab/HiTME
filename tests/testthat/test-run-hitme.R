# Load test data
test_data_path <- "lite_Seurat_test.rds"

test_that("Run.HiTME validates input objects", {
  # Test that NULL object raises error
  expect_error(
    Run.HiTME(object = NULL),
    "Please provide a Seurat object or a list of them"
  )
})

test_that("Run.HiTME validates input object class", {
  # Test that non-Seurat object raises error
  not_seurat <- list(data = data.frame(a = 1:5, b = 6:10))

  expect_error(
    Run.HiTME(object = not_seurat),
    "not Seurat objects"
  )
})

test_that("Run.HiTME runs with minimal parameters on real dataset", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  # Load lite test dataset
  obj <- readRDS(test_data_path)

  # Test with minimal parameters - no scGate, no ref.maps, no layer3
  result <- Run.HiTME(
    object = obj[, 1:50],
    scGate.model = NULL,
    ref.maps = NULL,
    additional.signatures = NULL,
    layer3 = NULL,
    ncores = 1,
    verbose = FALSE
  )

  # Check that it returns a Seurat object
  expect_s4_class(result, "Seurat")
  # Check that layer1 column was added
  expect_true("layer1" %in% names(result@meta.data))
})

test_that("Run.HiTME runs layer3", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  # Load lite test dataset
  obj <- readRDS(test_data_path)

  # Normalize data for scGate processing
  obj <- Seurat::NormalizeData(obj)

  # Test with minimal parameters - no scGate, no ref.maps, no layer3
  result <- Run.HiTME(
    object = obj[, 1:50],
    scGate.model = NULL,
    ref.maps = NULL,
    ncores = 1,
    verbose = FALSE
  )

  # Check that it returns a Seurat object
  expect_s4_class(result, "Seurat")
  # Check that layer1 column was added
  expect_true("layer3_annotation" %in% names(result@meta.data))
})

test_that("Run.HiTME handles scGate model parameter", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)

  # Normalize data for scGate processing
  obj <- Seurat::NormalizeData(obj)

  # Test with default scGate model - may require external downloads
  result <- Run.HiTME(
    object = obj,
    scGate.model = "default",
    ref.maps = NULL,
    layer3 = NULL,
    ncores = 1,
    species = "human",
    verbose = FALSE
  )

  # Either succeeds or gracefully fails (returns NULL)
  expect_true(is(result, "Seurat"))
})

test_that("Run.HiTME handles split.by parameter", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)

  # Add test metadata column that exists
  obj$test_sample <- sample(c("S1", "S2"), ncol(obj), replace = TRUE)

  # Test split.by with valid column
  result <- Run.HiTME(
    object = obj[, 1:50],
    split.by = "test_sample",
    remerge = TRUE,
    scGate.model = NULL,
    ref.maps = NULL,
    layer3 = NULL,
    ncores = 1,
    verbose = FALSE,
    progressbar = FALSE
  )

  expect_s4_class(result, "Seurat")
})

test_that("Run.HiTME rejects split.by with list input", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)
  obj_list <- list(obj1 = obj[, 1:25], obj2 = obj[, 26:50])

  # Test that split.by with list raises error
  expect_error(
    Run.HiTME(object = obj_list, split.by = "orig.ident"),
    "split.by only supported for a single Seurat object"
  )
})

test_that("Run.HiTME validates species parameter", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)

  # Test that invalid species raises error
  expect_error(
    Run.HiTME(object = obj[, 1:50], species = "zebrafish"),
    "supported species"
  )
})

test_that("Run.HiTME accepts valid species variants", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)

  # Test human variants
  result <- Run.HiTME(
    object = obj[, 1:50],
    species = "human",
    scGate.model = NULL,
    ref.maps = NULL,
    layer3 = NULL,
    ncores = 1,
    verbose = FALSE
  )

  expect_s4_class(result, "Seurat")

  # Test that "homo" variant works
  result2 <- Run.HiTME(
    object = obj[, 1:50],
    species = "homo",
    scGate.model = NULL,
    ref.maps = NULL,
    layer3 = NULL,
    ncores = 1,
    verbose = FALSE
  )

  expect_s4_class(result2, "Seurat")
})

test_that("Run.HiTME accepts NULL species with warning", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)

  # Test that NULL species generates warning
  expect_warning(
    Run.HiTME(
      object = obj[, 1:50],
      species = NULL,
      scGate.model = NULL,
      ref.maps = NULL,
      layer3 = NULL,
      ncores = 1,
      verbose = FALSE
    ),
    "Not using default scGate models"
  )
})

test_that("Run.HiTME handles ncores parameter", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)
  # Normalize data for scGate processing
  obj <- Seurat::NormalizeData(obj)

  # Test with ncores = 1
  result <- Run.HiTME(
    object = obj[, 1:50],
    ncores = 1,
    scGate.model = "default",
    ref.maps = NULL,
    layer3 = NULL,
    verbose = FALSE
  )

  expect_s4_class(result, "Seurat")
})

test_that("Run.HiTME handles ncores > 1", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)
  # Normalize data for scGate processing
  obj <- Seurat::NormalizeData(obj)

  # Test with ncores = 2
  result <- Run.HiTME(
    object = obj[, 1:50],
    ncores = 2,
    scGate.model = "default",
    ref.maps = NULL,
    layer3 = NULL,
    verbose = FALSE
  )

  expect_s4_class(result, "Seurat")
})

test_that("Run.HiTME accepts bparam SnowParam overriding ncores", {
  skip_if_not(file.exists(test_data_path), "Test data not available")
  skip_on_os("windows")

  obj <- readRDS(test_data_path)
  obj <- Seurat::NormalizeData(obj)

  snow_param <- BiocParallel::SnowParam(workers = 2, progressbar = FALSE)

  result <- Run.HiTME(
    object = obj[, 1:40],
    ncores = 1, # should be ignored in favor of bparam
    bparam = snow_param,
    scGate.model = "default",
    species = "human",
    ref.maps = NULL,
    layer3 = NULL,
    verbose = FALSE
  )

  expect_s4_class(result, "Seurat")
})

test_that("Run.HiTME accepts provided BiocParallel params for layer3 scoring", {
  skip_if_not(file.exists(test_data_path), "Test data not available")
  skip_on_os("windows")

  obj <- readRDS(test_data_path)
  obj <- Seurat::NormalizeData(obj)

  # Use a simple SnowParam to ensure BPPARAM is honored when scGate is disabled
  snow_param <- BiocParallel::SnowParam(workers = 2, progressbar = FALSE)

  result <- Run.HiTME(
    object = obj[, 1:40],
    ncores = 1,
    bparam = snow_param,
    scGate.model = NULL,
    layer3 = list(Sig3_test = rownames(obj)[1:10]),
    species = "human",
    ref.maps = NULL,
    verbose = FALSE
  )

  expect_s4_class(result, "Seurat")
  expect_true("layer3" %in% names(result@meta.data))
})

test_that("Run.HiTME respects remerge parameter with lists", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)
  # Normalize data for scGate processing
  obj <- Seurat::NormalizeData(obj)
  obj_list <- list(obj1 = obj[, 1:25], obj2 = obj[, 26:50])

  # Test remerge = TRUE returns single object
  result_merged <- Run.HiTME(
    object = obj_list,
    remerge = TRUE,
    scGate.model = "default",
    ref.maps = NULL,
    layer3 = NULL,
    ncores = 1,
    verbose = FALSE
  )

  expect_s4_class(result_merged, "Seurat")

  # Test remerge = FALSE returns list
  result_list <- Run.HiTME(
    object = obj_list,
    remerge = FALSE,
    scGate.model = "default",
    ref.maps = NULL,
    layer3 = NULL,
    ncores = 1,
    verbose = FALSE
  )

  expect_type(result_list, "list")
  expect_length(result_list, 2)
})

test_that("Run.HiTME handles scGate.model.branch parameter", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)

  # Normalize data for scGate processing
  obj <- Seurat::NormalizeData(obj)

  # Test with master branch - may require external downloads
  result <-
    Run.HiTME(
      object = obj[, 1:50],
      scGate.model = "default",
      scGate.model.branch = "master",
      species = "human",
      ref.maps = NULL,
      layer3 = NULL,
      ncores = 1,
      verbose = FALSE
    )

  expect_true(is(result, "Seurat"))

  result2 <-
    Run.HiTME(
      object = obj[, 1:50],
      scGate.model = "default",
      scGate.model.branch = "dev",
      species = "human",
      ref.maps = NULL,
      layer3 = NULL,
      ncores = 1,
      verbose = FALSE
    )

  expect_true(is(result2, "Seurat"))
})

test_that("Run.HiTME handles multi.asNA parameter", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)

  # Normalize data for scGate processing
  obj <- Seurat::NormalizeData(obj)

  # Test multi.asNA = TRUE
  result_na <-
    Run.HiTME(
      object = obj[, 1:50],
      scGate.model = "default",
      multi.asNA = TRUE,
      species = "human",
      ref.maps = NULL,
      layer3 = NULL,
      ncores = 1,
      verbose = FALSE
    )

  expect_true(is.null(result_na) || is(result_na, "Seurat"))

  # Test multi.asNA = FALSE
  result_multi <-
    Run.HiTME(
      object = obj[, 1:50],
      scGate.model = "default",
      multi.asNA = FALSE,
      species = "human",
      ref.maps = NULL,
      layer3 = NULL,
      ncores = 1,
      verbose = FALSE
    )

  expect_true(is.null(result_multi) || is(result_multi, "Seurat"))
})


test_that("Run.HiTME run layer2", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  # get ref.maps
  ref.maps <- ProjecTILs::get.reference.maps(collection = "human",
                                             as.list = F)

  obj <- readRDS(test_data_path)
  # Normalize data for scGate processing
  obj <- Seurat::NormalizeData(obj)

  result <- Run.HiTME(
    object = obj,
    ncores = 1,
    scGate.model = "default",
    ref.maps = ref.maps,
    layer3 = NULL,
    verbose = FALSE
  )

  expect_s4_class(result, "Seurat")
})

test_that("Run.HiTME merges functional.cluster into layer2", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  ref.maps <- ProjecTILs::get.reference.maps(collection = "human",
                                             as.list = FALSE)

  obj <- readRDS(test_data_path)
  obj <- Seurat::NormalizeData(obj)

  result <- Run.HiTME(
    object = obj,
    ncores = 1,
    scGate.model = "default",
    ref.maps = ref.maps,
    layer3 = NULL,
    verbose = FALSE
  )

  expect_s4_class(result, "Seurat")
  expect_true(all(c("functional.cluster", "functional.cluster.conf", "layer2") %in%
                    names(result@meta.data)))
  expected_layer2 <- dplyr::if_else(
    is.na(result$functional.cluster),
    as.character(result$layer1),
    as.character(result$functional.cluster)
  )
  expect_equal(as.character(result$layer2), expected_layer2)
})

test_that("Run.HiTME records layer2 metadata levels", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  ref.maps <- ProjecTILs::get.reference.maps(collection = "human",
                                             as.list = FALSE)

  obj <- readRDS(test_data_path)
  obj <- Seurat::NormalizeData(obj)

  result <- Run.HiTME(
    object = obj,
    ncores = 1,
    scGate.model = "default",
    ref.maps = ref.maps,
    layer3 = NULL,
    verbose = FALSE
  )

  ref_levels <- vapply(ref.maps, function(x) x@misc$layer1_link, character(1))
  misc_layer2 <- result@misc[["layer2_param"]][["functional.cluster"]]

  expect_true(is.list(misc_layer2))
  expect_equal(misc_layer2[["References_user_specified"]], names(ref.maps))
  expect_setequal(names(misc_layer2[["levels2_per_levels1"]]), ref_levels)
  expect_true(all(lengths(misc_layer2[["levels2_per_levels1"]]) > 0))
})

test_that("Run.HiTME falls back to layer1 when no layer2 mapping is applicable", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  ref.maps <- ProjecTILs::get.reference.maps(collection = "human",
                                             as.list = FALSE)
  ref.maps <- lapply(ref.maps, function(r) {
    r@misc$layer1_link <- "nonmatching_layer1"
    r
  })

  obj <- readRDS(test_data_path)
  obj <- Seurat::NormalizeData(obj)

  result <- Run.HiTME(
    object = obj,
    ncores = 1,
    scGate.model = "default",
    ref.maps = ref.maps,
    layer1_link = "layer1",
    layer3 = NULL,
    verbose = FALSE
  )

  expect_s4_class(result, "Seurat")
  expect_true(all(is.na(result$functional.cluster)))
  expect_equal(as.character(result$layer2), as.character(result$layer1))
})
