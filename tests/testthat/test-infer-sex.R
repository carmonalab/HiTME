# Load test data
test_data_path <- "lite_Seurat_test.rds"

test_that("infer.Sex validates Seurat object input", {
  # Test that NULL object raises error
  expect_error(
    infer.Sex(object = NULL),
    "Please provide a Seurat object, count matrix or a list of them"
  )
})

test_that("infer.Sex rejects non-Seurat/non-matrix objects", {
  # Test with invalid object type
  expect_error(
    infer.Sex(object = "not_a_seurat"),
    "Please provide a Seurat object, count matrix or a list of them"
  )
})

test_that("infer.Sex accepts Seurat object input", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)

  # Test with Seurat object
  result <- infer.Sex(object = obj[, 1:50], infer.level = "sample")

  # Should return Seurat object with sex inference metadata
  expect_s4_class(result, "Seurat")
})

test_that("infer.Sex accepts count matrix input", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)
  counts_matrix <- Seurat::GetAssayData(obj[, 1:50], layer = "counts")

  # Test with count matrix
  result <- infer.Sex(object = counts_matrix, infer.level = "sample")

  expect_true(!is.null(result))
})

test_that("infer.Sex accepts sparse matrix input", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)
  counts_matrix <- methods::as(Seurat::GetAssayData(obj[, 1:50], layer = "counts"), "dgCMatrix")

  # Test with sparse matrix
  result <- infer.Sex(object = counts_matrix, infer.level = "sample")

  expect_true(!is.null(result))
})

test_that("infer.Sex validates infer.level parameter", {
  # Test with invalid infer.level
  expect_error(
    infer.Sex(object = matrix(), infer.level = "invalid_level"),
    "Please check infer.level parameter"
  )
})

test_that("infer.Sex accepts infer.level='sample'", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)
  # Test with infer.level = sample
  result <- infer.Sex(object = obj[, 1:50], infer.level = "sample")

  expect_s4_class(result, "Seurat")
})

test_that("infer.Sex accepts infer.level='cell'", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)
  # Normalize data for scGate processing
  obj <- Seurat::NormalizeData(obj)
  # Test with infer.level = cell - lite dataset may not have all required layers
  # but function should accept the parameter and either succeed or fail gracefully
  result <- infer.Sex(object = obj[, 1:50], infer.level = "cell")

  expect_true(is(result, "Seurat"))
})

test_that("infer.Sex rejects split.by with list input", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)
  obj_list <- Seurat::SplitObject(obj, split.by = "sample_id")

  # Test that split.by with list raises error
  expect_error(
    infer.Sex(object = obj_list, split.by = "sample_id"),
    "split.by only supported for a single object"
  )
})

test_that("infer.Sex rejects split.by with matrix input", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)
  counts_matrix <- Seurat::GetAssayData(obj[, 1:50], layer = "counts")

  # Test that split.by with matrix raises error (no metadata)
  expect_error(
    infer.Sex(object = counts_matrix, split.by = "sample"),
    "When providing matrix or dgCMatrix"
  )
})

test_that("infer.Sex validates split.by metadata column", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)

  # Test with nonexistent split.by column
  expect_error(
    infer.Sex(object = obj[, 1:50], split.by = "nonexistent_column"),
    "not a metadata column"
  )
})

test_that("infer.Sex accepts valid split.by parameter", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)

  # Test with valid split.by using existing sample_id metadata
  result <- infer.Sex(object = obj[, 1:50], split.by = "sample_id")

  expect_s4_class(result, "Seurat")
})

test_that("infer.Sex accepts return.Seurat=TRUE for Seurat input", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)

  # Test with return.Seurat = TRUE
  result <- infer.Sex(object = obj[, 1:50], return.Seurat = TRUE)

  expect_s4_class(result, "Seurat")
})

test_that("infer.Sex accepts return.Seurat=FALSE", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)

  # Test with return.Seurat = FALSE
  result <- infer.Sex(object = obj[, 1:50], return.Seurat = FALSE)

  expect_true(is.data.frame(result))
})

test_that("infer.Sex accepts ncores parameter", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)

  # Test with ncores = 1
  result <- infer.Sex(object = obj[, 1:50], ncores = 1)

  expect_true(!is.null(result))
})

test_that("infer.Sex accepts ncores > 1", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)

  # Test with ncores = 2
  result <- infer.Sex(object = obj[, 1:50], ncores = 2)

  expect_true(!is.null(result))
})

test_that("infer.Sex accepts progressbar parameter", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)

  # Test with progressbar = TRUE
  result <- infer.Sex(object = obj[, 1:50], progressbar = TRUE)

  expect_true(!is.null(result))
})

test_that("infer.Sex accepts verbose parameter", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)

  # Test with verbose = TRUE
  result <- infer.Sex(object = obj[, 1:50], verbose = TRUE)

  expect_true(!is.null(result))
})
