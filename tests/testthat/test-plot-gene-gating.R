# Load test data
test_data_path <- "lite_Seurat_test.rds"

test_that("plot.geneGating validates input object", {
  # Test that NULL object raises error
  expect_error(
    plot.geneGating(object = NULL, scGate.model = "model1"),
    "Please provide a Seurat object or a list of them"
  )
})

test_that("plot.geneGating requires scGate model", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)

  # Test that NULL model raises error
  expect_error(
    plot.geneGating(object = obj[, 1:50], scGate.model = NULL),
    "Please provide a scGate model or list of them"
  )
})

test_that("plot.geneGating rejects list input with split.by", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)
  obj_list <- list(obj1 = obj[, 1:25], obj2 = obj[, 26:50])

  # Create a simple scGate model list
  model <- scGate::gating_model(name = "test_model",
                                signature = c("CD74"))

  # Test that split.by with list raises error
  expect_error(
    plot.geneGating(object = obj_list, scGate.model = model, split.by = "sample_id"),
    "split.by only supported for a single Seurat object"
  )
})

test_that("plot.geneGating validates group.by parameter", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)

  # Create a simple scGate model with available genes
  model <- scGate::gating_model(name = "test_model",
                                signature = c("CD74"))

  # Test with nonexistent group.by column
  expect_error(
    plot.geneGating(object = obj[, 1:50], scGate.model = model, group.by = "nonexistent"),
    "not a metadata column"
  )
})

test_that("plot.geneGating validates split.by parameter", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)

  # Create a scGate model with available genes
  model <- scGate::gating_model(name = "test_model",
                                signature = c("CD74"))

  # Test with nonexistent split.by column
  expect_error(
    plot.geneGating(object = obj[, 1:50], scGate.model = model, group.by = "sample_id", split.by = "nonexistent"),
    "not a metadata column"
  )
})

