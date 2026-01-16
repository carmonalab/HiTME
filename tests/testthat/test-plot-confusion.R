# Load test data
test_data_path <- "lite_Seurat_test.rds"

test_that("plot.confusion validates input object", {
  # Test that NULL object raises error
  expect_error(
    plot.confusion(object = NULL, var.1 = "col1", var.2 = "col2"),
    "Please provide a Seurat object or its metadata"
  )
})

test_that("plot.confusion accepts Seurat object from real data", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)

  # Add test metadata columns
  obj$celltype1 <- sample(c("TypeA", "TypeB"), ncol(obj), replace = TRUE)
  obj$celltype2 <- sample(c("TypeX", "TypeY"), ncol(obj), replace = TRUE)

  # Test that valid Seurat object produces ggplot
  result <- plot.confusion(
    object = obj,
    var.1 = "celltype1",
    var.2 = "celltype2",
    relative = FALSE,
    type = "tile"
  )

  expect_s3_class(result, "ggplot")
})

test_that("plot.confusion accepts data.frame metadata", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)

  # Create a data frame with cell type assignments
  metadata <- data.frame(
    celltype1 = sample(c("TypeA", "TypeB"), ncol(obj), replace = TRUE),
    celltype2 = sample(c("TypeX", "TypeY"), ncol(obj), replace = TRUE)
  )

  # Test that data.frame is accepted
  result <- plot.confusion(
    object = metadata,
    var.1 = "celltype1",
    var.2 = "celltype2"
  )

  expect_s3_class(result, "ggplot")
})

test_that("plot.confusion rejects invalid object types", {
  # Test with invalid object type
  invalid_obj <- list(a = 1, b = 2)

  expect_error(
    plot.confusion(object = invalid_obj, var.1 = "col1", var.2 = "col2"),
    "Please provide a Seurat object or its metadata"
  )
})

test_that("plot.confusion requires var.1 and var.2", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)
  obj$celltype <- sample(c("TypeA", "TypeB"), ncol(obj), replace = TRUE)

  # Test missing var.1 (will raise error when trying to access NULL column)
  expect_error(
    plot.confusion(object = obj, var.1 = NULL, var.2 = "celltype"),
    class = "error"
  )

  # Test missing var.2 (will raise error when trying to access NULL column)
  expect_error(
    plot.confusion(object = obj, var.1 = "celltype", var.2 = NULL),
    class = "error"
  )
})

test_that("plot.confusion validates metadata column existence", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)
  obj$celltype1 <- sample(c("TypeA", "TypeB"), ncol(obj), replace = TRUE)

  # Test with nonexistent column
  expect_error(
    plot.confusion(object = obj, var.1 = "celltype1", var.2 = "nonexistent"),
    "Not all classification variables"
  )
})

test_that("plot.confusion handles relative parameter", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)
  obj$celltype1 <- sample(c("TypeA", "TypeB"), ncol(obj), replace = TRUE)
  obj$celltype2 <- sample(c("TypeX", "TypeY"), ncol(obj), replace = TRUE)

  # Test with relative = FALSE (default)
  p1 <- plot.confusion(
    object = obj,
    var.1 = "celltype1",
    var.2 = "celltype2",
    relative = FALSE
  )
  expect_s3_class(p1, "ggplot")

  # Test with relative = TRUE
  p2 <- plot.confusion(
    object = obj,
    var.1 = "celltype1",
    var.2 = "celltype2",
    relative = TRUE
  )
  expect_s3_class(p2, "ggplot")
})

test_that("plot.confusion handles plot type options", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)
  obj$celltype1 <- sample(c("TypeA", "TypeB"), ncol(obj), replace = TRUE)
  obj$celltype2 <- sample(c("TypeX", "TypeY"), ncol(obj), replace = TRUE)

  # Test with tile plot (default)
  p_tile <- plot.confusion(
    object = obj,
    var.1 = "celltype1",
    var.2 = "celltype2",
    type = "tile"
  )
  expect_s3_class(p_tile, "ggplot")

  # Test with bar plot
  p_bar <- plot.confusion(
    object = obj,
    var.1 = "celltype1",
    var.2 = "celltype2",
    type = "bar"
  )
  expect_s3_class(p_bar, "ggplot")
})

test_that("plot.confusion handles custom labels", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)
  obj$celltype1 <- sample(c("TypeA", "TypeB"), ncol(obj), replace = TRUE)
  obj$celltype2 <- sample(c("TypeX", "TypeY"), ncol(obj), replace = TRUE)

  # Test with custom labels
  p <- plot.confusion(
    object = obj,
    var.1 = "celltype1",
    var.2 = "celltype2",
    xlab = "Classification Method 1",
    ylab = "Classification Method 2",
    plot.title = "Cell Type Confusion Matrix"
  )
  expect_s3_class(p, "ggplot")
})

test_that("plot.confusion handles NA values", {
  skip_if_not(file.exists(test_data_path), "Test data not available")

  obj <- readRDS(test_data_path)

  # Create metadata with NA values
  obj$celltype1 <- c("TypeA", "TypeB", NA, rep(c("TypeA", "TypeB"), length.out = ncol(obj) - 3))
  obj$celltype2 <- c("TypeX", "TypeY", "TypeX", rep(c("TypeX", "TypeY"), length.out = ncol(obj) - 3))

  # Test with useNA = "ifany" (default)
  p1 <- plot.confusion(
    object = obj,
    var.1 = "celltype1",
    var.2 = "celltype2",
    useNA = "ifany"
  )
  expect_s3_class(p1, "ggplot")

  # Test with useNA = "no"
  p2 <- plot.confusion(
    object = obj,
    var.1 = "celltype1",
    var.2 = "celltype2",
    useNA = "no"
  )
  expect_s3_class(p2, "ggplot")
})
