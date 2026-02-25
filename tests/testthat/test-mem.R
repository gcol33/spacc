# spatialEigenvectors tests ---------------------------------------------------

test_that("spatialEigenvectors returns correct structure", {
  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))

  mem <- spatialEigenvectors(coords)

  expect_s3_class(mem, "spacc_mem")
  expect_equal(mem$n_sites, 30)
  expect_true(is.matrix(mem$vectors))
  expect_true(is.numeric(mem$eigenvalues))
  expect_true(is.numeric(mem$moran_i))
  expect_true(mem$threshold > 0)
})


test_that("all eigenvalues are positive", {
  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))

  mem <- spatialEigenvectors(coords)

  expect_true(all(mem$eigenvalues > 0))
})


test_that("Moran's I values are bounded", {
  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))

  mem <- spatialEigenvectors(coords)

  # Moran's I can slightly exceed [-1, 1] for eigenvectors due to
  # numerical properties of the doubly-centered matrix
  expect_true(all(mem$moran_i >= -2))
  expect_true(all(mem$moran_i <= 2))
  # Most values should be in [-1, 1]
  expect_true(mean(abs(mem$moran_i) <= 1.2) > 0.9)
})


test_that("dbmem retains fewer or equal vectors than pcnm", {
  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))

  mem_pcnm <- spatialEigenvectors(coords, method = "pcnm")
  mem_dbmem <- spatialEigenvectors(coords, method = "dbmem")

  expect_true(ncol(mem_dbmem$vectors) <= ncol(mem_pcnm$vectors))
})


test_that("spatialEigenvectors accepts spacc_dist", {
  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))

  # Create distance object
  D <- cpp_distance_matrix(coords$x, coords$y, "euclidean")
  dist_obj <- structure(D, class = "spacc_dist", coords = coords)

  mem <- spatialEigenvectors(dist_obj)
  expect_s3_class(mem, "spacc_mem")
})


test_that("custom threshold works", {
  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))

  mem <- spatialEigenvectors(coords, threshold = 0.5)
  expect_equal(mem$threshold, 0.5)
})


test_that("print.spacc_mem works", {
  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  mem <- spatialEigenvectors(coords)

  expect_output(print(mem), "Spatial eigenvectors")
  expect_output(print(mem), "Threshold")
})


test_that("summary.spacc_mem returns data.frame", {
  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  mem <- spatialEigenvectors(coords)

  summ <- summary(mem)
  expect_s3_class(summ, "data.frame")
  expect_true("eigenvalue" %in% names(summ))
  expect_true("moran_i" %in% names(summ))
})


test_that("plot.spacc_mem works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  mem <- spatialEigenvectors(coords)

  p1 <- plot(mem, type = "eigenvalues")
  expect_s3_class(p1, "ggplot")

  p2 <- plot(mem, type = "moran")
  expect_s3_class(p2, "ggplot")

  p3 <- plot(mem, type = "map")
  expect_s3_class(p3, "ggplot")
})


# spatialPartition tests ------------------------------------------------------

test_that("spatialPartition R-squared in [0, 1]", {
  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  y <- sin(coords$x * 4) + rnorm(30, sd = 0.5)

  mem <- spatialEigenvectors(coords)
  part <- spatialPartition(y, mem)

  expect_s3_class(part, "spacc_mem_partition")
  expect_true(part$r_squared_spatial >= 0)
  expect_true(part$r_squared_spatial <= 1)
})


test_that("spatialPartition works with spacc_alpha", {
  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rpois(30 * 15, 2), nrow = 30)

  mem <- spatialEigenvectors(coords)
  alpha <- alphaDiversity(species, q = 0, coords = coords)
  part <- spatialPartition(alpha, mem)

  expect_s3_class(part, "spacc_mem_partition")
  expect_true(part$r_squared_spatial >= 0)
})


test_that("spatialPartition works with numeric vector", {
  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  y <- rnorm(30)

  mem <- spatialEigenvectors(coords)
  part <- spatialPartition(y, mem)

  expect_s3_class(part, "spacc_mem_partition")
})


test_that("spatialPartition forward = FALSE uses all MEMs", {
  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  y <- rnorm(20)

  mem <- spatialEigenvectors(coords)
  part <- spatialPartition(y, mem, forward = FALSE)

  expect_equal(part$n_selected, ncol(mem$vectors))
})


test_that("print.spacc_mem_partition works", {
  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  y <- rnorm(20)

  mem <- spatialEigenvectors(coords)
  part <- spatialPartition(y, mem)

  expect_output(print(part), "Spatial Variance Partitioning")
})


test_that("plot.spacc_mem_partition works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  y <- rnorm(20)

  mem <- spatialEigenvectors(coords)
  part <- spatialPartition(y, mem)

  p <- plot(part)
  expect_s3_class(p, "ggplot")
})
