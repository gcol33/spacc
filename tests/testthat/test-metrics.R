test_that("spaccMetrics returns correct structure", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  result <- spaccMetrics(species, coords,
                         metrics = c("slope_10", "half_richness", "auc"),
                         parallel = FALSE, progress = FALSE)

  expect_s3_class(result, "spacc_metrics")
  expect_equal(result$n_sites, 15)
  expect_equal(result$n_species, 8)
  expect_equal(result$metric_names, c("slope_10", "half_richness", "auc"))

  # Metrics data frame should have correct dimensions
  expect_equal(nrow(result$metrics), 15)
  expect_true("slope_10" %in% names(result$metrics))
  expect_true("half_richness" %in% names(result$metrics))
  expect_true("auc" %in% names(result$metrics))
  expect_true("x" %in% names(result$metrics))
  expect_true("y" %in% names(result$metrics))
})


test_that("slope metrics are positive", {
  skip_on_cran()

  set.seed(123)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  result <- spaccMetrics(species, coords,
                         metrics = c("slope_10", "slope_25"),
                         parallel = FALSE, progress = FALSE)

  # Slopes should generally be positive (species accumulate)
  expect_true(all(result$metrics$slope_10 >= 0))
  expect_true(all(result$metrics$slope_25 >= 0))
})


test_that("half_richness is bounded", {
  skip_on_cran()

  set.seed(456)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.4), nrow = 20)

  result <- spaccMetrics(species, coords,
                         metrics = c("half_richness"),
                         parallel = FALSE, progress = FALSE)

  # half_richness should be between 1 and n_sites
  expect_true(all(result$metrics$half_richness >= 1))
  expect_true(all(result$metrics$half_richness <= 20))
})


test_that("AUC is positive", {
  skip_on_cran()

  set.seed(789)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.35), nrow = 15)

  result <- spaccMetrics(species, coords,
                         metrics = c("auc"),
                         parallel = FALSE, progress = FALSE)

  # AUC should be positive
  expect_true(all(result$metrics$auc > 0))
})


test_that("different methods give different results", {
  skip_on_cran()

  set.seed(101)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  result_knn <- spaccMetrics(species, coords, metrics = c("auc"),
                              method = "knn", parallel = FALSE, progress = FALSE)
  result_kncn <- spaccMetrics(species, coords, metrics = c("auc"),
                               method = "kncn", parallel = FALSE, progress = FALSE)

  # Methods should produce different results
  expect_false(identical(result_knn$metrics$auc, result_kncn$metrics$auc))
})


test_that("print and summary methods work", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(10), y = runif(10))
  species <- matrix(rbinom(10 * 5, 1, 0.4), nrow = 10)

  result <- spaccMetrics(species, coords,
                         metrics = c("slope_10", "auc"),
                         parallel = FALSE, progress = FALSE)

  expect_output(print(result), "spacc_metrics")
  expect_output(summary(result), "Metric summary")
})


test_that("final_richness metric works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  result <- spaccMetrics(species, coords,
                         metrics = c("final_richness"),
                         parallel = FALSE, progress = FALSE)

  expect_true(all(result$metrics$final_richness > 0))
  expect_true(all(result$metrics$final_richness <= 8))
})


test_that("richness_75pct and richness_90pct work", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  result <- spaccMetrics(species, coords,
                         metrics = c("richness_75pct", "richness_90pct"),
                         parallel = FALSE, progress = FALSE)

  expect_true(all(result$metrics$richness_75pct >= 1))
  expect_true(all(result$metrics$richness_90pct >= result$metrics$richness_75pct))
})


test_that("spaccMetrics with spacc_dist coords", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  d <- distances(coords)
  result <- spaccMetrics(species, d,
                         metrics = c("slope_10", "auc"),
                         parallel = FALSE, progress = FALSE)

  expect_s3_class(result, "spacc_metrics")
})


test_that("spaccMetrics errors on unknown metric", {
  coords <- data.frame(x = 1:5, y = 1:5)
  species <- matrix(rbinom(5 * 3, 1, 0.5), nrow = 5)

  expect_error(
    spaccMetrics(species, coords, metrics = c("nonexistent"),
                 parallel = FALSE, progress = FALSE),
    "Unknown metric"
  )
})
