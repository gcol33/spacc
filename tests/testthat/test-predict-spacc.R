test_that("predict.spacc returns data frame with CI by default", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)

  result <- predict(sac, n = c(10, 20, 30))

  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 3)
  expect_true(all(c("n", "mean", "lower", "upper") %in% names(result)))
  expect_true(all(!is.na(result$mean)))
})


test_that("predict.spacc returns numeric vector with ci = FALSE", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)

  result <- predict(sac, n = c(10, 20, 30), ci = FALSE)

  expect_true(is.numeric(result))
  expect_length(result, 3)
  expect_true(all(!is.na(result)))
})


test_that("predict.spacc uses default n when not provided", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)

  result <- predict(sac)

  expect_true(is.data.frame(result))
  expect_true(nrow(result) >= 2)  # at least 2 unique default values
})


test_that("predict.spacc warns for out-of-range n", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)

  expect_warning(
    result <- predict(sac, n = c(10, 100)),
    "outside observed range"
  )
  expect_true(is.na(result$mean[result$n == 100]))
})


test_that("predict.spacc out-of-range with warn = FALSE suppresses warning", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)

  expect_no_warning(
    result <- predict(sac, n = c(10, 100), warn = FALSE)
  )
  expect_true(is.na(result$mean[result$n == 100]))
})


test_that("predict.spacc values are monotonically non-decreasing", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(50), y = runif(50))
  species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)

  sac <- spacc(species, coords, n_seeds = 10, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)

  result <- predict(sac, n = seq(5, 50, by = 5), ci = FALSE)
  expect_true(all(diff(result) >= -1e-10))  # allow tiny floating point diffs
})


test_that("predict.spacc CI bounds are consistent", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)

  result <- predict(sac, n = c(10, 20, 30))

  expect_true(all(result$lower <= result$mean))
  expect_true(all(result$upper >= result$mean))
})
