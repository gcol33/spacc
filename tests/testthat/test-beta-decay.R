test_that("betaDecay returns correct structure", {
  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 15, 1, 0.3), nrow = 20)

  result <- betaDecay(species, coords, progress = FALSE)

  expect_s3_class(result, "spacc_beta_decay")
  expect_true("pairs" %in% names(result))
  expect_true("fits" %in% names(result))
  expect_true("best_model" %in% names(result))
  expect_true("coefficients" %in% names(result))
})


test_that("betaDecay pair count is n*(n-1)/2", {
  set.seed(42)
  n <- 15
  coords <- data.frame(x = runif(n), y = runif(n))
  species <- matrix(rbinom(n * 10, 1, 0.3), nrow = n)

  result <- betaDecay(species, coords, progress = FALSE)

  expect_equal(nrow(result$pairs), n * (n - 1) / 2)
  expect_equal(result$n_pairs, n * (n - 1) / 2)
})


test_that("dissimilarity values in [0, 1]", {
  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 15, 1, 0.3), nrow = 20)

  result <- betaDecay(species, coords, progress = FALSE)

  expect_true(all(result$pairs$dissimilarity >= 0))
  expect_true(all(result$pairs$dissimilarity <= 1))
})


test_that("sorensen and jaccard both work", {
  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 10, 1, 0.4), nrow = 15)

  sor <- betaDecay(species, coords, index = "sorensen", progress = FALSE)
  jac <- betaDecay(species, coords, index = "jaccard", progress = FALSE)

  expect_equal(sor$index, "sorensen")
  expect_equal(jac$index, "jaccard")
  # Jaccard >= Sorensen for same pair
  expect_true(mean(jac$pairs$dissimilarity) >= mean(sor$pairs$dissimilarity) - 1e-10)
})


test_that("at least one model fits", {
  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 15, 1, 0.3), nrow = 20)

  result <- betaDecay(species, coords, model = "all", progress = FALSE)

  expect_true(nrow(result$coefficients) > 0)
  expect_true(!is.na(result$best_model))
})


test_that("single model fitting works", {
  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 15, 1, 0.3), nrow = 20)

  result <- betaDecay(species, coords, model = "linear", progress = FALSE)

  expect_equal(nrow(result$coefficients), 1)
  expect_equal(result$coefficients$model, "linear")
})


test_that("spacc_dist input works", {
  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 10, 1, 0.4), nrow = 15)

  d <- distances(coords)
  result <- betaDecay(species, d, progress = FALSE)

  expect_s3_class(result, "spacc_beta_decay")
})


test_that("print.spacc_beta_decay works", {
  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 10, 1, 0.4), nrow = 15)

  result <- betaDecay(species, coords, progress = FALSE)

  expect_output(print(result), "Beta distance-decay")
  expect_output(print(result), "Fitted models")
})


test_that("coef.spacc_beta_decay returns coefficients", {
  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 15, 1, 0.3), nrow = 20)

  result <- betaDecay(species, coords, progress = FALSE)

  cf <- coef(result)
  expect_true(is.numeric(cf))
  expect_true("a" %in% names(cf))
  expect_true("b" %in% names(cf))
})


test_that("as.data.frame.spacc_beta_decay returns pairs", {
  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 10, 1, 0.4), nrow = 15)

  result <- betaDecay(species, coords, progress = FALSE)
  df <- as.data.frame(result)

  expect_s3_class(df, "data.frame")
  expect_true("distance" %in% names(df))
  expect_true("dissimilarity" %in% names(df))
})


test_that("plot.spacc_beta_decay works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 15, 1, 0.3), nrow = 20)

  result <- betaDecay(species, coords, progress = FALSE)

  p <- plot(result)
  expect_s3_class(p, "ggplot")
})
