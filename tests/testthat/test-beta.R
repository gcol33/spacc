test_that("spaccBeta returns correct structure", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  result <- spaccBeta(species, coords, n_seeds = 5,
                      parallel = FALSE, progress = FALSE)

  expect_s3_class(result, "spacc_beta")
  expect_equal(result$n_sites, 20)
  expect_equal(result$n_seeds, 5)
  expect_equal(result$index, "sorensen")

  # Beta matrices should be n_seeds x (n_sites - 1)
  expect_equal(dim(result$beta_total), c(5, 19))
  expect_equal(dim(result$beta_turnover), c(5, 19))
  expect_equal(dim(result$beta_nestedness), c(5, 19))
})


test_that("Beta components sum to total", {
  skip_on_cran()

  set.seed(123)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  result <- spaccBeta(species, coords, n_seeds = 3,
                      parallel = FALSE, progress = FALSE)

  # turnover + nestedness should equal total (within tolerance)
  computed_total <- result$beta_turnover + result$beta_nestedness
  expect_equal(computed_total, result$beta_total, tolerance = 1e-10)
})


test_that("Beta diversity values are bounded 0-1", {
  skip_on_cran()

  set.seed(456)
  coords <- data.frame(x = runif(25), y = runif(25))
  species <- matrix(rbinom(25 * 12, 1, 0.35), nrow = 25)

  result <- spaccBeta(species, coords, n_seeds = 5,
                      parallel = FALSE, progress = FALSE)

  expect_true(all(result$beta_total >= 0 & result$beta_total <= 1))
  expect_true(all(result$beta_turnover >= 0 & result$beta_turnover <= 1))
  expect_true(all(result$beta_nestedness >= 0 & result$beta_nestedness <= 1))
})


test_that("Jaccard and Sorensen give different results", {
  skip_on_cran()

  set.seed(789)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  result_sor <- spaccBeta(species, coords, n_seeds = 3, index = "sorensen",
                          parallel = FALSE, progress = FALSE, seed = 1)
  result_jac <- spaccBeta(species, coords, n_seeds = 3, index = "jaccard",
                          parallel = FALSE, progress = FALSE, seed = 1)

  # Should be different
  expect_false(identical(result_sor$beta_total, result_jac$beta_total))
})
