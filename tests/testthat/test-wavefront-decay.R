test_that("wavefront returns correct structure", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  result <- wavefront(species, coords, n_seeds = 3, n_steps = 10,
                      progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc_wavefront")
  expect_equal(result$n_seeds, 3)
  expect_true(!is.null(result$n_steps))
})


test_that("wavefront print works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  result <- wavefront(species, coords, n_seeds = 3, n_steps = 8,
                      progress = FALSE, seed = 1)

  expect_output(print(result), "spacc")
})


test_that("distanceDecay returns correct structure", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  result <- distanceDecay(species, coords, n_seeds = 3,
                          progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc_decay")
  expect_equal(result$n_seeds, 3)
})


test_that("distanceDecay print works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  result <- distanceDecay(species, coords, n_seeds = 3,
                          progress = FALSE, seed = 1)

  expect_output(print(result), "spacc|decay|distance")
})
