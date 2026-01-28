test_that("dar returns correct structure", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 2), nrow = 20)

  result <- dar(species, coords, q = c(0, 1), n_seeds = 3,
                area_method = "count",
                parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc_dar")
  expect_equal(result$n_sites, 20)
  expect_equal(result$n_seeds, 3)
  expect_equal(result$q, c(0, 1))
  expect_equal(result$area_method, "count")
  expect_equal(dim(result$area), c(3, 20))
})


test_that("dar print works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 2), nrow = 15)

  result <- dar(species, coords, q = 0, n_seeds = 3,
                area_method = "count",
                parallel = FALSE, progress = FALSE, seed = 1)

  expect_output(print(result), "spacc DAR")
})


test_that("dar summary returns data.frame", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 2), nrow = 15)

  result <- dar(species, coords, q = c(0, 1), n_seeds = 3,
                area_method = "count",
                parallel = FALSE, progress = FALSE, seed = 1)

  summ <- summary(result)
  expect_s3_class(summ, "data.frame")
  expect_true("q" %in% names(summ))
  expect_true("mean_diversity" %in% names(summ))
})


test_that("sesars returns correct structure", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.4), nrow = 20)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)
  effort <- rpois(20, 10) + 1

  result <- sesars(sac, effort, model = "power")

  expect_s3_class(result, "spacc_sesars")
  expect_equal(result$model, "power")
  expect_true(!is.na(result$r_squared))
  expect_true(length(result$coef) == 3)
})


test_that("sesars additive model works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.4), nrow = 20)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)
  effort <- rpois(20, 10) + 1

  result <- sesars(sac, effort, model = "additive")

  expect_s3_class(result, "spacc_sesars")
  expect_equal(result$model, "additive")
})


test_that("sesars print works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.4), nrow = 20)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)
  effort <- rpois(20, 10) + 1

  result <- sesars(sac, effort)
  expect_output(print(result), "SESARS")
})


test_that("sfar returns correct structure", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.4), nrow = 20)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)
  patches <- rep(1:4, each = 5)

  result <- sfar(sac, patches)

  expect_s3_class(result, "spacc_sfar")
  expect_equal(result$n_patches, 4)
  expect_true(!is.na(result$r_squared))
  expect_true("z" %in% names(result$coef))
  expect_true("f" %in% names(result$coef))
})


test_that("sfar print works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.4), nrow = 20)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)
  patches <- rep(1:4, each = 5)

  result <- sfar(sac, patches)
  expect_output(print(result), "SFAR")
})


test_that("spaccEndemism returns correct structure", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  result <- spaccEndemism(species, coords, n_seeds = 3,
                           parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc_endemism")
  expect_equal(result$n_sites, 20)
  expect_equal(result$n_seeds, 3)
  expect_equal(dim(result$richness), c(3, 20))
  expect_equal(dim(result$endemism), c(3, 20))
})


test_that("spaccEndemism final endemism equals total richness", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  result <- spaccEndemism(species, coords, n_seeds = 3,
                           parallel = FALSE, progress = FALSE, seed = 1)

  # At the last step (all sites), all species are "endemic" (= present in accumulated = total)
  for (s in 1:3) {
    expect_equal(result$endemism[s, 15], result$richness[s, 15])
  }
})


test_that("spaccEndemism endemism <= richness", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  result <- spaccEndemism(species, coords, n_seeds = 3,
                           parallel = FALSE, progress = FALSE, seed = 1)

  expect_true(all(result$endemism <= result$richness))
})


test_that("spaccEndemism print works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  result <- spaccEndemism(species, coords, n_seeds = 3,
                           parallel = FALSE, progress = FALSE, seed = 1)

  expect_output(print(result), "spacc endemism")
})
