test_that("autoplot.spacc works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)
  sac <- spacc(species, coords, n_seeds = 3, progress = FALSE)

  p <- ggplot2::autoplot(sac)
  expect_s3_class(p, "ggplot")
})


test_that("autoplot.spacc_hill works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 3), nrow = 15)
  hill <- spaccHill(species, coords, q = c(0, 1), n_seeds = 3,
                     parallel = FALSE, progress = FALSE)

  p <- ggplot2::autoplot(hill)
  expect_s3_class(p, "ggplot")
})


test_that("autoplot.spacc_beta works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)
  beta <- spaccBeta(species, coords, n_seeds = 3,
                     parallel = FALSE, progress = FALSE)

  p <- ggplot2::autoplot(beta)
  expect_s3_class(p, "ggplot")
})


test_that("autoplot.spacc_rare works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)
  rare <- rarefy(species, n_boot = 10)

  p <- ggplot2::autoplot(rare)
  expect_s3_class(p, "ggplot")
})


test_that("autoplot.spacc_completeness works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  species <- matrix(rpois(50 * 30, 2), nrow = 50)
  comp <- completenessProfile(species, q = c(0, 1, 2))

  p <- ggplot2::autoplot(comp)
  expect_s3_class(p, "ggplot")
})


test_that("autoplot.spacc_zeta works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 15, 1, 0.4), nrow = 20)
  zeta <- zetaDiversity(species, coords, orders = 1:3,
                         n_samples = 20, progress = FALSE)

  p <- ggplot2::autoplot(zeta)
  expect_s3_class(p, "ggplot")
})
