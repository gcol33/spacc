test_that("zetaDiversity returns correct structure", {
  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 15, 1, 0.4), nrow = 20)

  zeta <- zetaDiversity(species, coords, orders = 1:5,
                         n_samples = 20, progress = FALSE)

  expect_s3_class(zeta, "spacc_zeta")
  expect_equal(length(zeta$zeta), 5)
  expect_equal(length(zeta$zeta_sd), 5)
  expect_equal(zeta$orders, 1:5)
  expect_equal(zeta$n_sites, 20)
  expect_equal(zeta$n_species, 15)
})


test_that("zeta values are monotonically non-increasing", {
  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 20, 1, 0.4), nrow = 30)

  zeta <- zetaDiversity(species, coords, orders = 1:5,
                         n_samples = 50, progress = FALSE)

  # Zeta should generally decrease with order
  # Allow small tolerance for sampling noise
  for (i in 2:length(zeta$zeta)) {
    expect_true(zeta$zeta[i] <= zeta$zeta[i - 1] + 1,
                info = sprintf("zeta_%d > zeta_%d", zeta$orders[i], zeta$orders[i - 1]))
  }
})


test_that("zeta_1 equals mean site richness", {
  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 15, 1, 0.4), nrow = 20)

  zeta <- zetaDiversity(species, coords, orders = 1:3,
                         n_samples = 20, progress = FALSE)

  expected_richness <- mean(rowSums((species > 0) * 1L))
  expect_equal(zeta$zeta[1], expected_richness)
})


test_that("random method works", {
  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 15, 1, 0.4), nrow = 20)

  zeta <- zetaDiversity(species, coords, orders = 1:3,
                         method = "random", n_samples = 20, progress = FALSE)

  expect_s3_class(zeta, "spacc_zeta")
  expect_equal(zeta$method, "random")
})


test_that("seed gives reproducible results", {
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 15, 1, 0.4), nrow = 20)

  z1 <- zetaDiversity(species, coords, orders = 1:3,
                       n_samples = 20, seed = 123, progress = FALSE)
  z2 <- zetaDiversity(species, coords, orders = 1:3,
                       n_samples = 20, seed = 123, progress = FALSE)

  expect_equal(z1$zeta, z2$zeta)
})


test_that("decline model fits are computed", {
  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 20, 1, 0.4), nrow = 30)

  zeta <- zetaDiversity(species, coords, orders = 1:8,
                         n_samples = 30, progress = FALSE)

  expect_s3_class(zeta$decline, "data.frame")
  expect_equal(nrow(zeta$decline), 2)
  expect_true("exponential" %in% zeta$decline$model)
  expect_true("power_law" %in% zeta$decline$model)
})


test_that("print works", {
  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 15, 1, 0.4), nrow = 20)

  zeta <- zetaDiversity(species, coords, orders = 1:3,
                         n_samples = 20, progress = FALSE)

  expect_output(print(zeta), "Zeta Diversity")
})


test_that("summary returns data.frame", {
  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 15, 1, 0.4), nrow = 20)

  zeta <- zetaDiversity(species, coords, orders = 1:3,
                         n_samples = 20, progress = FALSE)

  summ <- summary(zeta)
  expect_s3_class(summ, "data.frame")
  expect_true("zeta" %in% names(summ))
})


test_that("as.data.frame works", {
  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 15, 1, 0.4), nrow = 20)

  zeta <- zetaDiversity(species, coords, orders = 1:3,
                         n_samples = 20, progress = FALSE)

  df <- as.data.frame(zeta)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 3)
})


test_that("plot decline works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 15, 1, 0.4), nrow = 20)

  zeta <- zetaDiversity(species, coords, orders = 1:5,
                         n_samples = 20, progress = FALSE)

  p <- plot(zeta)
  expect_s3_class(p, "ggplot")

  p2 <- plot(zeta, log_scale = TRUE)
  expect_s3_class(p2, "ggplot")
})


test_that("plot ratio works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 15, 1, 0.4), nrow = 20)

  zeta <- zetaDiversity(species, coords, orders = 1:5,
                         n_samples = 20, progress = FALSE)

  p <- plot(zeta, type = "ratio")
  expect_s3_class(p, "ggplot")
})


test_that("orders clipped to n_sites", {
  set.seed(42)
  coords <- data.frame(x = runif(5), y = runif(5))
  species <- matrix(rbinom(5 * 10, 1, 0.5), nrow = 5)

  zeta <- zetaDiversity(species, coords, orders = 1:10,
                         n_samples = 10, progress = FALSE)

  expect_equal(max(zeta$orders), 5)
})
