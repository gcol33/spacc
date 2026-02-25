test_that("spaccHillBeta returns correct structure", {
  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rpois(30 * 15, 2), nrow = 30)

  hb <- spaccHillBeta(species, coords, n_seeds = 5, progress = FALSE)

  expect_s3_class(hb, "spacc_hill_beta")
  expect_equal(hb$n_sites, 30)
  expect_equal(hb$n_seeds, 5)
  expect_equal(length(hb$gamma), 3)
  expect_equal(length(hb$alpha), 3)
  expect_equal(length(hb$beta), 3)
})


test_that("beta * alpha approx gamma identity", {
  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rpois(30 * 15, 2), nrow = 30)

  hb <- spaccHillBeta(species, coords, q = c(0, 1, 2),
                       n_seeds = 5, progress = FALSE)

  for (qi in 1:3) {
    gamma_mean <- colMeans(hb$gamma[[qi]])
    alpha_mean <- colMeans(hb$alpha[[qi]])
    beta_mean <- colMeans(hb$beta[[qi]])

    # beta * alpha should approximately equal gamma
    # Check at final site (most stable)
    n <- hb$n_sites
    reconstructed <- beta_mean[n] * alpha_mean[n]
    expect_equal(reconstructed, gamma_mean[n], tolerance = 0.15)
  }
})


test_that("gamma is non-decreasing for q=0", {
  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rpois(30 * 15, 2), nrow = 30)

  hb <- spaccHillBeta(species, coords, q = 0, n_seeds = 5, progress = FALSE)

  # Each seed's gamma curve should be non-decreasing (richness)
  for (s in 1:5) {
    diffs <- diff(hb$gamma[[1]][s, ])
    expect_true(all(diffs >= -1e-10))
  }
})


test_that("beta >= 1", {
  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rpois(30 * 15, 2), nrow = 30)

  hb <- spaccHillBeta(species, coords, q = c(0, 1, 2),
                       n_seeds = 5, progress = FALSE)

  for (qi in seq_along(hb$q)) {
    expect_true(all(hb$beta[[qi]] >= 1 - 1e-10))
  }
})


test_that("beta = 1 when all sites identical", {
  coords <- data.frame(x = runif(10), y = runif(10))
  # All sites have identical species composition
  species <- matrix(rep(c(5, 3, 2, 1, 0), each = 10), nrow = 10, ncol = 5)

  hb <- spaccHillBeta(species, coords, q = c(0, 1, 2),
                       n_seeds = 3, progress = FALSE)

  for (qi in seq_along(hb$q)) {
    # All beta values should be 1 (or very close)
    expect_true(all(abs(hb$beta[[qi]] - 1) < 0.1))
  }
})


test_that("print.spacc_hill_beta works", {
  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 2), nrow = 20)

  hb <- spaccHillBeta(species, coords, n_seeds = 5, progress = FALSE)
  expect_output(print(hb), "spacc Hill beta")
  expect_output(print(hb), "mean final beta")
})


test_that("summary.spacc_hill_beta returns data.frame", {
  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 2), nrow = 20)

  hb <- spaccHillBeta(species, coords, n_seeds = 5, progress = FALSE)
  summ <- summary(hb)

  expect_s3_class(summ, "data.frame")
  expect_true("component" %in% names(summ))
  expect_true(all(c("gamma", "alpha", "beta") %in% summ$component))
})


test_that("plot.spacc_hill_beta works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 2), nrow = 20)

  hb <- spaccHillBeta(species, coords, n_seeds = 5, progress = FALSE)

  p1 <- plot(hb, component = "beta")
  expect_s3_class(p1, "ggplot")

  p2 <- plot(hb, component = "all")
  expect_s3_class(p2, "ggplot")
})
