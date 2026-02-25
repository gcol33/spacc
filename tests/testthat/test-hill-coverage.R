test_that("spaccHillCoverage returns correct structure", {
  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rpois(30 * 15, 2), nrow = 30)

  hc <- spaccHillCoverage(species, coords, n_seeds = 5, progress = FALSE)

  expect_s3_class(hc, "spacc_hill_coverage")
  expect_equal(hc$n_sites, 30)
  expect_equal(hc$n_seeds, 5)
  expect_equal(length(hc$hills), 3)  # q = 0, 1, 2
  expect_equal(nrow(hc$coverage), 5)
  expect_equal(ncol(hc$coverage), 30)
})


test_that("coverage values in [0, 1]", {
  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rpois(30 * 15, 2), nrow = 30)

  hc <- spaccHillCoverage(species, coords, n_seeds = 5, progress = FALSE)

  expect_true(all(hc$coverage >= 0))
  expect_true(all(hc$coverage <= 1))
})


test_that("q=0 curve is monotone non-decreasing", {
  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rpois(30 * 15, 2), nrow = 30)

  hc <- spaccHillCoverage(species, coords, q = 0, n_seeds = 5, progress = FALSE)

  # Check each seed's curve is non-decreasing (richness)
  for (s in 1:5) {
    diffs <- diff(hc$hills[[1]][s, ])
    expect_true(all(diffs >= -1e-10))
  }
})


test_that("q=2 <= q=1 <= q=0 at final step", {
  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rpois(30 * 15, 2), nrow = 30)

  hc <- spaccHillCoverage(species, coords, q = c(0, 1, 2),
                           n_seeds = 5, progress = FALSE)

  # Mean final values
  q0_final <- mean(hc$hills[[1]][, 30])
  q1_final <- mean(hc$hills[[2]][, 30])
  q2_final <- mean(hc$hills[[3]][, 30])

  expect_true(q0_final >= q1_final - 1e-10)
  expect_true(q1_final >= q2_final - 1e-10)
})


test_that("target_coverage returns standardized slot", {
  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rpois(30 * 15, 2), nrow = 30)

  hc <- spaccHillCoverage(species, coords, q = c(0, 1),
                           target_coverage = 0.9,
                           n_seeds = 5, progress = FALSE)

  expect_false(is.null(hc$standardized))
  expect_equal(length(hc$standardized), 2)  # one per q
})


test_that("print.spacc_hill_coverage works", {
  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 2), nrow = 20)

  hc <- spaccHillCoverage(species, coords, n_seeds = 5, progress = FALSE)
  expect_output(print(hc), "spacc Hill")
  expect_output(print(hc), "coverage")
})


test_that("summary.spacc_hill_coverage returns data.frame", {
  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 2), nrow = 20)

  hc <- spaccHillCoverage(species, coords, n_seeds = 5, progress = FALSE)
  summ <- summary(hc)

  expect_s3_class(summ, "data.frame")
  expect_true("mean_hill" %in% names(summ))
  expect_true("mean_coverage" %in% names(summ))
})


test_that("plot.spacc_hill_coverage works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 2), nrow = 20)

  hc <- spaccHillCoverage(species, coords, n_seeds = 5, progress = FALSE)

  p1 <- plot(hc, xaxis = "coverage")
  expect_s3_class(p1, "ggplot")

  p2 <- plot(hc, xaxis = "sites")
  expect_s3_class(p2, "ggplot")
})
