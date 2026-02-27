test_that("Coverage calculation is correct", {
  # Simple case: 10 species, each with abundance 10
  # No singletons, no doubletons -> coverage should be 1
  abundances <- rep(10, 10)
  expect_equal(spacc:::calc_coverage(abundances), 1)

  # All singletons: coverage should be low
  abundances_single <- rep(1, 10)
  cov <- spacc:::calc_coverage(abundances_single)
  expect_true(cov < 1)
  expect_true(cov > 0)
})


test_that("Coverage increases with accumulation", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rpois(30 * 15, 2), nrow = 30)

  result <- spaccCoverage(species, coords, n_seeds = 3,
                          parallel = FALSE, progress = FALSE)

  # Coverage should generally increase (allow small decreases due to singletons)
  for (seed in 1:3) {
    curve <- result$coverage[seed, ]
    final_cov <- curve[length(curve)]
    initial_cov <- curve[1]
    expect_true(final_cov >= initial_cov)
  }
})


test_that("spaccCoverage returns correct structure", {
  skip_on_cran()

  set.seed(123)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- spaccCoverage(species, coords, n_seeds = 5,
                          parallel = FALSE, progress = FALSE)

  expect_s3_class(result, "spacc_coverage")
  expect_equal(result$n_sites, 20)
  expect_equal(result$n_seeds, 5)

  expect_equal(dim(result$richness), c(5, 20))
  expect_equal(dim(result$individuals), c(5, 20))
  expect_equal(dim(result$coverage), c(5, 20))

  # Coverage should be between 0 and 1
  expect_true(all(result$coverage >= 0 & result$coverage <= 1))
})


test_that("interpolateCoverage works correctly", {
  skip_on_cran()

  set.seed(456)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rpois(30 * 15, 3), nrow = 30)

  cov_result <- spaccCoverage(species, coords, n_seeds = 5,
                               parallel = FALSE, progress = FALSE)

  interp <- interpolateCoverage(cov_result, target = c(0.8, 0.9, 0.95))

  expect_equal(ncol(interp), 3)
  expect_equal(nrow(interp), 5)

  # Interpolated values should increase with target coverage
  for (seed in 1:5) {
    if (!any(is.na(interp[seed, ]))) {
      expect_true(interp[seed, 1] <= interp[seed, 2])
      expect_true(interp[seed, 2] <= interp[seed, 3])
    }
  }
})


test_that("spaccCoverage with spacc_dist coords", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 3), nrow = 15)

  d <- distances(coords)
  result <- spaccCoverage(species, d, n_seeds = 3,
                          parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc_coverage")
  expect_equal(result$n_seeds, 3)
})


test_that("spaccCoverage print works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 3), nrow = 15)

  result <- spaccCoverage(species, coords, n_seeds = 3,
                          parallel = FALSE, progress = FALSE, seed = 1)

  expect_output(print(result), "spacc coverage")
})


test_that("spaccCoverage summary returns data.frame", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 3), nrow = 15)

  result <- spaccCoverage(species, coords, n_seeds = 3,
                          parallel = FALSE, progress = FALSE, seed = 1)

  summ <- summary(result)
  expect_s3_class(summ, "data.frame")
  expect_true("mean_coverage" %in% names(summ))
  expect_equal(nrow(summ), 15)
})


test_that("extrapolateCoverage returns correct structure", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  cov <- spaccCoverage(species, coords, n_seeds = 3,
                        parallel = FALSE, progress = FALSE, seed = 1)

  result <- extrapolateCoverage(cov, target_coverage = c(0.95, 0.99))

  expect_s3_class(result, "spacc_coverage_ext")
  expect_equal(result$q, 0)
  expect_equal(result$target_coverage, c(0.95, 0.99))
  expect_equal(dim(result$richness), c(3, 2))
})


test_that("extrapolateCoverage print works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  cov <- spaccCoverage(species, coords, n_seeds = 3,
                        parallel = FALSE, progress = FALSE, seed = 1)

  result <- extrapolateCoverage(cov, target_coverage = c(0.95, 0.99))

  expect_output(print(result), "Coverage-based extrapolation")
})


test_that("interpolateCoverage returns correct structure", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  cov <- spaccCoverage(species, coords, n_seeds = 3,
                        parallel = FALSE, progress = FALSE, seed = 1)

  result <- interpolateCoverage(cov, target = c(0.90, 0.95))

  expect_s3_class(result, "data.frame")
  expect_equal(ncol(result), 2)
  expect_equal(nrow(result), 3)
  expect_equal(names(result), c("C90", "C95"))
})


test_that("extrapolateCoverage returns correct structure", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  cov <- spaccCoverage(species, coords, n_seeds = 3,
                        parallel = FALSE, progress = FALSE, seed = 1)

  ext <- extrapolateCoverage(cov, target_coverage = c(0.95, 0.99))

  expect_s3_class(ext, "spacc_coverage_ext")
  expect_equal(ext$target_coverage, c(0.95, 0.99))
  expect_equal(nrow(ext$richness), 3)
  expect_equal(ncol(ext$richness), 2)
})


test_that("extrapolateCoverage with q=1 works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  cov <- spaccCoverage(species, coords, n_seeds = 3,
                        parallel = FALSE, progress = FALSE, seed = 1)

  ext <- extrapolateCoverage(cov, target_coverage = c(0.99), q = 1)

  expect_s3_class(ext, "spacc_coverage_ext")
  expect_equal(ext$q, 1)
})


test_that("extrapolateCoverage with q=2 works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  cov <- spaccCoverage(species, coords, n_seeds = 3,
                        parallel = FALSE, progress = FALSE, seed = 1)

  ext <- extrapolateCoverage(cov, target_coverage = c(0.99), q = 2)

  expect_s3_class(ext, "spacc_coverage_ext")
  expect_equal(ext$q, 2)
})


test_that("extrapolateCoverage print works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  cov <- spaccCoverage(species, coords, n_seeds = 3,
                        parallel = FALSE, progress = FALSE, seed = 1)

  ext <- extrapolateCoverage(cov, target_coverage = c(0.95, 0.99))

  expect_output(print(ext), "Coverage-based extrapolation")
})


test_that("extrapolateCoverage summary returns data.frame", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  cov <- spaccCoverage(species, coords, n_seeds = 3,
                        parallel = FALSE, progress = FALSE, seed = 1)

  ext <- extrapolateCoverage(cov, target_coverage = c(0.95, 0.99))

  summ <- summary(ext)
  expect_s3_class(summ, "data.frame")
  expect_true("mean_richness" %in% names(summ))
})


test_that("extrapolateCoverage plot returns ggplot", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  cov <- spaccCoverage(species, coords, n_seeds = 3,
                        parallel = FALSE, progress = FALSE, seed = 1)

  ext <- extrapolateCoverage(cov, target_coverage = c(0.95, 0.99))

  p <- plot(ext)
  expect_s3_class(p, "ggplot")
})


test_that("plot.spacc_coverage with different xaxis works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  cov <- spaccCoverage(species, coords, n_seeds = 3,
                        parallel = FALSE, progress = FALSE, seed = 1)

  p1 <- plot(cov, xaxis = "coverage")
  expect_s3_class(p1, "ggplot")

  p2 <- plot(cov, xaxis = "individuals")
  expect_s3_class(p2, "ggplot")
})


test_that("plot.spacc_coverage map errors without map data", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  cov <- spaccCoverage(species, coords, n_seeds = 3,
                        parallel = FALSE, progress = FALSE, seed = 1)

  expect_error(plot(cov, type = "map"), "map")
})


test_that("as_sf.spacc_coverage errors without map data", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  cov <- spaccCoverage(species, coords, n_seeds = 3,
                        parallel = FALSE, progress = FALSE, seed = 1)

  expect_error(as_sf(cov), "map")
})


# ============================================================================
# Chiu (2023) sample-based coverage estimator
# ============================================================================

test_that("calc_chiu_coverage returns correct values", {
  # Q1=2, Q2=1, G1=8, n+=100, T=10
  # a0_hat = (9/10) * 8 * 2 / (2*1) = 7.2
  # C = 100 / 107.2 ~ 0.9328
  abundances <- c(5, 3, 20, 30, 42)
  incidences <- as.integer(c(1, 1, 2, 5, 8))
  C <- spacc:::calc_chiu_coverage(abundances, incidences, 10L)
  expect_equal(C, 100 / 107.2, tolerance = 0.001)
})


test_that("calc_chiu_coverage handles edge cases", {
  # No uniques -> coverage = 1
  C <- spacc:::calc_chiu_coverage(c(10, 20), as.integer(c(3, 5)), 5L)
  expect_equal(C, 1.0)

  # Single site -> coverage = 0
  C <- spacc:::calc_chiu_coverage(c(5, 3), as.integer(c(1, 1)), 1L)
  expect_equal(C, 0.0)
})


test_that("spaccCoverage with coverage='chiu' returns correct structure", {
  skip_on_cran()

  set.seed(123)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- spaccCoverage(species, coords, n_seeds = 5,
                          coverage = "chiu",
                          parallel = FALSE, progress = FALSE)

  expect_s3_class(result, "spacc_coverage")
  expect_equal(result$coverage_type, "chiu")
  expect_equal(result$n_sites, 20)
  expect_equal(result$n_seeds, 5)
  expect_equal(dim(result$richness), c(5, 20))
  expect_equal(dim(result$coverage), c(5, 20))

  # Coverage should be between 0 and 1

  expect_true(all(result$coverage >= 0 & result$coverage <= 1))
})


test_that("spaccCoverage chiu coverage increases with accumulation", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rpois(30 * 15, 2), nrow = 30)

  result <- spaccCoverage(species, coords, n_seeds = 3,
                          coverage = "chiu",
                          parallel = FALSE, progress = FALSE)

  # Final coverage should be >= initial coverage
  for (s in 1:3) {
    curve <- result$coverage[s, ]
    expect_true(curve[length(curve)] >= curve[1])
  }
})


test_that("spaccCoverage default coverage is 'chao'", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 3), nrow = 15)

  result <- spaccCoverage(species, coords, n_seeds = 3,
                          parallel = FALSE, progress = FALSE)

  expect_equal(result$coverage_type, "chao")
})


test_that("spaccCoverage print shows coverage type", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 3), nrow = 15)

  result_chiu <- spaccCoverage(species, coords, n_seeds = 3,
                                coverage = "chiu",
                                parallel = FALSE, progress = FALSE)

  expect_output(print(result_chiu), "Chiu 2023")

  result_chao <- spaccCoverage(species, coords, n_seeds = 3,
                                coverage = "chao",
                                parallel = FALSE, progress = FALSE)

  expect_output(print(result_chao), "Chao-Jost 2012")
})


test_that("interpolateCoverage works with chiu coverage", {
  skip_on_cran()

  set.seed(456)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rpois(30 * 15, 3), nrow = 30)

  cov_result <- spaccCoverage(species, coords, n_seeds = 5,
                               coverage = "chiu",
                               parallel = FALSE, progress = FALSE)

  interp <- interpolateCoverage(cov_result, target = c(0.8, 0.9, 0.95))

  expect_equal(ncol(interp), 3)
  expect_equal(nrow(interp), 5)
})


test_that("spaccCoverage chiu works with parallel", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- spaccCoverage(species, coords, n_seeds = 5,
                          coverage = "chiu",
                          parallel = TRUE, n_cores = 2,
                          progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc_coverage")
  expect_true(all(result$coverage >= 0 & result$coverage <= 1))
})
