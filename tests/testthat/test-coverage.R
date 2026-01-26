test_that("Coverage calculation is correct", {
  # Simple case: 10 species, each with abundance 10
  # No singletons, no doubletons -> coverage should be 1
  abundances <- rep(10, 10)
  expect_equal(calc_coverage(abundances), 1)

  # All singletons: coverage should be low
  abundances_single <- rep(1, 10)
  cov <- calc_coverage(abundances_single)
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
