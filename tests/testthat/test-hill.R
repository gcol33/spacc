test_that("Hill number calculation is correct", {
  # Known values for simple community
  # 10 species with abundances 1:10
  abundances <- as.numeric(1:10)

  # q=0: Species richness = 10
  expect_equal(spacc:::calc_hill_number(abundances, 0), 10)

  # q=1: exp(Shannon)
  # Shannon H = -sum(p * log(p))
  N <- sum(abundances)
  p <- abundances / N
  H <- -sum(p * log(p))
  expect_equal(spacc:::calc_hill_number(abundances, 1), exp(H), tolerance = 1e-10)

  # q=2: Inverse Simpson
  # Simpson D = sum(p^2)
  D <- sum(p^2)
  expect_equal(spacc:::calc_hill_number(abundances, 2), 1/D, tolerance = 1e-10)
})


test_that("Hill number with equal abundances", {
  # All equal abundances: all Hill numbers should equal richness
  abundances <- rep(10, 5)

  expect_equal(spacc:::calc_hill_number(abundances, 0), 5)
  expect_equal(spacc:::calc_hill_number(abundances, 1), 5, tolerance = 1e-10)
  expect_equal(spacc:::calc_hill_number(abundances, 2), 5, tolerance = 1e-10)
})


test_that("Hill number handles edge cases", {
  # Empty community
  expect_equal(spacc:::calc_hill_number(numeric(0), 0), 0)
  expect_equal(spacc:::calc_hill_number(c(0, 0, 0), 0), 0)

  # Single species
  expect_equal(spacc:::calc_hill_number(c(10), 0), 1)
  expect_equal(spacc:::calc_hill_number(c(10), 1), 1)
  expect_equal(spacc:::calc_hill_number(c(10), 2), 1)
})


test_that("spaccHill returns correct structure", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 2), nrow = 20)

  result <- spaccHill(species, coords, q = c(0, 1, 2), n_seeds = 5,
                      parallel = FALSE, progress = FALSE)

  expect_s3_class(result, "spacc_hill")
  expect_equal(result$n_sites, 20)
  expect_equal(result$n_species, 10)
  expect_equal(result$n_seeds, 5)
  expect_equal(result$q, c(0, 1, 2))

  # Should have 3 curve matrices
  expect_equal(length(result$curves), 3)

  # Each matrix should be n_seeds x n_sites
  expect_equal(dim(result$curves[[1]]), c(5, 20))
})


test_that("Species richness (q=0) is monotonically non-decreasing", {
  skip_on_cran()

  set.seed(123)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rpois(30 * 15, 3), nrow = 30)

  result <- spaccHill(species, coords, q = c(0, 1, 2), n_seeds = 3,
                      parallel = FALSE, progress = FALSE)

  # Only q=0 (richness) is guaranteed to be monotonically non-decreasing.
  # For q>0, Hill numbers depend on evenness and can decrease when adding
  # sites that shift dominance.
  q0_idx <- which(result$q == 0)
  for (seed in 1:3) {
    curve <- result$curves[[q0_idx]][seed, ]
    diffs <- diff(curve)
    expect_true(all(diffs >= -1e-10),
                info = paste("q = 0, seed =", seed))
  }
})
