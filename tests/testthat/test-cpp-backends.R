# Tests to improve C++ code coverage

# ============================================================================
# Test kNCN backend (src/kncn.cpp)
# ============================================================================

test_that("spacc with kncn method works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  result <- spacc(species, coords, n_seeds = 5, method = "kncn",
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
  expect_equal(result$method, "kncn")
})


test_that("spacc with kncn method with multiple seeds", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  result <- spacc(species, coords, n_seeds = 5, method = "kncn",
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
  expect_equal(result$n_seeds, 5)
})


test_that("spacc with kncn method gives different results than knn", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(25), y = runif(25))
  species <- matrix(rbinom(25 * 12, 1, 0.3), nrow = 25)

  result_knn <- spacc(species, coords, n_seeds = 3, method = "knn",
                      parallel = FALSE, progress = FALSE, seed = 1)
  result_kncn <- spacc(species, coords, n_seeds = 3, method = "kncn",
                       parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result_knn, "spacc")
  expect_s3_class(result_kncn, "spacc")
  # kncn and knn should produce different curves
  expect_false(identical(result_knn$curves, result_kncn$curves))
})


# ============================================================================
# Test random backend (src/random.cpp)
# ============================================================================

test_that("spacc with random method works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  result <- spacc(species, coords, n_seeds = 5, method = "random",
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
  expect_equal(result$method, "random")
})


test_that("spacc with random method produces different curves than knn", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  result_knn <- spacc(species, coords, n_seeds = 5, method = "knn",
                      parallel = FALSE, progress = FALSE, seed = 1)
  result_rand <- spacc(species, coords, n_seeds = 5, method = "random",
                       parallel = FALSE, progress = FALSE, seed = 1)

  # Random should be different from knn (not spatially structured)
  expect_false(identical(result_knn$curves, result_rand$curves))
})


# ============================================================================
# Test spatial tree backends - kdtree and exact
# ============================================================================

test_that("spacc with kdtree backend works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  result <- spacc(species, coords, n_seeds = 3, method = "knn",
                  backend = "kdtree",
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
  expect_equal(result$backend, "kdtree")
})


test_that("spacc with exact backend works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  result <- spacc(species, coords, n_seeds = 3, method = "knn",
                  backend = "exact",
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
  expect_equal(result$backend, "exact")
})


test_that("spacc with auto backend works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  result <- spacc(species, coords, n_seeds = 3, method = "knn",
                  backend = "auto",
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
  expect_true(result$backend %in% c("exact", "kdtree"))
})


test_that("exact and kdtree backends give same results", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(25), y = runif(25))
  species <- matrix(rbinom(25 * 12, 1, 0.3), nrow = 25)

  result_exact <- spacc(species, coords, n_seeds = 3, method = "knn",
                        backend = "exact",
                        parallel = FALSE, progress = FALSE, seed = 1)
  result_tree <- spacc(species, coords, n_seeds = 3, method = "knn",
                       backend = "kdtree",
                       parallel = FALSE, progress = FALSE, seed = 1)

  # Results should be essentially identical
  expect_equal(result_exact$curves, result_tree$curves, tolerance = 1e-10)
})


# ============================================================================
# Test beta C++ (src/beta.cpp) - with jaccard index
# ============================================================================

test_that("spaccBeta with jaccard index works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  result <- spaccBeta(species, coords, n_seeds = 3, index = "jaccard",
                      parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc_beta")
  expect_equal(result$index, "jaccard")
})


# ============================================================================
# Test coverage C++ (src/coverage.cpp) - edge cases
# ============================================================================

test_that("spaccCoverage with sparse data works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  # Very sparse abundance data
  species <- matrix(rpois(20 * 10, 0.5), nrow = 20)

  result <- spaccCoverage(species, coords, n_seeds = 3,
                          parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc_coverage")
  expect_true(all(result$coverage >= 0 & result$coverage <= 1))
})


test_that("spaccCoverage with high abundance data works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  # High abundance data
  species <- matrix(rpois(20 * 10, 50), nrow = 20)

  result <- spaccCoverage(species, coords, n_seeds = 3,
                          parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc_coverage")
  # High abundance should give high coverage
  expect_true(mean(result$coverage[, 20]) > 0.9)
})


# ============================================================================
# Test knn with haversine distance (src/knn.cpp)
# ============================================================================

test_that("spacc with haversine distance works", {
  skip_on_cran()

  set.seed(42)
  # Lat/lon coordinates
  coords <- data.frame(x = runif(20, -5, 5), y = runif(20, 45, 50))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  result <- spacc(species, coords, n_seeds = 3, distance = "haversine",
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
  expect_equal(result$distance, "haversine")
})


# ============================================================================
# Test parallel execution (multiple seeds)
# ============================================================================

test_that("spacc with many seeds for parallel coverage", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(25), y = runif(25))
  species <- matrix(rbinom(25 * 12, 1, 0.3), nrow = 25)

  result <- spacc(species, coords, n_seeds = 10,
                  parallel = TRUE, n_cores = 2,
                  progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
  expect_equal(result$n_seeds, 10)
  expect_equal(nrow(result$curves), 10)
})


test_that("spaccHill with many seeds works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 2), nrow = 20)

  result <- spaccHill(species, coords, q = c(0, 1), n_seeds = 10,
                      parallel = TRUE, n_cores = 2,
                      progress = FALSE)

  expect_s3_class(result, "spacc_hill")
  expect_equal(result$n_seeds, 10)
})


test_that("spaccBeta with many seeds works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  result <- spaccBeta(species, coords, n_seeds = 10,
                      parallel = TRUE, n_cores = 2,
                      progress = FALSE)

  expect_s3_class(result, "spacc_beta")
  expect_equal(result$n_seeds, 10)
})
