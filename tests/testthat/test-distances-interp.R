test_that("distances.data.frame returns spacc_dist", {
  coords <- data.frame(x = runif(10), y = runif(10))

  result <- distances(coords)

  expect_s3_class(result, "spacc_dist")
  expect_equal(dim(result), c(10, 10))
  expect_equal(attr(result, "method"), "euclidean")
})


test_that("distances euclidean is symmetric and zero diagonal", {
  coords <- data.frame(x = runif(10), y = runif(10))

  d <- distances(coords, method = "euclidean")

  expect_equal(d, t(d))
  expect_equal(diag(d), rep(0, 10))
})


test_that("distances haversine works", {
  coords <- data.frame(x = runif(10, -5, 5), y = runif(10, 45, 50))

  d <- distances(coords, method = "haversine")

  expect_s3_class(d, "spacc_dist")
  expect_equal(attr(d, "method"), "haversine")
  expect_true(all(d >= 0))
})


test_that("distances euclidean matches base R dist", {
  set.seed(42)
  coords <- data.frame(x = runif(8), y = runif(8))

  d_spacc <- distances(coords, method = "euclidean")
  d_base <- as.matrix(dist(cbind(coords$x, coords$y)))

  expect_equal(as.numeric(d_spacc), as.numeric(d_base), tolerance = 1e-10)
})


test_that("distances with custom function", {
  coords <- data.frame(x = 1:5, y = rep(0, 5))

  # Manhattan distance
  manhattan <- function(x, y) {
    n <- length(x)
    m <- matrix(0, n, n)
    for (i in 1:n) for (j in 1:n) {
      m[i, j] <- abs(x[i] - x[j]) + abs(y[i] - y[j])
    }
    m
  }

  d <- distances(coords, fun = manhattan)

  expect_s3_class(d, "spacc_dist")
  expect_equal(d[1, 5], 4)
})


test_that("print.spacc_dist works", {
  coords <- data.frame(x = runif(10), y = runif(10))

  d <- distances(coords)

  expect_output(print(d), "spacc distance matrix")
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


test_that("spaccCoverage print works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 3), nrow = 15)

  cov <- spaccCoverage(species, coords, n_seeds = 3,
                        parallel = FALSE, progress = FALSE, seed = 1)

  expect_output(print(cov), "spacc coverage")
})


test_that("spaccCoverage summary returns data.frame", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 3), nrow = 15)

  cov <- spaccCoverage(species, coords, n_seeds = 3,
                        parallel = FALSE, progress = FALSE, seed = 1)

  summ <- summary(cov)
  expect_s3_class(summ, "data.frame")
  expect_true("mean_coverage" %in% names(summ))
  expect_equal(nrow(summ), 15)
})
