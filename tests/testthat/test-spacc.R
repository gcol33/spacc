test_that("spacc returns correct structure with knn", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.4), nrow = 20)

  result <- spacc(species, coords, n_seeds = 3, method = "knn",
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
  expect_equal(result$n_sites, 20)
  expect_equal(result$n_species, ncol(species))
  expect_equal(result$n_seeds, 3)
  expect_equal(result$method, "knn")
  expect_equal(dim(result$curves), c(3, 20))
})


test_that("spacc works with kncn method", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  result <- spacc(species, coords, n_seeds = 3, method = "kncn",
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
  expect_equal(result$method, "kncn")
  expect_equal(dim(result$curves), c(3, 15))
})


test_that("spacc works with random method", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  result <- spacc(species, coords, n_seeds = 3, method = "random",
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
  expect_equal(result$method, "random")
})


test_that("spacc works with gaussian method", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  result <- spacc(species, coords, n_seeds = 3, method = "gaussian",
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
  expect_equal(result$method, "gaussian")
})


test_that("spacc works with cone method", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  result <- spacc(species, coords, n_seeds = 3, method = "cone",
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
  expect_equal(result$method, "cone")
})


test_that("spacc works with collector method", {
  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  result <- spacc(species, coords, method = "collector",
                  parallel = FALSE, progress = FALSE)

  expect_s3_class(result, "spacc")
  expect_equal(result$method, "collector")
  # Collector produces a single curve
  expect_equal(result$n_seeds, 1)
})


test_that("spacc curves are monotonically non-decreasing", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.4), nrow = 20)

  result <- spacc(species, coords, n_seeds = 5, method = "knn",
                  parallel = FALSE, progress = FALSE, seed = 1)

  for (s in 1:5) {
    diffs <- diff(result$curves[s, ])
    expect_true(all(diffs >= 0), info = paste("Seed", s))
  }
})


test_that("spacc works with haversine distance", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15, -5, 5), y = runif(15, 45, 50))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  result <- spacc(species, coords, n_seeds = 3, method = "knn",
                  distance = "haversine",
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
  expect_equal(result$distance, "haversine")
})


test_that("spacc accepts spacc_dist object", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  d <- distances(coords)
  result <- spacc(species, d, n_seeds = 3, method = "knn",
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
})


test_that("spacc seed argument produces reproducible results", {
  skip_on_cran()

  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  r1 <- spacc(species, coords, n_seeds = 3, method = "knn",
              parallel = FALSE, progress = FALSE, seed = 42)
  r2 <- spacc(species, coords, n_seeds = 3, method = "knn",
              parallel = FALSE, progress = FALSE, seed = 42)

  expect_equal(r1$curves, r2$curves)
})


test_that("spacc works with abundance data", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 3), nrow = 15)

  result <- spacc(species, coords, n_seeds = 3, method = "knn",
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
})


test_that("spacc works with radius method", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  result <- spacc(species, coords, n_seeds = 3, method = "radius",
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
  expect_equal(result$method, "radius")
})


test_that("spacc works with kdtree backend", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.4), nrow = 20)

  result <- spacc(species, coords, n_seeds = 3, method = "knn",
                  backend = "kdtree",
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
  expect_equal(result$backend, "kdtree")
})


test_that("spacc works with exact backend", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  result <- spacc(species, coords, n_seeds = 3, method = "knn",
                  backend = "exact",
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
  expect_equal(result$backend, "exact")
})


test_that("spacc kncn with kdtree backend", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  result <- spacc(species, coords, n_seeds = 3, method = "kncn",
                  backend = "kdtree",
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
})


test_that("spacc with custom sigma for gaussian", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  result <- spacc(species, coords, n_seeds = 3, method = "gaussian",
                  sigma = 0.5,
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
  expect_equal(result$sigma, 0.5)
})


test_that("spacc with custom cone_width", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  result <- spacc(species, coords, n_seeds = 3, method = "cone",
                  cone_width = pi/6,
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
  expect_equal(result$cone_width, pi/6)
})


test_that("spacc with data.frame input converts to matrix", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- data.frame(matrix(rbinom(15 * 8, 1, 0.4), nrow = 15))

  result <- spacc(species, coords, n_seeds = 3, method = "knn",
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
})
