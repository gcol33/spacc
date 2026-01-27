test_that("time parameter produces valid spacc object", {
  set.seed(1)
  species <- matrix(rbinom(30 * 10, 1, 0.4), nrow = 30)
  coords <- data.frame(x = runif(30), y = runif(30))
  time <- rep(c(2020, 2021, 2022), each = 10)

  result <- spacc(species, coords, n_seeds = 5, method = "knn",
                  time = time, seed = 42, progress = FALSE)

  expect_s3_class(result, "spacc")
  expect_equal(nrow(result$curves), 5)
  expect_equal(ncol(result$curves), 30)
  expect_equal(result$time, time)
  expect_equal(result$w_space, 1)
  expect_equal(result$w_time, 1)
})


test_that("spatiotemporal differs from spatial-only", {
  set.seed(1)
  species <- matrix(rbinom(30 * 10, 1, 0.4), nrow = 30)
  coords <- data.frame(x = runif(30), y = runif(30))
  time <- rep(c(2020, 2021, 2022), each = 10)

  spatial <- spacc(species, coords, n_seeds = 5, method = "knn",
                   seed = 42, progress = FALSE)
  spatiotemporal <- spacc(species, coords, n_seeds = 5, method = "knn",
                          time = time, seed = 42, progress = FALSE)

  # Curves should differ (different distance matrices produce different orderings)
  expect_false(identical(spatial$curves, spatiotemporal$curves))
})


test_that("time parameter works with radius and gaussian methods", {
  set.seed(1)
  species <- matrix(rbinom(20 * 8, 1, 0.4), nrow = 20)
  coords <- data.frame(x = runif(20), y = runif(20))
  time <- rep(1:4, each = 5)

  result_radius <- spacc(species, coords, n_seeds = 3, method = "radius",
                         time = time, seed = 42, progress = FALSE)
  expect_s3_class(result_radius, "spacc")

  result_gauss <- spacc(species, coords, n_seeds = 3, method = "gaussian",
                        time = time, seed = 42, progress = FALSE)
  expect_s3_class(result_gauss, "spacc")
})


test_that("time parameter errors on unsupported methods", {
  species <- matrix(1, nrow = 10, ncol = 5)
  coords <- data.frame(x = 1:10, y = 1:10)
  time <- 1:10

  expect_error(
    spacc(species, coords, method = "kncn", time = time, progress = FALSE),
    "only supported for"
  )
  expect_error(
    spacc(species, coords, method = "random", time = time, progress = FALSE),
    "only supported for"
  )
  expect_error(
    spacc(species, coords, method = "cone", time = time, progress = FALSE),
    "only supported for"
  )
})


test_that("time parameter validates input", {
  species <- matrix(1, nrow = 10, ncol = 5)
  coords <- data.frame(x = 1:10, y = 1:10)

  expect_error(
    spacc(species, coords, time = 1:5, progress = FALSE),
    "length equal to nrow"
  )
  expect_error(
    spacc(species, coords, time = letters[1:10], progress = FALSE),
    "numeric"
  )
})


test_that("custom weights change results", {
  set.seed(1)
  species <- matrix(rbinom(30 * 10, 1, 0.4), nrow = 30)
  coords <- data.frame(x = runif(30), y = runif(30))
  time <- rep(c(2020, 2021, 2022), each = 10)

  w1 <- spacc(species, coords, n_seeds = 5, method = "knn",
              time = time, w_space = 1, w_time = 0.01, seed = 42, progress = FALSE)
  w2 <- spacc(species, coords, n_seeds = 5, method = "knn",
              time = time, w_space = 0.01, w_time = 1, seed = 42, progress = FALSE)

  # Different weights should produce different curves
  expect_false(identical(w1$curves, w2$curves))
})


test_that("time works with groups", {
  set.seed(1)
  species <- matrix(rbinom(30 * 10, 1, 0.4), nrow = 30)
  coords <- data.frame(x = runif(30), y = runif(30))
  time <- rep(c(2020, 2021, 2022), each = 10)
  groups <- rep(c("A", "B"), each = 5)

  result <- spacc(species, coords, n_seeds = 5, method = "knn",
                  time = time, groups = groups, seed = 42, progress = FALSE)

  expect_s3_class(result, "spacc")
  expect_true(spacc:::is_grouped(result))
  expect_equal(length(result$group_names), 2)
  expect_equal(result$time, time)
})
