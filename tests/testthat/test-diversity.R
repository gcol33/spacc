test_that("alphaDiversity returns correct structure", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- alphaDiversity(species, q = c(0, 1, 2))

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 20)
  expect_equal(ncol(result), 3)
})


test_that("alphaDiversity q=0 equals site richness", {
  set.seed(42)
  species <- matrix(rpois(15 * 8, 3), nrow = 15)

  result <- alphaDiversity(species, q = 0)

  # q=0 should equal number of species present at each site
  expected <- rowSums(species > 0)
  expect_equal(as.numeric(result[, 1]), expected)
})


test_that("alphaDiversity q=0 >= q=1 >= q=2", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- alphaDiversity(species, q = c(0, 1, 2))

  # Hill number inequality: q0 >= q1 >= q2
  for (i in 1:nrow(result)) {
    expect_true(result[i, 1] >= result[i, 2] - 1e-10,
                info = paste("Row", i, "q0 >= q1"))
    expect_true(result[i, 2] >= result[i, 3] - 1e-10,
                info = paste("Row", i, "q1 >= q2"))
  }
})


test_that("alphaDiversity with coords returns spacc_alpha", {
  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 3), nrow = 15)

  result <- alphaDiversity(species, q = c(0, 1, 2), coords = coords)

  expect_s3_class(result, "spacc_alpha")
})


test_that("gammaDiversity returns named vector", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- gammaDiversity(species, q = c(0, 1, 2))

  expect_true(is.numeric(result))
  expect_length(result, 3)
})


test_that("gammaDiversity q=0 equals total species count", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- gammaDiversity(species, q = 0)

  expected <- sum(colSums(species) > 0)
  expect_equal(as.numeric(result), expected)
})


test_that("gammaDiversity q=0 >= q=1 >= q=2", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- gammaDiversity(species, q = c(0, 1, 2))

  expect_true(result[1] >= result[2] - 1e-10)
  expect_true(result[2] >= result[3] - 1e-10)
})


test_that("diversityPartition returns correct structure", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- diversityPartition(species, q = c(0, 1, 2))

  expect_s3_class(result, "spacc_partition")
})


test_that("diversityPartition gamma = alpha * beta", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- diversityPartition(species, q = c(0, 1, 2))

  # Multiplicative partitioning: gamma = alpha * beta
  for (i in seq_along(result$q)) {
    computed_gamma <- result$alpha[i] * result$beta[i]
    expect_equal(computed_gamma, result$gamma[i], tolerance = 1e-6,
                 info = paste("q =", result$q[i]))
  }
})


test_that("diversityPartition beta >= 1", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- diversityPartition(species, q = c(0, 1, 2))

  expect_true(all(result$beta >= 1 - 1e-10))
})
