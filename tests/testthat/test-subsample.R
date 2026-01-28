test_that("subsample random returns correct number of sites", {
  coords <- data.frame(x = runif(50), y = runif(50))

  result <- subsample(coords, n = 20, method = "random", seed = 42)

  expect_length(result, 20)
  expect_true(all(result >= 1 & result <= 50))
  expect_equal(length(unique(result)), 20)
})


test_that("subsample grid returns indices", {
  coords <- data.frame(x = runif(50), y = runif(50))

  result <- subsample(coords, method = "grid", cell_size = 0.3, seed = 42)

  expect_true(length(result) > 0)
  expect_true(length(result) <= 50)
  expect_true(all(result >= 1 & result <= 50))
})


test_that("subsample grid with n target", {
  coords <- data.frame(x = runif(100), y = runif(100))

  result <- subsample(coords, n = 25, method = "grid", seed = 42)

  expect_true(length(result) > 0)
  expect_true(length(result) <= 100)
})


test_that("subsample thinning enforces minimum distance", {
  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))

  result <- subsample(coords, method = "thinning", min_dist = 0.2, seed = 42)

  # Check minimum distance among retained sites
  kept <- coords[result, ]
  if (length(result) > 1) {
    d <- as.matrix(dist(cbind(kept$x, kept$y)))
    diag(d) <- Inf
    expect_true(min(d) >= 0.2)
  }
})


test_that("subsample errors without required arguments", {
  coords <- data.frame(x = runif(20), y = runif(20))

  expect_error(subsample(coords, method = "random"), "n must be specified")
  expect_error(subsample(coords, method = "thinning"), "min_dist must be specified")
})


test_that("subsample errors on bad coords", {
  coords <- data.frame(a = 1:5, b = 1:5)

  expect_error(subsample(coords, n = 3, method = "random"), "x and y")
})


test_that("subsample random is reproducible with seed", {
  coords <- data.frame(x = runif(30), y = runif(30))

  r1 <- subsample(coords, n = 10, method = "random", seed = 42)
  r2 <- subsample(coords, n = 10, method = "random", seed = 42)

  expect_equal(r1, r2)
})
