test_that("coleman returns correct structure", {
  set.seed(42)
  species <- matrix(rbinom(20 * 10, 1, 0.4), nrow = 20)

  result <- coleman(species)

  expect_s3_class(result, "data.frame")
  expect_equal(names(result), c("sites", "expected", "sd"))
  expect_equal(nrow(result), 20)
  expect_equal(result$sites, 1:20)
})


test_that("coleman expected richness is monotonically non-decreasing", {
  set.seed(42)
  species <- matrix(rbinom(30 * 15, 1, 0.4), nrow = 30)

  result <- coleman(species)

  diffs <- diff(result$expected)
  expect_true(all(diffs >= -1e-10))
})


test_that("coleman at N sites equals observed richness", {
  set.seed(42)
  species <- matrix(rbinom(20 * 10, 1, 0.4), nrow = 20)

  result <- coleman(species)
  observed <- sum(colSums(species) > 0)

  expect_equal(result$expected[20], observed, tolerance = 1e-10)
})


test_that("mao_tau returns correct structure", {
  set.seed(42)
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  result <- mao_tau(species)

  expect_s3_class(result, "data.frame")
  expect_equal(names(result), c("sites", "expected", "sd", "lower", "upper"))
  expect_equal(nrow(result), 15)
})


test_that("mao_tau is monotonically non-decreasing", {
  set.seed(42)
  species <- matrix(rbinom(20 * 10, 1, 0.4), nrow = 20)

  result <- mao_tau(species)

  diffs <- diff(result$expected)
  expect_true(all(diffs >= -1e-10))
})


test_that("mao_tau at N sites equals observed richness", {
  set.seed(42)
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  result <- mao_tau(species)
  observed <- sum(colSums(species) > 0)

  expect_equal(result$expected[15], observed, tolerance = 1e-10)
})


test_that("collector returns correct structure", {
  set.seed(42)
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  result <- collector(species)

  expect_s3_class(result, "data.frame")
  expect_equal(names(result), c("sites", "species"))
  expect_equal(nrow(result), 15)
})


test_that("collector is monotonically non-decreasing", {
  set.seed(42)
  species <- matrix(rbinom(20 * 10, 1, 0.4), nrow = 20)

  result <- collector(species)

  diffs <- diff(result$species)
  expect_true(all(diffs >= 0))
})


test_that("collector final richness equals observed", {
  set.seed(42)
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  result <- collector(species)
  observed <- sum(colSums(species) > 0)

  expect_equal(result$species[15], observed)
})


test_that("spatialRarefaction returns correct structure", {
  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  result <- spatialRarefaction(species, coords, n_perm = 10)

  expect_s3_class(result, "data.frame")
  expect_equal(names(result), c("sites", "mean", "sd", "lower", "upper"))
  expect_equal(nrow(result), 15)
})


test_that("spatialRarefaction mean is non-decreasing", {
  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.4), nrow = 20)

  result <- spatialRarefaction(species, coords, n_perm = 20)

  # Mean curve should be approximately non-decreasing (stochastic, allow small dips)
  diffs <- diff(result$mean)
  expect_true(sum(diffs < -0.5) == 0)
})
