test_that("rarefy returns correct structure", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- rarefy(species, n_boot = 10)

  expect_s3_class(result, "spacc_rare")
  expect_true(length(result$n) > 0)
  expect_equal(length(result$expected), length(result$n))
  expect_equal(length(result$sd), length(result$n))
  expect_equal(length(result$lower), length(result$n))
  expect_equal(length(result$upper), length(result$n))
  expect_equal(result$n_boot, 10)
})


test_that("rarefy expected richness is monotonically non-decreasing", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- rarefy(species, n_boot = 10)

  diffs <- diff(result$expected)
  expect_true(all(diffs >= -1e-10))
})


test_that("rarefy at total individuals returns observed richness", {
  set.seed(42)
  species <- matrix(rpois(15 * 8, 3), nrow = 15)

  result <- rarefy(species, n_individuals = sum(species), n_boot = 10)

  observed <- sum(colSums(species) > 0)
  expect_equal(result$expected[length(result$expected)], observed)
})


test_that("rarefy with custom n_individuals", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- rarefy(species, n_individuals = c(10, 50, 100), n_boot = 10)

  expect_equal(result$n, c(10, 50, 100))
  expect_length(result$expected, 3)
})


test_that("rarefy q=1 Shannon diversity", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- rarefy(species, q = 1, n_boot = 10)

  expect_s3_class(result, "spacc_rare")
  expect_equal(result$q, 1)
  # Shannon diversity should be >= 1 for non-trivial data
  expect_true(all(result$expected >= 1))
  # Should be monotonically non-decreasing
  expect_true(all(diff(result$expected) >= -0.01))
})


test_that("rarefy q=2 Simpson diversity", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- rarefy(species, q = 2, n_boot = 10)

  expect_s3_class(result, "spacc_rare")
  expect_equal(result$q, 2)
  # Simpson diversity should be >= 1
  expect_true(all(result$expected >= 1))
})


test_that("rarefy q=0 backward compatible", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- rarefy(species, n_boot = 10)
  result_q0 <- rarefy(species, q = 0, n_boot = 10)

  # Default should behave same as q=0
  # (bootstrap introduces randomness, but expected should match)
  expect_equal(result$expected, result_q0$expected)
})


test_that("Hill ordering: D0 >= D1 >= D2 at full sample", {
  set.seed(42)
  species <- matrix(rpois(30 * 15, 3), nrow = 30)
  n_total <- sum(colSums(species)[colSums(species) > 0])

  r0 <- rarefy(species, n_individuals = n_total, q = 0, n_boot = 5)
  r1 <- rarefy(species, n_individuals = n_total, q = 1, n_boot = 5)
  r2 <- rarefy(species, n_individuals = n_total, q = 2, n_boot = 5)

  expect_true(r0$expected >= r1$expected - 0.1)
  expect_true(r1$expected >= r2$expected - 0.1)
})


test_that("print.spacc_rare works", {
  set.seed(42)
  species <- matrix(rpois(15 * 8, 3), nrow = 15)

  result <- rarefy(species, n_boot = 10)

  expect_output(print(result), "Individual-based rarefaction")
})


test_that("as.data.frame.spacc_rare works", {
  set.seed(42)
  species <- matrix(rpois(15 * 8, 3), nrow = 15)

  result <- rarefy(species, n_boot = 10)
  df <- as.data.frame(result)

  expect_s3_class(df, "data.frame")
  expect_equal(names(df), c("n", "expected", "sd", "lower", "upper"))
  expect_equal(nrow(df), length(result$n))
})
