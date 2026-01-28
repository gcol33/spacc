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
