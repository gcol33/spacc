test_that("chao1 returns correct structure", {
  set.seed(42)
  species <- matrix(rpois(50 * 30, 2), nrow = 50)

  result <- chao1(species)

  expect_s3_class(result, "spacc_estimate")
  expect_equal(result$estimator, "chao1")
  expect_true(is.numeric(result$estimate))
  expect_true(is.numeric(result$se))
  expect_true(is.numeric(result$lower))
  expect_true(is.numeric(result$upper))
  expect_true(is.numeric(result$S_obs))
})


test_that("chao1 estimate >= S_obs", {
  set.seed(42)
  species <- matrix(rpois(50 * 30, 2), nrow = 50)

  result <- chao1(species)
  expect_true(result$estimate >= result$S_obs)
})


test_that("chao1 handles f2 = 0 edge case", {
  # Create data where all species have abundance >= 2
  species <- matrix(0, nrow = 10, ncol = 5)
  species[1:3, 1] <- 1  # f1 = 1 (singleton)
  species[1:5, 2] <- 3
  species[1:8, 3] <- 4
  species[1:10, 4] <- 5
  species[1:10, 5] <- 6
  # Pool: sp1 = 3, sp2 = 15, sp3 = 32, sp4 = 50, sp5 = 60 -> f2 = 0

  result <- chao1(species)
  expect_s3_class(result, "spacc_estimate")
  expect_true(result$estimate >= result$S_obs)
})


test_that("chao2 returns correct structure", {
  set.seed(42)
  species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)

  result <- chao2(species)

  expect_s3_class(result, "spacc_estimate")
  expect_equal(result$estimator, "chao2")
  expect_true(result$estimate >= result$S_obs)
})


test_that("chao2 uses incidence data", {
  set.seed(42)
  species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)

  result <- chao2(species)
  expect_true("Q1" %in% names(result$details))
  expect_true("Q2" %in% names(result$details))
})


test_that("ace returns correct structure", {
  set.seed(42)
  species <- matrix(rpois(50 * 30, 2), nrow = 50)

  result <- ace(species)

  expect_s3_class(result, "spacc_estimate")
  expect_equal(result$estimator, "ace")
  expect_true(result$estimate >= result$S_obs)
})


test_that("ace handles no rare species", {
  # All species have high abundance
  species <- matrix(rpois(50 * 10, 20), nrow = 50)

  result <- ace(species, threshold = 10)
  expect_s3_class(result, "spacc_estimate")
  expect_equal(result$estimate, result$S_obs)
})


test_that("jackknife order 1 returns correct structure", {
  set.seed(42)
  species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)

  result <- jackknife(species, order = 1)

  expect_s3_class(result, "spacc_estimate")
  expect_equal(result$estimator, "jackknife1")
  expect_true(result$estimate >= result$S_obs)
})


test_that("jackknife order 2 returns correct structure", {
  set.seed(42)
  species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)

  result <- jackknife(species, order = 2)

  expect_s3_class(result, "spacc_estimate")
  expect_equal(result$estimator, "jackknife2")
  expect_true(result$estimate >= result$S_obs)
})


test_that("jackknife rejects invalid order", {
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)
  expect_error(jackknife(species, order = 3), "order must be 1 or 2")
})


test_that("bootstrap_richness returns correct structure", {
  set.seed(42)
  species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)

  result <- bootstrap_richness(species, n_boot = 50)

  expect_s3_class(result, "spacc_estimate")
  expect_equal(result$estimator, "bootstrap")
  expect_true(result$estimate >= result$S_obs)
})


test_that("single-site data works", {
  species <- matrix(c(5, 3, 0, 1, 2), nrow = 1)

  result <- chao1(species)
  expect_s3_class(result, "spacc_estimate")
  expect_true(result$estimate >= result$S_obs)
})


test_that("print.spacc_estimate works", {
  set.seed(42)
  species <- matrix(rpois(50 * 30, 2), nrow = 50)
  result <- chao1(species)

  expect_output(print(result), "Richness Estimator")
  expect_output(print(result), "Observed species")
})


test_that("summary.spacc_estimate returns data.frame", {
  set.seed(42)
  species <- matrix(rpois(50 * 30, 2), nrow = 50)
  result <- chao1(species)

  summ <- summary(result)
  expect_s3_class(summ, "data.frame")
  expect_true("estimate" %in% names(summ))
})


test_that("as.data.frame.spacc_estimate works", {
  set.seed(42)
  species <- matrix(rpois(50 * 30, 2), nrow = 50)
  result <- chao1(species)

  df <- as.data.frame(result)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 1)
})


test_that("plot.spacc_estimate works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  species <- matrix(rpois(50 * 30, 2), nrow = 50)
  result <- chao1(species)

  p <- plot(result)
  expect_s3_class(p, "ggplot")
})


test_that("CI lower bound >= S_obs", {
  set.seed(42)
  species <- matrix(rpois(50 * 30, 2), nrow = 50)

  r1 <- chao1(species)
  r2 <- chao2(species)
  r3 <- ace(species)
  r4 <- jackknife(species)
  r5 <- bootstrap_richness(species, n_boot = 50)

  expect_true(r1$lower >= r1$S_obs - 1e-10)
  expect_true(r2$lower >= r2$S_obs - 1e-10)
  expect_true(r3$lower >= r3$S_obs - 1e-10)
  expect_true(r4$lower >= r4$S_obs - 1e-10)
  expect_true(r5$lower >= r5$S_obs - 1e-10)
})
