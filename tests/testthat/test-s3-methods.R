# Tests for S3 methods (as.data.frame, print, summary, plot, [, c)
# Verifies that all custom spacc classes have working S3 methods

# =============================================================================
# AS.DATA.FRAME METHODS
# =============================================================================

test_that("as.data.frame.spacc works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  result <- spacc(species, coords, n_seeds = 3, parallel = FALSE, progress = FALSE)
  df <- as.data.frame(result)

  expect_s3_class(df, "data.frame")
  expect_true("sites" %in% names(df))
  expect_true("mean" %in% names(df))
})


test_that("as.data.frame.spacc_fit works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  sac <- spacc(species, coords, n_seeds = 5, parallel = FALSE, progress = FALSE)
  fit <- extrapolate(sac)
  df <- as.data.frame(fit)

  expect_s3_class(df, "data.frame")
  expect_true("sites" %in% names(df))
})


test_that("as.data.frame.spacc_hill works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 3), nrow = 15)

  result <- spaccHill(species, coords, n_seeds = 3, q = c(0, 1), progress = FALSE)
  df <- as.data.frame(result)

  expect_s3_class(df, "data.frame")
  expect_true("q" %in% names(df))
})


test_that("as.data.frame.spacc_beta works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  result <- spaccBeta(species, coords, n_seeds = 3, parallel = FALSE, progress = FALSE)
  df <- as.data.frame(result)

  expect_s3_class(df, "data.frame")
  expect_true("beta_total" %in% names(df))
})


# =============================================================================
# PRINT AND SUMMARY METHODS
# =============================================================================

test_that("print.spacc works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  result <- spacc(species, coords, n_seeds = 3, parallel = FALSE, progress = FALSE)

  expect_output(print(result), "spacc:")
})


test_that("summary.spacc works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  result <- spacc(species, coords, n_seeds = 3, parallel = FALSE, progress = FALSE)
  summ <- summary(result)

  expect_s3_class(summ, "summary.spacc")
  expect_true(!is.null(summ$mean))
})


# =============================================================================
# SUBSETTING AND COMBINING
# =============================================================================

test_that("[ subsetting works for spacc", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  result <- spacc(species, coords, n_seeds = 5, parallel = FALSE, progress = FALSE)
  subset <- result[1:3]

  expect_s3_class(subset, "spacc")
  expect_equal(subset$n_seeds, 3)
})


test_that("c.spacc combines objects", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  r1 <- spacc(species, coords, n_seeds = 3, parallel = FALSE, progress = FALSE)
  r2 <- spacc(species, coords, n_seeds = 3, parallel = FALSE, progress = FALSE)

  combined <- c(group1 = r1, group2 = r2)

  expect_s3_class(combined, "spacc")
  expect_true(spacc:::is_grouped(combined))
})
