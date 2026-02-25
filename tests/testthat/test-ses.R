test_that("ses returns correct structure for spacc", {
  set.seed(42)
  coords <- data.frame(x = runif(10), y = runif(10))
  species <- matrix(rbinom(10 * 8, 1, 0.4), nrow = 10)

  sac <- spacc(species, coords, n_seeds = 5, progress = FALSE)
  result <- ses(sac, species, n_perm = 9, progress = FALSE)

  expect_s3_class(result, "spacc_ses")
  expect_true(is.numeric(result$ses))
  expect_true(is.numeric(result$p_value))
  expect_true(is.numeric(result$null_mean))
  expect_true(is.numeric(result$null_sd))
  expect_equal(result$null_model, "frequency")
  expect_equal(result$input_class, "spacc")
})


test_that("p-values in [0, 1]", {
  set.seed(42)
  coords <- data.frame(x = runif(10), y = runif(10))
  species <- matrix(rbinom(10 * 8, 1, 0.4), nrow = 10)

  sac <- spacc(species, coords, n_seeds = 5, progress = FALSE)
  result <- ses(sac, species, n_perm = 9, progress = FALSE)

  expect_true(all(result$p_value >= 0))
  expect_true(all(result$p_value <= 1))
})


test_that("null models preserve row marginals (richness)", {
  set.seed(42)
  m <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)
  row_sums_orig <- rowSums(m)

  null_m <- spacc:::randomize_matrix(m, "richness")
  expect_equal(rowSums(null_m), row_sums_orig)
})


test_that("null models preserve column marginals (frequency)", {
  set.seed(42)
  m <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)
  col_sums_orig <- colSums(m)

  null_m <- spacc:::randomize_matrix(m, "frequency")
  expect_equal(colSums(null_m), col_sums_orig)
})


test_that("independent swap preserves both marginals", {
  set.seed(42)
  m <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)
  row_sums_orig <- rowSums(m)
  col_sums_orig <- colSums(m)

  null_m <- spacc:::randomize_matrix(m, "both")
  expect_equal(rowSums(null_m), row_sums_orig)
  expect_equal(colSums(null_m), col_sums_orig)
})


test_that("curveball preserves both marginals", {
  set.seed(42)
  m <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)
  row_sums_orig <- rowSums(m)
  col_sums_orig <- colSums(m)

  null_m <- spacc:::randomize_matrix(m, "curveball")
  expect_equal(rowSums(null_m), row_sums_orig)
  expect_equal(colSums(null_m), col_sums_orig)
})


test_that("ses works with spacc_alpha", {
  set.seed(42)
  coords <- data.frame(x = runif(10), y = runif(10))
  species <- matrix(rpois(10 * 8, 3), nrow = 10)

  alpha <- alphaDiversity(species, q = c(0, 1), coords = coords)
  result <- ses(alpha, species, n_perm = 9, progress = FALSE)

  expect_s3_class(result, "spacc_ses")
  expect_equal(result$input_class, "spacc_alpha")
})


test_that("ses rejects unsupported class", {
  x <- list()
  class(x) <- "spacc_fake"
  species <- matrix(1, 5, 5)

  expect_error(ses(x, species), "does not support class")
})


test_that("print.spacc_ses works", {
  set.seed(42)
  coords <- data.frame(x = runif(10), y = runif(10))
  species <- matrix(rbinom(10 * 8, 1, 0.4), nrow = 10)

  sac <- spacc(species, coords, n_seeds = 5, progress = FALSE)
  result <- ses(sac, species, n_perm = 9, progress = FALSE)

  expect_output(print(result), "Standardized Effect Size")
  expect_output(print(result), "Null model")
})


test_that("summary.spacc_ses returns data.frame", {
  set.seed(42)
  coords <- data.frame(x = runif(10), y = runif(10))
  species <- matrix(rbinom(10 * 8, 1, 0.4), nrow = 10)

  sac <- spacc(species, coords, n_seeds = 5, progress = FALSE)
  result <- ses(sac, species, n_perm = 9, progress = FALSE)

  summ <- summary(result)
  expect_s3_class(summ, "data.frame")
  expect_true("ses" %in% names(summ))
  expect_true("p_value" %in% names(summ))
})


test_that("as.data.frame.spacc_ses works", {
  set.seed(42)
  coords <- data.frame(x = runif(10), y = runif(10))
  species <- matrix(rbinom(10 * 8, 1, 0.4), nrow = 10)

  sac <- spacc(species, coords, n_seeds = 5, progress = FALSE)
  result <- ses(sac, species, n_perm = 9, progress = FALSE)

  df <- as.data.frame(result)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), length(result$ses))
})


test_that("plot.spacc_ses curve works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(10), y = runif(10))
  species <- matrix(rbinom(10 * 8, 1, 0.4), nrow = 10)

  sac <- spacc(species, coords, n_seeds = 5, progress = FALSE)
  result <- ses(sac, species, n_perm = 9, progress = FALSE)

  p <- plot(result)
  expect_s3_class(p, "ggplot")
})


test_that("plot.spacc_ses histogram works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(10), y = runif(10))
  species <- matrix(rbinom(10 * 8, 1, 0.4), nrow = 10)

  sac <- spacc(species, coords, n_seeds = 5, progress = FALSE)
  result <- ses(sac, species, n_perm = 9, progress = FALSE)

  p <- plot(result, type = "histogram")
  expect_s3_class(p, "ggplot")
})


test_that("ses seed argument gives reproducible results", {
  coords <- data.frame(x = runif(10), y = runif(10))
  species <- matrix(rbinom(10 * 8, 1, 0.4), nrow = 10)

  set.seed(42)
  sac <- spacc(species, coords, n_seeds = 5, progress = FALSE)

  r1 <- ses(sac, species, n_perm = 9, seed = 123, progress = FALSE)
  r2 <- ses(sac, species, n_perm = 9, seed = 123, progress = FALSE)

  expect_equal(r1$ses, r2$ses)
})


test_that("ses stores coords", {
  set.seed(42)
  coords <- data.frame(x = runif(10), y = runif(10))
  species <- matrix(rbinom(10 * 8, 1, 0.4), nrow = 10)

  sac <- spacc(species, coords, n_seeds = 5, progress = FALSE)
  result <- ses(sac, species, n_perm = 9, progress = FALSE)

  expect_false(is.null(result$coords))
  expect_equal(nrow(result$coords), 10)
})


test_that("plot.spacc_ses map works for per-site input", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(10), y = runif(10))
  species <- matrix(rpois(10 * 8, 3), nrow = 10)

  alpha <- alphaDiversity(species, q = c(0, 1), coords = coords)
  result <- ses(alpha, species, n_perm = 9, progress = FALSE)

  p <- plot(result, type = "map")
  expect_s3_class(p, "ggplot")

  p2 <- plot(result, type = "map", metric_map = "p_value")
  expect_s3_class(p2, "ggplot")
})


test_that("plot.spacc_ses map errors without coords", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  species <- matrix(rpois(10 * 8, 3), nrow = 10)

  alpha <- alphaDiversity(species, q = c(0, 1))
  # alpha without coords => plain data.frame, not spacc_alpha
  # Manually create a spacc_ses without coords for testing
  result <- structure(
    list(ses = rnorm(10), p_value = runif(10), coords = NULL,
         metric = "test", null_model = "frequency", n_valid = 9,
         input_class = "spacc_alpha"),
    class = "spacc_ses"
  )

  expect_error(plot(result, type = "map"), "No coordinates")
})


test_that("as_sf.spacc_ses works for per-site input", {
  skip_if_not_installed("sf")

  set.seed(42)
  coords <- data.frame(x = runif(10), y = runif(10))
  species <- matrix(rpois(10 * 8, 3), nrow = 10)

  alpha <- alphaDiversity(species, q = c(0, 1), coords = coords)
  result <- ses(alpha, species, n_perm = 9, progress = FALSE)

  sf_obj <- as_sf(result)
  expect_s3_class(sf_obj, "sf")
  expect_true("ses" %in% names(sf_obj))
  expect_true("p_value" %in% names(sf_obj))
})


test_that("as_sf.spacc_ses errors for curve-based input", {
  set.seed(42)
  coords <- data.frame(x = runif(10), y = runif(10))
  species <- matrix(rbinom(10 * 8, 1, 0.4), nrow = 10)

  sac <- spacc(species, coords, n_seeds = 5, progress = FALSE)
  result <- ses(sac, species, n_perm = 9, progress = FALSE)

  expect_error(as_sf(result), "per-site")
})
