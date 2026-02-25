test_that("completenessProfile returns correct structure", {
  set.seed(42)
  species <- matrix(rpois(50 * 30, 2), nrow = 50)

  comp <- completenessProfile(species)

  expect_s3_class(comp, "spacc_completeness")
  expect_equal(length(comp$completeness), length(comp$q))
  expect_equal(length(comp$observed), length(comp$q))
  expect_equal(length(comp$estimated), length(comp$q))
  expect_null(comp$per_site)
  expect_null(comp$coords)
  expect_equal(comp$n_sites, 50)
  expect_equal(comp$n_species, 30)
})


test_that("completeness values are between 0 and 1", {
  set.seed(42)
  species <- matrix(rpois(50 * 30, 2), nrow = 50)

  comp <- completenessProfile(species, q = c(0, 1, 2))

  expect_true(all(comp$completeness >= 0))
  expect_true(all(comp$completeness <= 1))
})


test_that("completeness increases with q", {
  set.seed(42)
  species <- matrix(rpois(100 * 50, 2), nrow = 100)

  comp <- completenessProfile(species, q = c(0, 1, 2))

  # Generally completeness increases with q (dominant species detected early)
  expect_true(comp$completeness[3] >= comp$completeness[1] - 0.1)
})


test_that("estimated >= observed for all q", {
  set.seed(42)
  species <- matrix(rpois(50 * 30, 2), nrow = 50)

  comp <- completenessProfile(species, q = c(0, 1, 2))

  # Asymptotic estimates should be >= observed
  expect_true(all(comp$estimated >= comp$observed - 0.01))
})


test_that("per-site completeness computed with coords", {
  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 15, 2), nrow = 20)

  comp <- completenessProfile(species, q = c(0, 1, 2), coords = coords)

  expect_false(is.null(comp$per_site))
  expect_equal(nrow(comp$per_site), 20)
  expect_equal(ncol(comp$per_site), 3)
  expect_true(all(comp$per_site >= 0, na.rm = TRUE))
  expect_true(all(comp$per_site <= 1, na.rm = TRUE))
})


test_that("print.spacc_completeness works", {
  set.seed(42)
  species <- matrix(rpois(50 * 30, 2), nrow = 50)
  comp <- completenessProfile(species, q = c(0, 1, 2))

  expect_output(print(comp), "Sample Completeness")
})


test_that("summary returns data.frame", {
  set.seed(42)
  species <- matrix(rpois(50 * 30, 2), nrow = 50)
  comp <- completenessProfile(species, q = c(0, 1, 2))

  summ <- summary(comp)
  expect_s3_class(summ, "data.frame")
  expect_true("completeness" %in% names(summ))
})


test_that("as.data.frame works", {
  set.seed(42)
  species <- matrix(rpois(50 * 30, 2), nrow = 50)
  comp <- completenessProfile(species, q = c(0, 1, 2))

  df <- as.data.frame(comp)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 3)
})


test_that("plot profile works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  species <- matrix(rpois(50 * 30, 2), nrow = 50)
  comp <- completenessProfile(species, q = c(0, 1, 2))

  p <- plot(comp)
  expect_s3_class(p, "ggplot")
})


test_that("plot comparison works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  species <- matrix(rpois(50 * 30, 2), nrow = 50)
  comp <- completenessProfile(species, q = c(0, 1, 2))

  p <- plot(comp, type = "comparison")
  expect_s3_class(p, "ggplot")
})


test_that("plot map works with coords", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 15, 2), nrow = 20)
  comp <- completenessProfile(species, q = c(0, 1, 2), coords = coords)

  p <- plot(comp, type = "map")
  expect_s3_class(p, "ggplot")
})


test_that("plot map errors without coords", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  species <- matrix(rpois(50 * 30, 2), nrow = 50)
  comp <- completenessProfile(species, q = c(0, 1, 2))

  expect_error(plot(comp, type = "map"), "per-site")
})


test_that("as_sf works with coords", {
  skip_if_not_installed("sf")

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 15, 2), nrow = 20)
  comp <- completenessProfile(species, q = c(0, 1, 2), coords = coords)

  sf_obj <- as_sf(comp)
  expect_s3_class(sf_obj, "sf")
  expect_true("completeness_mean" %in% names(sf_obj))
})


test_that("as_sf errors without coords", {
  set.seed(42)
  species <- matrix(rpois(50 * 30, 2), nrow = 50)
  comp <- completenessProfile(species, q = c(0, 1, 2))

  expect_error(as_sf(comp), "per-site")
})


test_that("edge case: all species abundant", {
  set.seed(42)
  species <- matrix(rpois(50 * 10, 20), nrow = 50)

  comp <- completenessProfile(species, q = c(0, 1, 2))

  # High abundance = high completeness
  expect_true(all(comp$completeness > 0.5))
})
