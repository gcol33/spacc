# Tests for error handling and input validation
# Ensures proper error messages for invalid inputs

# =============================================================================
# SPACC INPUT VALIDATION
# =============================================================================

test_that("spacc rejects non-numeric species matrix", {
  coords <- data.frame(x = c(0, 1), y = c(0, 1))
  species <- matrix(c("a", "b", "c", "d"), nrow = 2)

  expect_error(
    spacc(species, coords, parallel = FALSE, progress = FALSE),
    "numeric"
  )
})


test_that("spacc rejects mismatched dimensions", {
  coords <- data.frame(x = c(0, 1, 2), y = c(0, 1, 2))  # 3 sites

species <- matrix(1:4, nrow = 2)  # 2 sites

  expect_error(
    spacc(species, coords, parallel = FALSE, progress = FALSE),
    "rows|sites|dimension"
  )
})


test_that("spacc rejects invalid method", {
  coords <- data.frame(x = c(0, 1), y = c(0, 1))
  species <- matrix(c(1, 0, 1, 1), nrow = 2)

  expect_error(
    spacc(species, coords, method = "invalid_method", parallel = FALSE, progress = FALSE),
    "should be one of"
  )
})


test_that("spacc rejects negative n_seeds", {
  coords <- data.frame(x = c(0, 1), y = c(0, 1))
  species <- matrix(c(1, 0, 1, 1), nrow = 2)

  expect_error(
    spacc(species, coords, n_seeds = -1, parallel = FALSE, progress = FALSE),
    "n_seeds|positive"
  )
})


# =============================================================================
# SPACCHILL INPUT VALIDATION
# =============================================================================

test_that("spaccHill rejects invalid q values", {
  coords <- data.frame(x = runif(10), y = runif(10))
  species <- matrix(rpois(10 * 5, 3), nrow = 10)

  expect_error(
    spaccHill(species, coords, q = -1, progress = FALSE),
    "q|negative"
  )
})


# =============================================================================
# SPACCBETA INPUT VALIDATION
# =============================================================================

test_that("spaccBeta rejects invalid index", {
  coords <- data.frame(x = runif(10), y = runif(10))
  species <- matrix(rbinom(10 * 5, 1, 0.5), nrow = 10)

  expect_error(
    spaccBeta(species, coords, index = "invalid", parallel = FALSE, progress = FALSE),
    "should be one of"
  )
})


# =============================================================================
# EXTRAPOLATE INPUT VALIDATION
# =============================================================================

test_that("extrapolate rejects invalid model", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(10), y = runif(10))
  species <- matrix(rbinom(10 * 5, 1, 0.5), nrow = 10)

  sac <- spacc(species, coords, n_seeds = 3, parallel = FALSE, progress = FALSE)

  expect_error(
    extrapolate(sac, model = "invalid_model"),
    "should be one of"
  )
})


# =============================================================================
# COMPARE INPUT VALIDATION
# =============================================================================

test_that("compare rejects non-spacc objects", {
  expect_error(
    spacc::compare(list(a = 1), list(b = 2)),
    "spacc"
  )
})
