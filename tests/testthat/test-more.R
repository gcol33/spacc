more_fixture <- function() {
  # 10 sites; frequencies: 4 uniques, 2 duplicates, 3 species at 5 sites
  freq <- c(1, 1, 1, 1, 2, 2, 5, 5, 5)
  m <- matrix(0L, 10, length(freq))
  for (j in seq_along(freq)) m[seq_len(freq[j]), j] <- 1L
  m
}

test_that("chao2 with n_units = NULL or Inf is the with-replacement estimator", {
  x <- more_fixture()
  a <- chao2(x)
  b <- chao2(x, n_units = Inf)
  expect_identical(a$estimate, b$estimate)
  expect_identical(a$estimator, "chao2")
  expect_equal(a$estimate, 9 + 9 / 10 * 16 / 4)
})

test_that("MoRE matches Chiu (2023) Eq. 5 by hand", {
  x <- more_fixture()
  w <- 10 / 40
  # Q1 = 4 > Q2 = 2: denominator Q2 + w (Q1 - Q2)
  expected <- 9 + 9 / 10 * (1 - w) * 16 / (2 + w * 2)
  e <- chao2(x, n_units = 40)
  expect_identical(e$estimator, "MoRE")
  expect_equal(e$estimate, expected)
  expect_true(e$lower > e$S_obs && e$upper > e$estimate)
})

test_that("MoRE uses (A)^+ = 1 when Q1 <= Q2 and the Q2 = 0 branch", {
  # Q1 = 1, Q2 = 2
  x <- matrix(0L, 6, 4)
  x[1, 1] <- 1L; x[1:2, 2] <- 1L; x[3:4, 3] <- 1L; x[, 4] <- 1L
  w <- 6 / 12
  expect_equal(chao2(x, n_units = 12)$estimate, 4 + 5 / 6 * (1 - w) * 1 / (2 + w))

  # Q1 = 3, Q2 = 0
  y <- matrix(0L, 6, 4)
  y[1, 1] <- 1L; y[2, 2] <- 1L; y[3, 3] <- 1L; y[, 4] <- 1L
  expect_equal(chao2(y, n_units = 12)$estimate, 4 + 5 / 6 * (1 - w) / w * 2)
})

test_that("MoRE equals observed richness when every unit is sampled", {
  x <- more_fixture()
  e <- chao2(x, n_units = nrow(x))
  expect_equal(e$estimate, e$S_obs)
  expect_equal(e$se, 0, tolerance = 1e-12)
})

test_that("MoRE decreases towards S_obs as the sampling fraction grows", {
  x <- more_fixture()
  est <- vapply(c(1e4, 100, 40, 20, 12), function(n) chao2(x, n_units = n)$estimate,
                numeric(1))
  expect_true(all(diff(est) < 0))
})

test_that("n_units smaller than the number of sites is an error", {
  expect_error(chao2(more_fixture(), n_units = 5), "at least the number of sampled units")
})

test_that("extrapolate(n_units = ) compares the asymptote with MoRE", {
  set.seed(2)
  coords <- data.frame(x = runif(40), y = runif(40))
  species <- matrix(rbinom(40 * 25, 1, 0.2), nrow = 40)
  sac <- spacc(species, coords, n_seeds = 5, progress = FALSE)
  fit <- suppressWarnings(extrapolate(sac, interval = "none", n_units = 50))
  expect_equal(fit$compare$MoRE, chao2(species, n_units = 50)$estimate)
  expect_output(print(fit), "MoRE")
  fit0 <- suppressWarnings(extrapolate(sac, interval = "none"))
  expect_null(fit0$compare$MoRE)
})
