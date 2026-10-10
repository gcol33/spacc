# Accumulation rate curve estimator (Shestopaloff 2025)

# One row whose species have the frequency counts y = (y_1..y_20, y_21+).
arc_fixture <- function(y) matrix(rep(seq_along(y), y), nrow = 1)

test_that("arc reproduces the author's reference code", {
  # reference values from github.com/kshestop/ARC (point_est.r, K = 20,
  # portions 0-0.30, 1000 points) on a sparse Poisson-Gamma sample and a
  # Steam catalogue snapshot, and on a Gamma(1, 1) sample where the weights
  # end at the log-logistic minimum of Eq. 8 rather than the global one
  sparse <- c(229, 69, 28, 7, 3, 1, rep(0, 15))
  steam <- c(1064, 621, 396, 290, 193, 161, 151, 106, 99, 81, 90, 69, 59, 54,
             55, 44, 53, 43, 38, 41, 1536)
  intermediate <- c(249, 129, 66, 30, 19, 9, 7, 3, 1, 0, 1, rep(0, 10))
  expect_equal(arc(arc_fixture(sparse), n_boot = 0)$estimate, 1076.444973, tolerance = 1e-7)
  expect_equal(arc(arc_fixture(steam), n_boot = 0)$estimate, 7765.742462, tolerance = 1e-7)
  expect_equal(arc(arc_fixture(intermediate), n_boot = 0)$estimate, 1013.949861, tolerance = 1e-7)
})

test_that("arc returns a spacc_estimate", {
  x <- arc_fixture(c(229, 69, 28, 7, 3, 1, rep(0, 15)))
  e <- arc(x, n_boot = 0)
  expect_s3_class(e, "spacc_estimate")
  expect_equal(e$estimator, "ARC")
  expect_equal(e$S_obs, 337)
  expect_true(is.na(e$lower) && is.na(e$upper) && is.na(e$se))
  expect_named(e$details$undetected, c("lognormal", "loglogistic", "weibull"))
  expect_true(all(e$details$weights >= 0 & e$details$weights <= 1))
  expect_equal(e$estimate, 337 + sum(e$details$weights * e$details$undetected),
               tolerance = 0.05)
  expect_output(print(e), "ARC")
})

test_that("incidence data use the number of occupied sites", {
  set.seed(4)
  x <- matrix(rpois(60 * 200, rep(rgamma(200, 0.3, 4), each = 60)), nrow = 60)
  occupied <- colSums(x > 0)
  expect_equal(arc(x, data = "incidence", n_boot = 0)$estimate,
               arc(matrix(occupied, nrow = 1), n_boot = 0)$estimate)
})

test_that("the bootstrap interval is reproducible and brackets the estimate", {
  skip_on_cran()
  x <- arc_fixture(c(229, 69, 28, 7, 3, 1, rep(0, 15)))
  set.seed(1)
  a <- arc(x, n_boot = 40, parallel = FALSE)
  set.seed(1)
  b <- arc(x, n_boot = 40, n_cores = 2)
  expect_equal(a$lower, b$lower)
  expect_equal(a$upper, b$upper)
  expect_true(a$lower < a$estimate && a$estimate < a$upper)
  expect_gt(a$se, 0)
  expect_length(a$details$boot, 40)
  expect_gt(a$details$penalty, 0)
})

test_that("adjust applies the smoothing adjustment", {
  x <- arc_fixture(c(229, 69, 28, 7, 3, 1, rep(0, 15)))
  expect_false(isTRUE(all.equal(arc(x, n_boot = 0, adjust = TRUE)$estimate,
                                arc(x, n_boot = 0)$estimate)))
})

test_that("arc validates its input", {
  expect_error(arc(matrix(c(30, 40, 50), nrow = 1), n_boot = 0), "too few species")
  expect_error(arc(matrix(1:5, nrow = 1), K = 1), "K must be")
})

test_that("arc recovers richness in sparse Poisson-Gamma samples", {
  skip_on_cran()
  set.seed(31)
  # 1000 species with Gamma(1, 2) rates (Shestopaloff 2025, scenario 1)
  r <- t(replicate(20, {
    x <- matrix(rpois(1000, rgamma(1000, 1, 2)), nrow = 1)
    c(arc = arc(x, n_boot = 0, parallel = FALSE)$estimate, iChao1 = iChao1(x)$estimate)
  }))
  rel <- (r - 1000) / 1000
  # the relative RMSE of ARC is about 0.18 here, so the median of 20 draws
  # has a standard error near 0.05
  expect_lt(abs(median(rel[, "arc"])), 0.12)
  expect_lt(sqrt(mean(rel[, "arc"]^2)), sqrt(mean(rel[, "iChao1"]^2)))
})
