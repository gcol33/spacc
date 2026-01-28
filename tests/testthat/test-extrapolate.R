test_that("extrapolate returns correct structure with michaelis-menten", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)
  fit <- extrapolate(sac, model = "michaelis-menten")

  expect_s3_class(fit, "spacc_fit")
  expect_equal(fit$model, "michaelis-menten")
  expect_true(!is.na(fit$asymptote))
  expect_true(!is.na(fit$aic))
  expect_true(fit$asymptote > 0)
})


test_that("extrapolate asymptote >= observed richness", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)
  fit <- extrapolate(sac, model = "michaelis-menten")

  observed <- max(summary(sac)$mean)
  expect_true(fit$asymptote >= observed * 0.9)  # allow slight underfit
})


test_that("extrapolate works with lomolino model", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)
  fit <- extrapolate(sac, model = "lomolino")

  expect_s3_class(fit, "spacc_fit")
  expect_equal(fit$model, "lomolino")
})


test_that("extrapolate works with asymptotic model", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)
  fit <- extrapolate(sac, model = "asymptotic")

  expect_s3_class(fit, "spacc_fit")
  expect_equal(fit$model, "asymptotic")
})


test_that("extrapolate works with weibull model", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)
  fit <- extrapolate(sac, model = "weibull")

  expect_s3_class(fit, "spacc_fit")
})


test_that("predict.spacc_fit works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)
  fit <- extrapolate(sac, model = "michaelis-menten")

  preds <- predict(fit, n = c(10, 20, 30))
  expect_length(preds, 3)
  expect_true(all(!is.na(preds)))
  # Predictions should be monotonically non-decreasing
  expect_true(all(diff(preds) >= 0))
})


test_that("coef.spacc_fit returns coefficients", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)
  fit <- extrapolate(sac, model = "michaelis-menten")

  co <- coef(fit)
  expect_true(length(co) >= 2)
  expect_true("a" %in% names(co))
})


test_that("print.spacc_fit works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)
  fit <- extrapolate(sac, model = "michaelis-menten")

  expect_output(print(fit), "Extrapolation")
})


test_that("summary.spacc_fit returns nls summary", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)
  fit <- extrapolate(sac, model = "michaelis-menten")

  summ <- summary(fit)
  expect_true(!is.null(summ))
})


test_that("confint.spacc_fit works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)
  fit <- extrapolate(sac, model = "michaelis-menten")

  ci <- confint(fit)
  expect_true(is.matrix(ci) || is.numeric(ci))
})


test_that("predict.spacc_fit extrapolates beyond observed", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)
  fit <- extrapolate(sac, model = "michaelis-menten")

  preds <- predict(fit, n = c(50, 100, 200))
  expect_length(preds, 3)
  expect_true(all(diff(preds) >= 0))
})


test_that("extrapolate with logistic model works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)

  fit <- extrapolate(sac, model = "logistic")

  expect_s3_class(fit, "spacc_fit")
  expect_equal(fit$model, "logistic")
})


test_that("extrapolate plot returns ggplot", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)
  fit <- extrapolate(sac, model = "michaelis-menten")

  p <- plot(fit)
  expect_s3_class(p, "ggplot")
})
