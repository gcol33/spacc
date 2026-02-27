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


# --- compareModels tests ---

test_that("compareModels returns correct structure", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)

  cm <- compareModels(sac, models = c("michaelis-menten", "asymptotic", "logistic"))

  expect_s3_class(cm, "spacc_model_compare")
  expect_true(is.data.frame(cm$table))
  expect_true("AIC" %in% names(cm$table))
  expect_true("delta_AIC" %in% names(cm$table))
  expect_true("AIC_weight" %in% names(cm$table))
  expect_true("BIC" %in% names(cm$table))
  expect_true("converged" %in% names(cm$table))
  expect_true(is.list(cm$fits))
  expect_true(!is.na(cm$best_model))
})


test_that("compareModels table is sorted by AIC", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)

  cm <- compareModels(sac, models = c("michaelis-menten", "asymptotic"))
  aic_vals <- cm$table$AIC[cm$table$converged]
  expect_true(all(diff(aic_vals) >= 0))
})


test_that("compareModels Akaike weights sum to 1", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)

  cm <- compareModels(sac, models = c("michaelis-menten", "asymptotic", "logistic"))

  weights <- cm$table$AIC_weight[cm$table$converged]
  expect_equal(sum(weights), 1.0, tolerance = 1e-6)
})


test_that("compareModels S3 methods work", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)

  cm <- compareModels(sac, models = c("michaelis-menten", "asymptotic"))

  # print
  expect_output(print(cm), "SAR Model Comparison")

  # coef
  co <- coef(cm)
  expect_true(is.numeric(co))
  expect_true("a" %in% names(co))

  # predict
  preds <- predict(cm, n = c(10, 20, 30))
  expect_length(preds, 3)
  expect_true(all(!is.na(preds)))

  # as.data.frame
  df <- as.data.frame(cm)
  expect_true(is.data.frame(df))
  expect_equal(nrow(df), 2)
})


test_that("compareModels plot returns ggplot", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)

  cm <- compareModels(sac, models = c("michaelis-menten", "asymptotic"))

  p <- plot(cm)
  expect_s3_class(p, "ggplot")
})


test_that("compareModels best_model has delta_AIC = 0", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)

  cm <- compareModels(sac, models = c("michaelis-menten", "asymptotic"))

  best_row <- cm$table[cm$table$model == cm$best_model, ]
  expect_equal(best_row$delta_AIC, 0)
})
