# ---- SAR model registry --------------------------------------------------
# Each entry defines an asymptotic species-accumulation model S(x). `fn` is the
# mean-value function; the asymptote is always the parameter named "a". `start`
# builds nls starting values from the data summary (y_max, x_half); "port"
# entries add box constraints via `lower`/`upper`. Both fitting (nls) and
# prediction read this single definition, so adding a model is one entry.
.sar_registry <- list(
  "michaelis-menten" = list(
    fn = function(x, a, b) a * x / (b + x),
    start = function(ym, xh) list(a = ym * 1.2, b = xh),
    algorithm = "default",
    control = stats::nls.control()
  ),
  "lomolino" = list(
    fn = function(x, a, b, c) a / (1 + b^(log(c / x))),
    start = function(ym, xh) list(a = ym * 1.2, b = 2, c = xh),
    algorithm = "default",
    control = stats::nls.control(maxiter = 200)
  ),
  "asymptotic" = list(
    fn = function(x, a, b) a * (1 - exp(-b * x)),
    start = function(ym, xh) list(a = ym * 1.2, b = 0.01),
    algorithm = "default",
    control = stats::nls.control()
  ),
  "weibull" = list(
    fn = function(x, a, b, c) a * (1 - exp(-(x / b)^c)),
    start = function(ym, xh) list(a = ym * 1.2, b = xh, c = 1),
    algorithm = "default",
    control = stats::nls.control(maxiter = 200)
  ),
  "logistic" = list(
    fn = function(x, a, b, c) a / (1 + exp(-b * (x - c))),
    start = function(ym, xh) list(a = ym * 1.2, b = 0.1, c = xh),
    algorithm = "default",
    control = stats::nls.control()
  ),
  "evt" = list(
    # EVT-inspired triphasic SAR (Borda-de-Agua et al. 2025): mixture of two
    # Weibull components capturing the small/intermediate/large-scale phases.
    fn = function(x, a, w, b1, c1, b2, c2)
      a * (w * (1 - exp(-(x / b1)^c1)) + (1 - w) * (1 - exp(-(x / b2)^c2))),
    start = function(ym, xh) list(a = ym * 1.3, w = 0.7,
                                  b1 = max(1, xh * 0.3), c1 = 1.5,
                                  b2 = max(1, xh * 2), c2 = 0.5),
    lower = function(ym, xh) c(ym * 0.5, 0.01, 0.1, 0.1, 0.1, 0.1),
    upper = function(ym, xh) c(ym * 5, 0.99, xh * 10, 10, xh * 20, 10),
    algorithm = "port",
    control = stats::nls.control(maxiter = 500, warnOnly = TRUE)
  )
)

.sar_params <- function(model) {
  setdiff(names(formals(.sar_registry[[model]]$fn)), "x")
}

# Evaluate a fitted model at arbitrary effort from its coefficients. Uses the
# same `fn` as the nls formula, so prediction and fitting never diverge.
.sar_eval <- function(model, coefs, x) {
  fn <- .sar_registry[[model]]$fn
  p <- .sar_params(model)
  do.call(fn, c(list(x = x), as.list(coefs[p])))
}

# Fit one SAR model to a data frame with columns x, y. Returns the nls fit or
# NULL on failure. `quiet = TRUE` (used inside the bootstrap) suppresses both
# the warning and nls's own convergence chatter.
.fit_sar <- function(df, model, quiet = TRUE, ...) {
  if (model == "evt" && nrow(df) < 10) {
    if (!quiet) warning("EVT model requires >= 10 data points; falling back to Weibull",
                        call. = FALSE)
    return(.fit_sar(df, "weibull", quiet = quiet, ...))
  }

  spec <- .sar_registry[[model]]
  fn <- spec$fn
  ym <- max(df$y)
  xh <- df$x[which.min(abs(df$y - ym / 2))]
  if (!length(xh) || is.na(xh) || xh <= 0) xh <- max(stats::median(df$x), 1)

  p <- setdiff(names(formals(fn)), "x")
  rhs <- as.call(c(list(quote(fn), quote(x)), lapply(p, as.name)))
  form <- stats::as.formula(call("~", quote(y), rhs))
  environment(form) <- list2env(list(fn = fn), parent = environment())

  args <- list(formula = form, data = df, start = spec$start(ym, xh),
               control = spec$control)
  if (identical(spec$algorithm, "port")) {
    args$algorithm <- "port"
    args$lower <- spec$lower(ym, xh)
    args$upper <- spec$upper(ym, xh)
  }
  args <- utils::modifyList(args, list(...))  # user-supplied nls args win

  tryCatch(
    suppressWarnings(do.call(stats::nls, args)),
    error = function(e) {
      if (!quiet) warning("Model fitting failed: ", e$message, call. = FALSE)
      NULL
    }
  )
}

# Across-seed bootstrap: refit `model` to individual seed curves drawn with
# replacement, returning a matrix of coefficients (one converged fit per row).
# Resampling seeds (rather than re-averaging) keeps the spread from collapsing
# as the number of seeds grows -- it is the path/across-seed uncertainty, not
# the (shrinking) uncertainty of the mean curve.
.extrapolate_boot <- function(curves, model, R, ...) {
  S <- nrow(curves)
  n <- ncol(curves)
  if (is.null(S) || S < 2L) return(NULL)

  idx <- sample.int(S, size = R, replace = TRUE)
  coef_list <- vector("list", R)
  for (r in seq_len(R)) {
    dfr <- data.frame(x = seq_len(n), y = curves[idx[r], ])
    fr <- .fit_sar(dfr, model, quiet = TRUE, ...)
    if (!is.null(fr)) coef_list[[r]] <- stats::coef(fr)
  }
  coef_list <- Filter(Negate(is.null), coef_list)
  if (!length(coef_list)) return(NULL)

  coefs <- do.call(rbind, coef_list)
  list(coefs = coefs, n_ok = nrow(coefs), n_fail = R - nrow(coefs))
}


#' Extrapolate Total Species Richness
#'
#' Fit an asymptotic model to the spatial accumulation curve and estimate total
#' species richness beyond the observed sampling effort.
#'
#' @param object A `spacc` object.
#' @param model Character. Model to fit: `"michaelis-menten"` (default),
#'   `"lomolino"`, `"asymptotic"`, `"weibull"`, `"logistic"`, or `"evt"`
#'   (Extreme Value Theory, Borda-de-Agua et al. 2025).
#' @param interval Character. How to build the asymptote confidence interval:
#'   `"bootstrap"` (default) refits the model across resampled seed curves and
#'   takes percentile bounds; `"profile"` uses the (over-confident) `nls`
#'   profile interval on the mean-curve fit; `"none"` skips it.
#' @param R Integer. Number of bootstrap refits when `interval = "bootstrap"`.
#'   Default 200.
#' @param level Numeric. Confidence level for the interval. Default 0.95.
#' @param compare Logical. If `TRUE` (default) and the object carries incidence
#'   frequencies, compare the asymptote to the nonparametric `chao2()` /
#'   `iChao2()` estimates and flag large disagreement.
#' @param warn_ratio Numeric. Warn when the fitted asymptote exceeds the
#'   observed richness by more than this factor. Default 2. Set to `Inf` to
#'   silence.
#' @param ... Additional arguments passed to [stats::nls()].
#'
#' @return An object of class `spacc_fit` containing:
#'   \item{asymptote}{Estimated total species richness (asymptote of the model)}
#'   \item{asymptote_ci}{Confidence interval for the asymptote}
#'   \item{model}{Model name}
#'   \item{interval}{Interval method used}
#'   \item{fit}{The `nls` fit object}
#'   \item{aic}{AIC of the model}
#'   \item{gof}{List with residual `rmse` and relative `rmse_rel` over the
#'     observed range}
#'   \item{compare}{Nonparametric `chao2` / `iChao2` estimates, if available}
#'   \item{boot}{Bootstrap coefficient draws, if `interval = "bootstrap"`}
#'
#' @section Bias caveat:
#' The asymptote is a *parametric extrapolation* of the accumulation curve. On
#' clustered or strongly under-sampled presence-absence data it tends to
#' overestimate true richness, sometimes substantially, because the curve is far
#' from saturation. The bootstrap interval quantifies the uncertainty of the
#' fitted asymptote (curve-fit and across-seed variability) and is much wider
#' than the `nls` profile interval, but it is *not* a calibrated interval for
#' true total richness: it is centred on a possibly biased point estimate. For
#' calibrated total-richness estimates prefer the nonparametric estimators
#' [chao2()] / [iChao2()], which are unbiased on the same data. The printout
#' shows their values alongside the asymptote for comparison.
#'
#' @references
#' Lomolino, M.V. (2000). Ecology's most general, yet protean pattern: the
#' species-area relationship. Journal of Biogeography, 27, 17-26.
#'
#' Flather, C.H. (1996). Fitting species-accumulation functions and assessing
#' regional land use impacts on avian diversity. Journal of Biogeography,
#' 23, 155-168.
#'
#' Borda-de-Agua, L., Whittaker, R.J., Cardoso, P., et al. (2025). Extreme
#' value theory explains the topography and scaling of the species-area
#' relationship. Nature Communications, 16, 5346.
#'
#' @seealso [extrapolateArea()] for area-based extrapolation to a region larger
#'   than the one sampled; [chao2()], [iChao2()] for nonparametric richness.
#'
#' @examples
#' \donttest{
#' coords <- data.frame(x = runif(50), y = runif(50))
#' species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)
#' sac <- spacc(species, coords)
#' fit <- extrapolate(sac)
#' print(fit)
#' }
#'
#' @export
extrapolate <- function(object,
                        model = c("michaelis-menten", "lomolino", "asymptotic",
                                  "weibull", "logistic", "evt"),
                        interval = c("bootstrap", "profile", "none"),
                        R = 200L,
                        level = 0.95,
                        compare = TRUE,
                        warn_ratio = 2,
                        ...) {

  model <- match.arg(model)
  interval <- match.arg(interval)

  stopifnot("object must be a spacc object" = inherits(object, "spacc"))
  if (is_grouped(object)) {
    stop("extrapolate() does not support grouped spacc objects; extrapolate each group separately.",
         call. = FALSE)
  }

  summ <- summary(object)
  df <- data.frame(x = summ$sites, y = summ$mean)

  fit <- .fit_sar(df, model, quiet = FALSE, ...)

  if (is.null(fit)) {
    return(structure(
      list(asymptote = NA_real_, asymptote_ci = c(NA_real_, NA_real_),
           model = model, interval = interval, fit = NULL, aic = NA_real_,
           gof = NULL, compare = NULL, boot = NULL, level = level,
           observed = max(df$y), max_effort = max(df$x),
           data = df, spacc = object),
      class = "spacc_fit"
    ))
  }

  coefs <- stats::coef(fit)
  asymptote <- unname(coefs["a"])
  alpha <- (1 - level) / 2

  # ---- asymptote interval ----
  boot <- NULL
  if (interval == "bootstrap") {
    boot <- .extrapolate_boot(object$curves, model, R, ...)
    if (is.null(boot)) {
      warning("Bootstrap interval unavailable (need >= 2 seeds and some converged refits); ",
              "falling back to profile interval.", call. = FALSE)
      interval <- "profile"
    } else {
      ci <- unname(stats::quantile(boot$coefs[, "a"], c(alpha, 1 - alpha), na.rm = TRUE))
    }
  }
  if (interval == "profile") {
    ci <- tryCatch(suppressMessages(stats::confint(fit, parm = "a", level = level)),
                   error = function(e) c(NA_real_, NA_real_))
    ci <- unname(as.numeric(ci))
  }
  if (interval == "none") ci <- c(NA_real_, NA_real_)

  # ---- goodness of fit over the observed range ----
  resid <- df$y - .sar_eval(model, coefs, df$x)
  rmse <- sqrt(mean(resid^2))
  gof <- list(rmse = rmse, rmse_rel = rmse / mean(df$y), n_points = nrow(df))

  observed <- max(df$y)
  max_effort <- max(df$x)

  # ---- nonparametric comparison ----
  cmp <- NULL
  if (isTRUE(compare) && !is.null(object$incidence_freq)) {
    c2 <- .chao2_from_freq(object$incidence_freq, object$n_sites)
    ic2 <- .ichao2_from_freq(object$incidence_freq, object$n_sites)
    cmp <- list(chao2 = c2$estimate, iChao2 = ic2$estimate, S_obs = c2$S_obs)
  }

  # ---- guardrail warnings (issue #5) ----
  if (is.finite(warn_ratio) && asymptote > warn_ratio * observed) {
    warning(sprintf(
      "Asymptote (%.0f) exceeds observed richness (%.0f) by more than %.1fx; parametric extrapolation may be unreliable on under-sampled data. See ?extrapolate.",
      asymptote, observed, warn_ratio), call. = FALSE)
  }
  if (!is.null(cmp) && is.finite(cmp$chao2) && cmp$chao2 > 0) {
    rel <- asymptote / cmp$chao2
    if (rel > 1.5 || rel < 1 / 1.5) {
      warning(sprintf(
        "Asymptote (%.0f) disagrees with nonparametric chao2 (%.0f) by %.0f%%; prefer chao2()/iChao2() for total-richness intervals.",
        asymptote, cmp$chao2, 100 * abs(rel - 1)), call. = FALSE)
    }
  }

  structure(
    list(
      asymptote = asymptote,
      asymptote_ci = ci,
      model = model,
      interval = interval,
      fit = fit,
      aic = stats::AIC(fit),
      gof = gof,
      compare = cmp,
      boot = boot,
      level = level,
      observed = observed,
      max_effort = max_effort,
      data = df,
      spacc = object
    ),
    class = "spacc_fit"
  )
}


#' @export
print.spacc_fit <- function(x, ...) {
  cat("Extrapolation:", x$model, "\n")
  cat(strrep("-", 38), "\n")
  if (is.na(x$asymptote)) {
    cat("Model fitting failed\n")
    return(invisible(x))
  }

  cat(sprintf("Estimated asymptote: %.1f species\n", x$asymptote))
  if (all(is.finite(x$asymptote_ci))) {
    cat(sprintf("%.0f%% CI (%s):       %.1f - %.1f\n",
                100 * (x$level %||% 0.95), x$interval,
                x$asymptote_ci[1], x$asymptote_ci[2]))
  }
  cat(sprintf("Observed:            %.1f species (%.0f%% of estimated)\n",
              x$observed, 100 * x$observed / x$asymptote))
  cat(sprintf("AIC: %.1f", x$aic))
  if (!is.null(x$gof)) {
    cat(sprintf("   RMSE: %.2f (%.1f%% of mean)\n",
                x$gof$rmse, 100 * x$gof$rmse_rel))
  } else cat("\n")
  cat(sprintf("Reliable to ~%.0f sites (2.5x sampled effort of %d)\n",
              2.5 * x$max_effort, x$max_effort))

  if (!is.null(x$compare)) {
    cat(strrep("-", 38), "\n")
    cat(sprintf("Nonparametric:  chao2 = %.1f   iChao2 = %.1f\n",
                x$compare$chao2, x$compare$iChao2))
    cat("(prefer these for calibrated total-richness intervals; see ?extrapolate)\n")
  }
  invisible(x)
}


#' @export
summary.spacc_fit <- function(object, ...) {
  if (!is.null(object$fit)) {
    summary(object$fit)
  } else {
    print(object)
  }
}


#' @export
coef.spacc_fit <- function(object, ...) {
  if (!is.null(object$fit)) {
    stats::coef(object$fit)
  } else {
    NA
  }
}


#' Confidence interval for an extrapolated asymptote
#'
#' Returns the calibrated bootstrap interval stored on the fit by default. With
#' `method = "profile"` the (over-confident) `nls` profile interval is computed
#' instead. The interval is for the asymptote *parameter*; see the bias caveat
#' in [extrapolate()].
#'
#' @param object A `spacc_fit` object.
#' @param parm Ignored; the asymptote (`"a"`) is always returned.
#' @param level Confidence level. Recomputed from stored bootstrap draws.
#' @param method `"stored"` (default) or `"profile"`.
#' @param ... Passed to [stats::confint()] when `method = "profile"`.
#' @return A 1-row matrix with the lower and upper bound for `a`.
#' @export
confint.spacc_fit <- function(object, parm, level = 0.95, method = c("stored", "profile"), ...) {
  method <- match.arg(method)
  if (is.null(object$fit)) return(NA)

  if (method == "profile") {
    ci <- tryCatch(suppressMessages(stats::confint(object$fit, parm = "a", level = level, ...)),
                   error = function(e) c(NA_real_, NA_real_))
    ci <- as.numeric(ci)
  } else if (!is.null(object$boot)) {
    alpha <- (1 - level) / 2
    ci <- unname(stats::quantile(object$boot$coefs[, "a"], c(alpha, 1 - alpha), na.rm = TRUE))
  } else {
    ci <- as.numeric(object$asymptote_ci)
  }

  out <- matrix(ci, nrow = 1,
                dimnames = list("a", c(sprintf("%.1f %%", 100 * (1 - level) / 2),
                                       sprintf("%.1f %%", 100 * (1 - (1 - level) / 2)))))
  out
}


#' Predict richness at a given sampling effort
#'
#' @param object A `spacc_fit` object.
#' @param n Integer vector of site counts. Defaults to the observed range.
#' @param interval `"none"` (default) returns point predictions; `"bootstrap"`
#'   returns a data frame with percentile bounds from the stored bootstrap draws.
#' @param level Confidence level for the interval. Defaults to the fit's level.
#' @param ... Unused.
#' @return A numeric vector (point) or a data frame with `n`, `fit`, `lower`,
#'   `upper` (interval).
#' @export
predict.spacc_fit <- function(object, n = NULL, interval = c("none", "bootstrap"),
                              level = NULL, ...) {
  if (is.null(object$fit)) {
    return(NA)
  }
  interval <- match.arg(interval)
  if (is.null(n)) n <- object$data$x

  if (!is.null(object$max_effort) && any(n > 2.5 * object$max_effort)) {
    warning(sprintf(
      "Predicting beyond ~2.5x the sampled effort (max sampled = %d sites); extrapolation is unreliable this far out.",
      object$max_effort), call. = FALSE)
  }

  point <- .sar_eval(object$model, stats::coef(object$fit), n)

  if (interval == "none" || is.null(object$boot)) {
    return(point)
  }

  if (is.null(level)) level <- object$level %||% 0.95
  alpha <- (1 - level) / 2
  pm <- vapply(seq_len(nrow(object$boot$coefs)),
               function(r) .sar_eval(object$model, object$boot$coefs[r, ], n),
               numeric(length(n)))
  pm <- matrix(pm, nrow = length(n))
  data.frame(
    n = n,
    fit = point,
    lower = apply(pm, 1, stats::quantile, probs = alpha, na.rm = TRUE),
    upper = apply(pm, 1, stats::quantile, probs = 1 - alpha, na.rm = TRUE)
  )
}


#' @export
plot.spacc_fit <- function(x, extrapolate_to = NULL, interval = TRUE, ...) {
  check_suggests("ggplot2")

  df <- x$data

  if (is.null(extrapolate_to)) {
    extrapolate_to <- nrow(df) * 2
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(color = "#2E7D32", alpha = 0.5) +
    ggplot2::geom_line(color = "#2E7D32")

  if (!is.null(x$fit)) {
    pred_x <- seq(1, extrapolate_to, length.out = 200)

    if (isTRUE(interval) && !is.null(x$boot)) {
      pb <- suppressWarnings(predict(x, n = pred_x, interval = "bootstrap"))
      p <- p + ggplot2::geom_ribbon(
        data = pb,
        ggplot2::aes(x = n, ymin = lower, ymax = upper),
        inherit.aes = FALSE, fill = "#FF9800", alpha = 0.2
      )
    }

    pred_y <- suppressWarnings(predict(x, n = pred_x))
    pred_df <- data.frame(x = pred_x, y = pred_y)

    p <- p +
      ggplot2::geom_line(
        data = pred_df,
        ggplot2::aes(x = x, y = y),
        color = "#FF9800", linewidth = 1, linetype = "dashed"
      ) +
      ggplot2::geom_hline(
        yintercept = x$asymptote,
        linetype = "dotted", color = "#F44336"
      ) +
      ggplot2::annotate(
        "text",
        x = extrapolate_to * 0.9,
        y = x$asymptote * 1.02,
        label = sprintf("Asymptote: %.0f", x$asymptote),
        color = "#F44336", size = 3.5, hjust = 1
      )
  }

  p +
    ggplot2::labs(
      x = "Sites sampled",
      y = "Cumulative species",
      title = "Species Accumulation with Extrapolation",
      subtitle = sprintf("Model: %s, AIC: %.1f", x$model, x$aic)
    ) +
    spacc_theme()
}


#' Compare Multiple SAR Models
#'
#' Fit all (or a subset of) asymptotic species-area models and compare them
#' using AIC, BIC, delta-AIC, and Akaike weights.
#'
#' @param object A `spacc` object.
#' @param models Character vector of models to fit. Defaults to all six:
#'   `"michaelis-menten"`, `"lomolino"`, `"asymptotic"`, `"weibull"`,
#'   `"logistic"`, `"evt"`.
#' @param ... Additional arguments passed to [extrapolate()].
#'
#' @return An object of class `spacc_model_compare` containing:
#'   \item{table}{Data frame with model comparison statistics}
#'   \item{fits}{Named list of `spacc_fit` objects}
#'   \item{best_model}{Name of the best model by AIC}
#'   \item{data}{Mean-curve data frame used for fitting}
#'   \item{spacc}{Original spacc object}
#'
#' @examples
#' \donttest{
#' coords <- data.frame(x = runif(50), y = runif(50))
#' species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)
#' sac <- spacc(species, coords)
#' cm <- compareModels(sac)
#' print(cm)
#' }
#'
#' @export
compareModels <- function(object,
                          models = c("michaelis-menten", "lomolino",
                                     "asymptotic", "weibull",
                                     "logistic", "evt"),
                          ...) {

  stopifnot("object must be a spacc object" = inherits(object, "spacc"))
  models <- match.arg(models, several.ok = TRUE)

  # Model comparison is about AIC, not intervals: skip the bootstrap and
  # guardrail warnings here (callers can request them per-model via extrapolate).
  dots <- list(...)
  if (is.null(dots$interval)) dots$interval <- "none"
  if (is.null(dots$compare)) dots$compare <- FALSE
  if (is.null(dots$warn_ratio)) dots$warn_ratio <- Inf

  fits <- stats::setNames(
    lapply(models, function(m) {
      do.call(extrapolate, c(list(object = object, model = m), dots))
    }),
    models
  )

  tbl <- data.frame(
    model = models,
    asymptote = vapply(fits, function(f) f$asymptote, numeric(1)),
    AIC = vapply(fits, function(f) {
      if (!is.null(f$fit)) stats::AIC(f$fit) else NA_real_
    }, numeric(1)),
    BIC = vapply(fits, function(f) {
      if (!is.null(f$fit)) stats::BIC(f$fit) else NA_real_
    }, numeric(1)),
    n_params = vapply(fits, function(f) {
      if (!is.null(f$fit)) length(stats::coef(f$fit)) else NA_integer_
    }, numeric(1)),
    converged = vapply(fits, function(f) !is.null(f$fit), logical(1)),
    stringsAsFactors = FALSE
  )

  aic_vals <- tbl$AIC
  min_aic <- min(aic_vals, na.rm = TRUE)
  tbl$delta_AIC <- aic_vals - min_aic
  raw_weights <- exp(-0.5 * tbl$delta_AIC)
  tbl$AIC_weight <- raw_weights / sum(raw_weights, na.rm = TRUE)

  tbl <- tbl[order(tbl$AIC, na.last = TRUE), ]
  rownames(tbl) <- NULL

  best <- if (any(tbl$converged)) tbl$model[1] else NA_character_

  summ <- summary(object)
  df <- data.frame(x = summ$sites, y = summ$mean)

  structure(
    list(
      table = tbl,
      fits = fits,
      best_model = best,
      data = df,
      spacc = object
    ),
    class = "spacc_model_compare"
  )
}


#' @export
print.spacc_model_compare <- function(x, ...) {
  cat("SAR Model Comparison\n")
  cat(strrep("-", 30), "\n")

  tbl <- x$table
  for (i in seq_len(nrow(tbl))) {
    star <- if (!is.na(x$best_model) && tbl$model[i] == x$best_model) " *" else ""
    if (tbl$converged[i]) {
      cat(sprintf("  %-20s AIC: %8.1f  dAIC: %6.1f  w: %.3f  S_max: %6.1f%s\n",
                  tbl$model[i], tbl$AIC[i], tbl$delta_AIC[i],
                  tbl$AIC_weight[i], tbl$asymptote[i], star))
    } else {
      cat(sprintf("  %-20s (failed to converge)\n", tbl$model[i]))
    }
  }
  cat(strrep("-", 30), "\n")
  if (!is.na(x$best_model)) {
    cat(sprintf("Best model: %s\n", x$best_model))
  }
  invisible(x)
}


#' @export
summary.spacc_model_compare <- function(object, ...) {
  best <- object$best_model
  if (!is.na(best) && !is.null(object$fits[[best]]$fit)) {
    summary(object$fits[[best]]$fit)
  } else {
    print(object)
  }
}


#' @export
coef.spacc_model_compare <- function(object, ...) {
  best <- object$best_model
  if (!is.na(best) && !is.null(object$fits[[best]]$fit)) {
    stats::coef(object$fits[[best]]$fit)
  } else {
    NA
  }
}


#' @export
predict.spacc_model_compare <- function(object, n = NULL, ...) {
  best <- object$best_model
  if (!is.na(best)) {
    predict(object$fits[[best]], n = n, ...)
  } else {
    NA
  }
}


#' @export
as.data.frame.spacc_model_compare <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$table
}


#' @export
plot.spacc_model_compare <- function(x, extrapolate_to = NULL, ...) {
  check_suggests("ggplot2")

  df_obs <- x$data

  if (is.null(extrapolate_to)) {
    extrapolate_to <- max(df_obs$x) * 1.5
  }

  pred_x <- seq(1, extrapolate_to, length.out = 200)

  pred_list <- lapply(x$fits, function(fit) {
    if (is.null(fit$fit)) return(NULL)
    data.frame(
      x = pred_x,
      y = suppressWarnings(predict(fit, n = pred_x)),
      model = fit$model,
      stringsAsFactors = FALSE
    )
  })
  pred_df <- do.call(rbind, Filter(Negate(is.null), pred_list))

  if (is.null(pred_df) || nrow(pred_df) == 0) {
    stop("No models converged")
  }

  pred_df$best <- pred_df$model == x$best_model

  ggplot2::ggplot() +
    ggplot2::geom_point(
      data = df_obs,
      ggplot2::aes(x = x, y = y),
      color = "#2E7D32", alpha = 0.5
    ) +
    ggplot2::geom_line(
      data = df_obs,
      ggplot2::aes(x = x, y = y),
      color = "#2E7D32"
    ) +
    ggplot2::geom_line(
      data = pred_df[!pred_df$best, ],
      ggplot2::aes(x = x, y = y, color = model),
      linetype = "dashed", alpha = 0.5, linewidth = 0.5
    ) +
    ggplot2::geom_line(
      data = pred_df[pred_df$best, ],
      ggplot2::aes(x = x, y = y, color = model),
      linewidth = 1.2
    ) +
    ggplot2::labs(
      x = "Sites sampled",
      y = "Cumulative species",
      title = "SAR Model Comparison",
      subtitle = sprintf("Best model: %s", x$best_model),
      color = "Model"
    ) +
    spacc_theme()
}
