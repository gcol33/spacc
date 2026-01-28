#' Extrapolate Total Species Richness
#'
#' Fit an asymptotic model to estimate total species richness beyond
#' the observed sampling effort.
#'
#' @param object A `spacc` object.
#' @param model Character. Model to fit: `"michaelis-menten"` (default),
#'   `"lomolino"`, `"asymptotic"`, `"weibull"`, `"logistic"`, or `"evt"`
#'   (Extreme Value Theory, Borda-de-Agua et al. 2025).
#' @param ... Additional arguments passed to [stats::nls()].
#'
#' @return An object of class `spacc_fit` containing:
#'   \item{asymptote}{Estimated total species richness}
#'   \item{asymptote_ci}{Confidence interval for asymptote}
#'   \item{model}{Model name}
#'   \item{fit}{The nls fit object}
#'   \item{aic}{AIC of the model}
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
#' @examples
#' \dontrun{
#' sac <- spacc(species, coords)
#' fit <- extrapolate(sac, model = "lomolino")
#'
#' print(fit)
#' plot(fit)
#' predict(fit, n = c(100, 500, 1000))
#' }
#'
#' @export
extrapolate <- function(object, model = c("michaelis-menten", "lomolino", "asymptotic", "weibull", "logistic", "evt"), ...) {

  model <- match.arg(model)

  stopifnot("object must be a spacc object" = inherits(object, "spacc"))

  summ <- summary(object)
  df <- data.frame(
    x = summ$sites,
    y = summ$mean
  )

  # Starting values and formulas for each model
  y_max <- max(df$y)
  x_half <- df$x[which.min(abs(df$y - y_max / 2))]

  fit <- tryCatch({
    switch(model,
      "michaelis-menten" = stats::nls(
        y ~ a * x / (b + x),
        data = df,
        start = list(a = y_max * 1.2, b = x_half),
        ...
      ),
      "lomolino" = stats::nls(
        y ~ a / (1 + b^(log(c/x))),
        data = df,
        start = list(a = y_max * 1.2, b = 2, c = x_half),
        control = stats::nls.control(maxiter = 200),
        ...
      ),
      "asymptotic" = stats::nls(
        y ~ a * (1 - exp(-b * x)),
        data = df,
        start = list(a = y_max * 1.2, b = 0.01),
        ...
      ),
      "weibull" = stats::nls(
        y ~ a * (1 - exp(-(x/b)^c)),
        data = df,
        start = list(a = y_max * 1.2, b = x_half, c = 1),
        control = stats::nls.control(maxiter = 200),
        ...
      ),
      "logistic" = stats::nls(
        y ~ a / (1 + exp(-b * (x - c))),
        data = df,
        start = list(a = y_max * 1.2, b = 0.1, c = x_half),
        ...
      ),
      "evt" = {
        # EVT-inspired three-phase SAR model (Borda-de-Agua et al. 2025)
        # Uses a mixture of two Weibull components to capture triphasic pattern
        # S(x) = a * (w * (1 - exp(-(x/b1)^c1)) + (1-w) * (1 - exp(-(x/b2)^c2)))
        if (nrow(df) < 10) {
          warning("EVT model requires >= 10 data points; falling back to Weibull")
          stats::nls(
            y ~ a * (1 - exp(-(x/b)^c)),
            data = df,
            start = list(a = y_max * 1.2, b = x_half, c = 1),
            control = stats::nls.control(maxiter = 200),
            ...
          )
        } else {
          stats::nls(
            y ~ a * (w * (1 - exp(-(x/b1)^c1)) + (1 - w) * (1 - exp(-(x/b2)^c2))),
            data = df,
            start = list(a = y_max * 1.3, w = 0.7,
                          b1 = max(1, x_half * 0.3), c1 = 1.5,
                          b2 = max(1, x_half * 2), c2 = 0.5),
            lower = c(y_max * 0.5, 0.01, 0.1, 0.1, 0.1, 0.1),
            upper = c(y_max * 5, 0.99, x_half * 10, 10, x_half * 20, 10),
            algorithm = "port",
            control = stats::nls.control(maxiter = 500, warnOnly = TRUE),
            ...
          )
        }
      }
    )
  }, error = function(e) {
    warning("Model fitting failed: ", e$message)
    NULL
  })

  if (is.null(fit)) {
    return(structure(
      list(
        asymptote = NA,
        asymptote_ci = c(NA, NA),
        model = model,
        fit = NULL,
        aic = NA,
        data = df,
        spacc = object
      ),
      class = "spacc_fit"
    ))
  }

  # Extract asymptote (parameter 'a' in all models)
  coefs <- stats::coef(fit)
  asymptote <- coefs["a"]

  # Confidence interval for asymptote
  ci <- tryCatch({
    stats::confint(fit)["a", ]
  }, error = function(e) c(NA, NA))

  structure(
    list(
      asymptote = asymptote,
      asymptote_ci = ci,
      model = model,
      fit = fit,
      aic = stats::AIC(fit),
      data = df,
      spacc = object
    ),
    class = "spacc_fit"
  )
}


#' @export
print.spacc_fit <- function(x, ...) {
  cat("Extrapolation:", x$model, "\n")
  cat(strrep("-", 30), "\n")
  if (is.na(x$asymptote)) {
    cat("Model fitting failed\n")
  } else {
    cat(sprintf("Estimated asymptote: %.1f species\n", x$asymptote))
    cat(sprintf("95%% CI: %.1f - %.1f\n", x$asymptote_ci[1], x$asymptote_ci[2]))
    cat(sprintf("AIC: %.1f\n", x$aic))
    cat(sprintf("Observed: %.1f species (%.0f%% of estimated)\n",
                max(x$data$y),
                100 * max(x$data$y) / x$asymptote))
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


#' @export
confint.spacc_fit <- function(object, parm, level = 0.95, ...) {
  if (!is.null(object$fit)) {
    stats::confint(object$fit, parm = parm, level = level, ...)
  } else {
    NA
  }
}


#' @export
predict.spacc_fit <- function(object, n = NULL, ...) {
  if (is.null(object$fit)) {
    return(NA)
  }

  if (is.null(n)) {
    n <- object$data$x
  }

  stats::predict(object$fit, newdata = data.frame(x = n))
}


#' @export
plot.spacc_fit <- function(x, extrapolate_to = NULL, ...) {
  check_suggests("ggplot2")

  df <- x$data

  # Default extrapolation
  if (is.null(extrapolate_to)) {
    extrapolate_to <- nrow(df) * 2
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(color = "#2E7D32", alpha = 0.5) +
    ggplot2::geom_line(color = "#2E7D32")

  # Add fitted curve
  if (!is.null(x$fit)) {
    pred_x <- seq(1, extrapolate_to, length.out = 200)
    pred_y <- predict(x, n = pred_x)
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
    ggplot2::theme_minimal(base_size = 12)
}
