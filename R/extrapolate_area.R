# ---- convex hull area (self-contained, no dependencies) ------------------
# Andrew's monotone-chain hull + shoelace area. Returns the area of the convex
# hull of the points; 0 for fewer than 3 non-collinear points.
.hull_area <- function(px, py) {
  pts <- unique(cbind(px, py))
  n <- nrow(pts)
  if (n < 3) return(0)
  px <- pts[, 1]; py <- pts[, 2]
  o <- order(px, py)
  px <- px[o]; py <- py[o]

  cross <- function(o1, o2, o3) {
    (px[o2] - px[o1]) * (py[o3] - py[o1]) - (py[o2] - py[o1]) * (px[o3] - px[o1])
  }
  build <- function(idx) {
    h <- integer(0)
    for (i in idx) {
      while (length(h) >= 2 &&
             cross(h[length(h) - 1], h[length(h)], i) <= 0) {
        h <- h[-length(h)]
      }
      h <- c(h, i)
    }
    h
  }
  lower <- build(seq_len(n))
  upper <- build(rev(seq_len(n)))
  hull <- c(lower[-length(lower)], upper[-length(upper)])
  if (length(hull) < 3) return(0)

  hx <- px[hull]; hy <- py[hull]
  abs(sum(hx * c(hy[-1], hy[1]) - c(hx[-1], hx[1]) * hy)) / 2
}

# Median-bias-corrected percentile interval. The site bootstrap resamples with
# replacement, which loses ~37% of unique sites and biases every replicate's
# richness downward; this keeps the bootstrap spread but recentres it on the
# full-sample point estimate so the interval brackets the point.
.centered_ci <- function(boot_vals, point, level) {
  boot_vals <- boot_vals[is.finite(boot_vals)]
  if (length(boot_vals) < 2) return(c(NA_real_, NA_real_))
  alpha <- (1 - level) / 2
  shift <- point - stats::median(boot_vals)
  unname(stats::quantile(boot_vals, c(alpha, 1 - alpha)) + shift)
}

# Combinations of `a` subareas out of H: enumerate when there are few, sample
# otherwise. Returns a list of integer index vectors.
.choose_combos <- function(H, a, n_combos) {
  if (choose(H, a) <= n_combos) {
    utils::combn(H, a, simplify = FALSE)
  } else {
    unique(lapply(seq_len(n_combos), function(.) sort(sample.int(H, a))))
  }
}

# Build the Ugland total-species (T-S) curve: partition sites into spatial
# subareas, then for each number of subareas average union richness and union
# (convex-hull) area over random combinations.
.ts_curve <- function(pa, cx, cy, n_subareas, n_combos) {
  cl <- stats::kmeans(scale(cbind(cx, cy)), centers = n_subareas,
                      nstart = 5, iter.max = 50)$cluster
  subs <- sort(unique(cl))
  H <- length(subs)

  area <- numeric(H); rich <- numeric(H); rsd <- numeric(H)
  for (a in seq_len(H)) {
    combos <- .choose_combos(H, a, n_combos)
    ar <- numeric(length(combos)); ri <- numeric(length(combos))
    for (j in seq_along(combos)) {
      mem <- which(cl %in% subs[combos[[j]]])
      ri[j] <- sum(colSums(pa[mem, , drop = FALSE]) > 0)
      ar[j] <- .hull_area(cx[mem], cy[mem])
    }
    area[a] <- mean(ar)
    rich[a] <- mean(ri)
    rsd[a]  <- if (length(ri) > 1) stats::sd(ri) else 0
  }
  data.frame(n_sub = seq_len(H), area = area, richness = rich, richness_sd = rsd)
}


#' Extrapolate Richness to a Larger Area (Ugland T-S Curve)
#'
#' Estimate species richness for a spatial extent larger than the one sampled,
#' using the total-species (T-S) curve of Ugland, Gray & Ellingsen (2003).
#' Sites are partitioned into spatial subareas; the T-S curve is the expected
#' total richness of random combinations of subareas plotted against their
#' (convex-hull) area; an asymptotic model is fitted to that curve and
#' extrapolated to a target area, with a site-bootstrap confidence interval.
#'
#' Unlike [extrapolate()], which extends the curve in *sample count*, this
#' answers "how many species in a region of a given *area*", including areas
#' beyond the sampled extent.
#'
#' @param x A site-by-species matrix (presence/absence or abundance; binarized).
#' @param coords A two-column matrix or data frame of site coordinates (x, y).
#' @param target_area Numeric. Area to extrapolate to, in the squared units of
#'   `coords`. Defaults to twice the observed (convex-hull) area.
#' @param n_subareas Integer. Number of spatial subareas to partition sites into
#'   (k-means on coordinates). Default 10; must be between 4 and the number of
#'   sites.
#' @param model Character. Area model to fit: `"michaelis-menten"` (default),
#'   `"asymptotic"`, `"weibull"`, `"lomolino"`, or `"logistic"`.
#' @param n_combos Integer. Number of random subarea combinations averaged per
#'   T-S point (all combinations are used when there are fewer). Default 30.
#' @param R Integer. Number of site-bootstrap replicates for the confidence
#'   interval. Default 200.
#' @param level Numeric. Confidence level. Default 0.95.
#' @param seed Optional integer for reproducibility.
#' @param progress Logical. Show a progress message. Default `TRUE`.
#'
#' @return An object of class `spacc_area` with components:
#'   \item{ts_curve}{Data frame of the T-S curve (`n_sub`, `area`, `richness`,
#'     `richness_sd`)}
#'   \item{fit}{The fitted `nls` area model}
#'   \item{observed_area, observed_richness}{Convex-hull area and richness of the
#'     full sample}
#'   \item{target_area}{Area extrapolated to}
#'   \item{predicted_richness, predicted_ci}{Extrapolated richness and its
#'     bootstrap interval}
#'   \item{asymptote, asymptote_ci}{Model asymptote and its bootstrap interval}
#'
#' @references
#' Ugland, K.I., Gray, J.S. & Ellingsen, K.E. (2003). The species-accumulation
#' curve and estimation of species richness. Journal of Animal Ecology, 72,
#' 888-897. \doi{10.1046/j.1365-2656.2003.00748.x}
#'
#' @seealso [extrapolate()] for sample-count extrapolation; [chao2()] for
#'   nonparametric richness.
#'
#' @examples
#' \donttest{
#' coords <- data.frame(x = runif(150), y = runif(150))
#' species <- matrix(rbinom(150 * 60, 1, 0.15), nrow = 150)
#' area_fit <- extrapolateArea(species, coords, n_subareas = 8, progress = FALSE)
#' print(area_fit)
#' }
#'
#' @export
extrapolateArea <- function(x, coords,
                            target_area = NULL,
                            n_subareas = 10L,
                            model = c("michaelis-menten", "asymptotic",
                                      "weibull", "lomolino", "logistic"),
                            n_combos = 30L,
                            R = 200L,
                            level = 0.95,
                            seed = NULL,
                            progress = TRUE) {

  model <- match.arg(model)
  if (!is.null(seed)) set.seed(seed)

  x <- as.matrix(x)
  pa <- (x > 0) * 1L
  coords <- as.matrix(coords[, 1:2])
  cx <- as.numeric(coords[, 1]); cy <- as.numeric(coords[, 2])
  n_sites <- nrow(pa)

  stopifnot(
    "coords must have one row per site" = length(cx) == n_sites,
    "n_subareas must be >= 4" = n_subareas >= 4,
    "n_subareas cannot exceed the number of sites" = n_subareas <= n_sites
  )
  n_subareas <- as.integer(n_subareas)

  A_obs <- .hull_area(cx, cy)
  if (A_obs <= 0) stop("Coordinates are degenerate (zero hull area).", call. = FALSE)
  S_obs <- sum(colSums(pa) > 0)
  if (is.null(target_area)) target_area <- 2 * A_obs

  if (isTRUE(progress)) {
    cli_info(sprintf("Building T-S curve (%d subareas) and bootstrapping (R=%d)",
                     n_subareas, R))
  }

  fit_ts <- function(ts) {
    df <- ts[ts$area > 0, c("area", "richness")]
    names(df) <- c("x", "y")
    .fit_sar(df, model, quiet = TRUE)
  }

  ts <- .ts_curve(pa, cx, cy, n_subareas, n_combos)
  fit <- fit_ts(ts)
  if (is.null(fit)) {
    stop("Area model failed to converge; try a different model or fewer subareas.",
         call. = FALSE)
  }
  coefs <- stats::coef(fit)
  asymptote <- unname(coefs["a"])
  pred <- unname(.sar_eval(model, coefs, target_area))

  # ---- site bootstrap ----
  param_n <- length(.sar_params(model))
  boot_coefs <- matrix(NA_real_, nrow = R, ncol = param_n,
                       dimnames = list(NULL, .sar_params(model)))
  for (b in seq_len(R)) {
    idx <- sample.int(n_sites, replace = TRUE)
    tb <- tryCatch(.ts_curve(pa[idx, , drop = FALSE], cx[idx], cy[idx],
                             n_subareas, n_combos),
                   error = function(e) NULL)
    if (is.null(tb)) next
    fb <- fit_ts(tb)
    if (!is.null(fb)) boot_coefs[b, ] <- stats::coef(fb)
  }
  boot_coefs <- boot_coefs[stats::complete.cases(boot_coefs), , drop = FALSE]

  if (nrow(boot_coefs) >= 2) {
    asym_ci <- .centered_ci(boot_coefs[, "a"], asymptote, level)
    boot_pred <- apply(boot_coefs, 1, function(cf) .sar_eval(model, cf, target_area))
    pred_ci <- .centered_ci(boot_pred, pred, level)
  } else {
    asym_ci <- c(NA_real_, NA_real_)
    pred_ci <- c(NA_real_, NA_real_)
  }

  structure(
    list(
      ts_curve = ts,
      fit = fit,
      model = model,
      observed_area = A_obs,
      observed_richness = S_obs,
      target_area = target_area,
      predicted_richness = pred,
      predicted_ci = pred_ci,
      asymptote = asymptote,
      asymptote_ci = asym_ci,
      level = level,
      n_subareas = nrow(ts),
      boot = list(coefs = boot_coefs, n_ok = nrow(boot_coefs), n_fail = R - nrow(boot_coefs))
    ),
    class = "spacc_area"
  )
}


#' @export
print.spacc_area <- function(x, ...) {
  cat("Area Extrapolation (Ugland T-S curve)\n")
  cat(strrep("-", 38), "\n")
  cat(sprintf("Model:               %s\n", x$model))
  cat(sprintf("Subareas:            %d\n", x$n_subareas))
  cat(sprintf("Observed:            %d species over area %.4g\n",
              x$observed_richness, x$observed_area))
  cat(sprintf("Asymptote:           %.1f species", x$asymptote))
  if (all(is.finite(x$asymptote_ci))) {
    cat(sprintf("  (%.0f%% CI %.1f - %.1f)", 100 * x$level,
                x$asymptote_ci[1], x$asymptote_ci[2]))
  }
  cat("\n")
  cat(sprintf("Predicted at area %.4g: %.1f species", x$target_area,
              x$predicted_richness))
  if (all(is.finite(x$predicted_ci))) {
    cat(sprintf("  (%.0f%% CI %.1f - %.1f)", 100 * x$level,
                x$predicted_ci[1], x$predicted_ci[2]))
  }
  cat("\n")
  invisible(x)
}


#' @export
summary.spacc_area <- function(object, ...) {
  data.frame(
    model = object$model,
    n_subareas = object$n_subareas,
    observed_area = object$observed_area,
    observed_richness = object$observed_richness,
    target_area = object$target_area,
    predicted_richness = object$predicted_richness,
    lower = object$predicted_ci[1],
    upper = object$predicted_ci[2]
  )
}


#' @export
as.data.frame.spacc_area <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$ts_curve
}


#' @export
coef.spacc_area <- function(object, ...) stats::coef(object$fit)


#' Predict richness at a given area
#'
#' @param object A `spacc_area` object.
#' @param area Numeric vector of areas. Defaults to the observed T-S areas.
#' @param interval `"none"` (default) or `"bootstrap"` for percentile bounds.
#' @param level Confidence level. Defaults to the fit's level.
#' @param ... Unused.
#' @return A numeric vector (point) or a data frame with `area`, `fit`, `lower`,
#'   `upper`.
#' @export
predict.spacc_area <- function(object, area = NULL, interval = c("none", "bootstrap"),
                               level = NULL, ...) {
  interval <- match.arg(interval)
  if (is.null(area)) area <- object$ts_curve$area
  point <- .sar_eval(object$model, stats::coef(object$fit), area)

  if (interval == "none" || is.null(object$boot) || object$boot$n_ok < 2) {
    return(point)
  }
  if (is.null(level)) level <- object$level
  pm <- vapply(seq_len(nrow(object$boot$coefs)),
               function(r) .sar_eval(object$model, object$boot$coefs[r, ], area),
               numeric(length(area)))
  pm <- matrix(pm, nrow = length(area))
  ci <- vapply(seq_along(area), function(i) .centered_ci(pm[i, ], point[i], level),
               numeric(2))
  data.frame(
    area = area,
    fit = point,
    lower = ci[1, ],
    upper = ci[2, ]
  )
}


#' @export
plot.spacc_area <- function(x, extrapolate_to = NULL, ...) {
  check_suggests("ggplot2")

  ts <- x$ts_curve
  if (is.null(extrapolate_to)) extrapolate_to <- x$target_area
  grid <- seq(min(ts$area[ts$area > 0]), extrapolate_to, length.out = 200)

  band <- suppressWarnings(predict(x, area = grid, interval = "bootstrap"))
  if (!is.data.frame(band)) band <- data.frame(area = grid, fit = band, lower = NA, upper = NA)

  p <- ggplot2::ggplot()
  if (all(!is.na(band$lower))) {
    p <- p + ggplot2::geom_ribbon(
      data = band, ggplot2::aes(x = area, ymin = lower, ymax = upper),
      fill = "#FF9800", alpha = 0.2
    )
  }
  p +
    ggplot2::geom_line(data = band, ggplot2::aes(x = area, y = fit),
                       color = "#FF9800", linewidth = 1, linetype = "dashed") +
    ggplot2::geom_point(data = ts, ggplot2::aes(x = area, y = richness),
                        color = "#2E7D32", size = 2) +
    ggplot2::geom_vline(xintercept = x$observed_area, linetype = "dotted",
                        color = "#78909C") +
    ggplot2::annotate("point", x = x$target_area, y = x$predicted_richness,
                      color = "#F44336", size = 3) +
    ggplot2::labs(
      x = "Area", y = "Total species",
      title = "Area-Based Richness Extrapolation",
      subtitle = sprintf("Ugland T-S curve, model: %s", x$model)
    ) +
    spacc_theme()
}
