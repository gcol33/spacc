#' Beta Distance-Decay
#'
#' Compute the distance-decay of community similarity by fitting models to
#' pairwise beta dissimilarity as a function of geographic distance.
#'
#' @param x A site-by-species matrix (presence/absence or abundance), or a
#'   `spacc` object (will use its stored coordinates).
#' @param coords A data.frame with columns `x` and `y`, or a `spacc_dist`
#'   object.
#' @param index Character. Dissimilarity index: `"sorensen"` (default) or
#'   `"jaccard"`.
#' @param model Character. Decay model to fit: `"exponential"`, `"power"`,
#'   `"linear"`, or `"all"` (default). When `"all"`, fits all three and selects
#'   the best by AIC.
#' @param distance Character. Distance method: `"euclidean"` (default) or
#'   `"haversine"`.
#' @param n_bins Integer. Number of distance bins for binned means in plots.
#'   Default 50.
#' @param progress Logical. Show progress messages? Default `TRUE`.
#'
#' @return An object of class `spacc_beta_decay` containing:
#'   \item{pairs}{Data.frame with columns `distance`, `dissimilarity`, `site_i`,
#'     `site_j`}
#'   \item{fits}{Named list of fitted model objects}
#'   \item{best_model}{Name of best model by AIC}
#'   \item{half_life}{Distance at which similarity halves (exponential model only)}
#'   \item{coefficients}{Data.frame of model coefficients and AIC}
#'   \item{index}{Dissimilarity index used}
#'   \item{n_sites}{Number of sites}
#'   \item{n_pairs}{Number of pairwise comparisons}
#'
#' @details
#' Three decay models are available:
#'
#' - **Exponential**: \eqn{\beta = 1 - a \cdot e^{-b \cdot d}}
#' - **Power**: \eqn{\beta = a \cdot d^b}
#' - **Linear**: \eqn{\beta = a + b \cdot d}
#'
#' The half-life (exponential model) is the distance at which similarity
#' decays to half its initial value: \eqn{d_{1/2} = \ln(2) / b}.
#'
#' @references
#' Nekola, J.C. & White, P.S. (1999). The distance decay of similarity in
#' biogeography and ecology. Journal of Biogeography, 26, 867-878.
#'
#' Morlon, H., Chuyong, G., Condit, R., et al. (2008). A general framework
#' for the distance-decay of similarity in ecological communities. Ecology
#' Letters, 11, 904-917.
#'
#' @seealso [spaccBeta()] for spatial beta diversity accumulation,
#'   [distanceDecay()] for species similarity decay
#'
#' @examples
#' \donttest{
#' coords <- data.frame(x = runif(30), y = runif(30))
#' species <- matrix(rbinom(30 * 20, 1, 0.3), nrow = 30)
#'
#' bd <- betaDecay(species, coords)
#' print(bd)
#' plot(bd)
#' }
#'
#' @export
betaDecay <- function(x,
                      coords,
                      index = c("sorensen", "jaccard"),
                      model = c("all", "exponential", "power", "linear"),
                      distance = c("euclidean", "haversine"),
                      n_bins = 50L,
                      progress = TRUE) {

  index <- match.arg(index)
  model <- match.arg(model)
  distance <- match.arg(distance)

  x <- as.matrix(x)
  pa <- (x > 0) * 1L

  # Handle coords
  if (inherits(coords, "spacc_dist")) {
    dist_mat <- as.matrix(coords)
    coord_data <- attr(coords, "coords")
  } else {
    stopifnot("coords must have x and y columns" = all(c("x", "y") %in% names(coords)))
    coord_data <- coords
    if (progress) cli_info("Computing distance matrix")
    dist_mat <- cpp_distance_matrix(coord_data$x, coord_data$y, distance)
  }

  stopifnot("x and coords must have same rows" = nrow(x) == nrow(coord_data))

  n_sites <- nrow(pa)
  n_pairs <- n_sites * (n_sites - 1) / 2

  if (n_sites > 1000 && progress) {
    cli_info(sprintf("Large dataset: %d sites -> %.0f pairs. This may take a moment.",
                      n_sites, n_pairs))
  }

  # Compute all pairwise beta dissimilarities
  if (progress) cli_info(sprintf("Computing %s pairwise dissimilarities (%s)",
                                  format(n_pairs, big.mark = ","), index))

  pair_dist <- numeric(n_pairs)
  pair_beta <- numeric(n_pairs)
  pair_i <- integer(n_pairs)
  pair_j <- integer(n_pairs)
  idx <- 0L

  for (i in seq_len(n_sites - 1)) {
    sp_i <- pa[i, ]
    for (j in (i + 1):n_sites) {
      idx <- idx + 1L
      sp_j <- pa[j, ]

      a <- sum(sp_i & sp_j)  # shared
      b <- sum(sp_i & !sp_j) # only in i
      c_count <- sum(!sp_i & sp_j) # only in j

      if (index == "sorensen") {
        denom <- 2 * a + b + c_count
        pair_beta[idx] <- if (denom > 0) (b + c_count) / denom else 0
      } else {
        denom <- a + b + c_count
        pair_beta[idx] <- if (denom > 0) (b + c_count) / denom else 0
      }

      pair_dist[idx] <- dist_mat[i, j]
      pair_i[idx] <- i
      pair_j[idx] <- j
    }
  }

  pairs <- data.frame(
    distance = pair_dist,
    dissimilarity = pair_beta,
    site_i = pair_i,
    site_j = pair_j
  )

  # Fit models
  models_to_fit <- if (model == "all") c("exponential", "power", "linear") else model
  fits <- list()
  coef_df <- data.frame(
    model = character(0),
    a = numeric(0),
    b = numeric(0),
    AIC = numeric(0),
    stringsAsFactors = FALSE
  )

  if (progress) cli_info("Fitting decay models")

  for (m in models_to_fit) {
    fit <- tryCatch(fit_decay_model(pairs, m), error = function(e) NULL)
    if (!is.null(fit)) {
      fits[[m]] <- fit$fit
      coef_df <- rbind(coef_df, data.frame(
        model = m, a = fit$a, b = fit$b, AIC = fit$aic,
        stringsAsFactors = FALSE
      ))
    }
  }

  # Select best model
  best_model <- if (nrow(coef_df) > 0) coef_df$model[which.min(coef_df$AIC)] else NA_character_

  # Half-life for exponential
  half_life <- NA_real_
  if ("exponential" %in% names(fits)) {
    b_exp <- coef_df$b[coef_df$model == "exponential"]
    if (!is.na(b_exp) && b_exp > 0) {
      half_life <- log(2) / b_exp
    }
  }

  if (progress) cli_success("Done")

  structure(
    list(
      pairs = pairs,
      fits = fits,
      best_model = best_model,
      half_life = half_life,
      coefficients = coef_df,
      index = index,
      distance_type = distance,
      n_sites = n_sites,
      n_pairs = n_pairs,
      n_bins = n_bins,
      call = match.call()
    ),
    class = "spacc_beta_decay"
  )
}


# Fit a single decay model to pairwise data
fit_decay_model <- function(pairs, model) {
  d <- pairs$distance
  beta <- pairs$dissimilarity

  # Remove zero-distance pairs
  valid <- d > 0
  d <- d[valid]
  beta <- beta[valid]

  if (model == "exponential") {
    # beta = 1 - a * exp(-b * d)
    fit <- stats::nls(
      beta ~ 1 - a * exp(-b * d),
      start = list(a = 1, b = 1 / stats::median(d)),
      control = stats::nls.control(maxiter = 100, warnOnly = TRUE),
      algorithm = "port",
      lower = c(0, 0),
      upper = c(1, Inf)
    )
    cf <- stats::coef(fit)
    return(list(fit = fit, a = cf["a"], b = cf["b"],
                aic = stats::AIC(fit)))

  } else if (model == "power") {
    # log(beta) = log(a) + b * log(d)
    # Only use pairs with beta > 0
    pos <- beta > 0
    if (sum(pos) < 3) return(NULL)
    fit <- stats::lm(log(beta[pos]) ~ log(d[pos]))
    cf <- stats::coef(fit)
    return(list(fit = fit, a = exp(cf[1]), b = cf[2],
                aic = stats::AIC(fit)))

  } else if (model == "linear") {
    fit <- stats::lm(beta ~ d)
    cf <- stats::coef(fit)
    return(list(fit = fit, a = cf[1], b = cf[2],
                aic = stats::AIC(fit)))
  }
}


#' @export
print.spacc_beta_decay <- function(x, ...) {
  cat(sprintf("Beta distance-decay: %d sites, %s pairs\n",
              x$n_sites, format(x$n_pairs, big.mark = ",")))
  cat(sprintf("Index: %s, Distance: %s\n", x$index, x$distance_type))

  if (nrow(x$coefficients) > 0) {
    cat("\nFitted models:\n")
    for (i in seq_len(nrow(x$coefficients))) {
      row <- x$coefficients[i, ]
      best <- if (!is.na(x$best_model) && row$model == x$best_model) " *" else ""
      cat(sprintf("  %s: a=%.4f, b=%.6f, AIC=%.1f%s\n",
                  row$model, row$a, row$b, row$AIC, best))
    }
    if (!is.na(x$best_model)) cat(sprintf("\nBest model: %s\n", x$best_model))
    if (!is.na(x$half_life)) cat(sprintf("Half-life (exp): %.2f\n", x$half_life))
  } else {
    cat("No models fitted successfully.\n")
  }
  invisible(x)
}


#' @export
summary.spacc_beta_decay <- function(object, ...) {
  dist_range <- range(object$pairs$distance)
  beta_range <- range(object$pairs$dissimilarity)

  list(
    n_sites = object$n_sites,
    n_pairs = object$n_pairs,
    distance_range = dist_range,
    dissimilarity_range = beta_range,
    mean_dissimilarity = mean(object$pairs$dissimilarity),
    coefficients = object$coefficients,
    best_model = object$best_model,
    half_life = object$half_life,
    index = object$index
  )
}


#' @export
coef.spacc_beta_decay <- function(object, ...) {
  if (nrow(object$coefficients) == 0) return(NULL)
  row <- object$coefficients[object$coefficients$model == object$best_model, ]
  c(a = row$a, b = row$b)
}


#' @export
as.data.frame.spacc_beta_decay <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$pairs
}


#' @export
plot.spacc_beta_decay <- function(x, bins = TRUE, model = NULL, ...) {
  check_suggests("ggplot2")

  df <- x$pairs

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[["distance"]],
                                          y = .data[["dissimilarity"]])) +
    ggplot2::geom_point(alpha = 0.1, size = 0.5, color = "grey50")

  # Binned means
  if (bins) {
    breaks <- seq(min(df$distance), max(df$distance), length.out = x$n_bins + 1)
    df$bin <- cut(df$distance, breaks = breaks, include.lowest = TRUE)
    bin_means <- stats::aggregate(
      cbind(dissimilarity, distance) ~ bin,
      data = df[!is.na(df$bin), ],
      FUN = mean
    )

    p <- p + ggplot2::geom_point(
      data = bin_means,
      ggplot2::aes(x = .data[["distance"]], y = .data[["dissimilarity"]]),
      color = "#FF9800", size = 2
    )
  }

  # Fitted curve
  show_model <- if (!is.null(model)) model else x$best_model
  if (!is.na(show_model) && show_model %in% names(x$fits)) {
    d_seq <- seq(min(df$distance[df$distance > 0]), max(df$distance), length.out = 200)

    if (show_model == "exponential") {
      cf <- stats::coef(x$fits[["exponential"]])
      pred <- 1 - cf["a"] * exp(-cf["b"] * d_seq)
    } else if (show_model == "power") {
      cf <- x$coefficients[x$coefficients$model == "power", ]
      pred <- cf$a * d_seq^cf$b
    } else if (show_model == "linear") {
      cf <- stats::coef(x$fits[["linear"]])
      pred <- cf[1] + cf[2] * d_seq
    }

    pred_df <- data.frame(distance = d_seq, dissimilarity = pred)
    p <- p + ggplot2::geom_line(
      data = pred_df,
      ggplot2::aes(x = .data[["distance"]], y = .data[["dissimilarity"]]),
      color = "#C62828", linewidth = 1
    )
  }

  p +
    ggplot2::labs(
      x = "Geographic distance",
      y = paste0("Dissimilarity (", x$index, ")"),
      title = "Beta Distance-Decay",
      subtitle = if (!is.na(show_model)) paste("Model:", show_model) else NULL
    ) +
    spacc_theme()
}
