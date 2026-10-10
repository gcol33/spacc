#' Accumulation Rate Curve (ARC) Richness Estimator
#'
#' Estimate total species richness from the rate at which new species appear
#' along a random reordering of the observations (Shestopaloff 2025). The rate
#' curve is extrapolated with three decay curves whose tails give the number
#' of species not yet seen, combined with bounded weights.
#'
#' @param x A site-by-species matrix. For `data = "abundance"` the columns are
#'   pooled into species abundances; for `data = "incidence"` the matrix is
#'   binarized and each species counts the sites it occupies.
#' @param data Character. `"abundance"` (default) treats every individual as
#'   one observation; `"incidence"` treats every species-site occurrence as one.
#' @param K Integer. Species seen more than `K` times are left out of the rate
#'   curve. Default 20.
#' @param adjust Logical. Apply the smoothing adjustment of Eq. 11 with a
#'   log-linear fit to the frequency counts? Default `FALSE`. See Details.
#' @param n_boot Integer. Parametric bootstrap draws for the interval. Default
#'   2000. `0` returns the point estimate only.
#' @param parallel Logical. Use parallel processing? Default `TRUE`.
#' @param n_cores Integer. Number of cores.
#'
#' @return An object of class `spacc_estimate` with `estimator = "ARC"`. Its
#'   `details` hold `undetected` (undetected species from each curve) and
#'   `weights` (their weights), both averaged over the curve portions, `K`,
#'   `data`, `penalty` (the selected Poisson-mixture penalty) and `boot` (the
#'   bootstrap estimates).
#'
#' @details
#' Let \eqn{y_k} be the number of species observed \eqn{k} times and \eqn{N}
#' the number of observations of the species with \eqn{k \le K}. Drawing the
#' observations in random order, the expected number of these species seen
#' after \eqn{n} draws is
#' \deqn{S(n) = \sum_{k=1}^{K} y_k \left[1 - (1 - n/N)^k\right],}
#' and the accumulation rate is \eqn{s_n = S(n) - S(n-1)}, evaluated at up to
#' 1000 points. Three curves are fitted to the tail of \eqn{s_n} by least
#' squares with weights proportional to the squared point index: the
#' log-normal \eqn{1 - \Phi((\log n - \mu)/\sigma)}, the log-logistic
#' \eqn{1/(1 + \alpha n^\gamma)} and the Weibull \eqn{\exp(-\lambda n^\gamma)},
#' the last on the log scale. Each curve's integral from \eqn{N} to infinity
#' is its number of undetected species; the log-logistic integral is finite
#' only for \eqn{\gamma > 1}.
#'
#' The estimate combines the three with weights fitted by the squared misfit
#' to \eqn{s_n} plus a penalty on the product of the log-logistic and Weibull
#' weights (Eq. 8), each bounded by caps built from the three estimates and the
#' observed richness (Eq. 7). The undetected species are averaged over fits to
#' the last 100, 95, ..., 70 percent of the curve.
#'
#' The penalty is there to make the fit prefer either the log-logistic or the
#' Weibull, and gives the loss one local minimum for each. The weights are the
#' end point of the projected gradient flow of the loss started at zero
#' weights, which `arc()` solves in closed form face by face of the box of
#' admissible weights. This is the minimum the author's reference code reaches
#' (L-BFGS within an augmented Lagrangian, started at zero). In Poisson-Gamma
#' samples of intermediate sparsity it is the log-logistic one, with median
#' bias between -1 and +4 percent; the global minimum of the loss is the
#' Weibull one there, with median bias between -9 and -17 percent
#' (`vignette("richness-estimation")`).
#'
#' With `adjust = TRUE` the counts \eqn{y_1, \dots, y_{k_1}} (the first run of
#' positive counts) are mixed half and half with a log-linear fit
#' \eqn{\log y_k = b_0 + b_1 k}, and fully replaced by it; the adjusted
#' estimate is \eqn{2 S_{0.5} - S_1}. The adjustment assumes the ratios
#' \eqn{(k+1) y_{k+1} / y_k} are linear in \eqn{k}.
#'
#' The interval comes from a parametric bootstrap of the counts. Their
#' covariance is that of a mixture of 30 Poisson distributions whose weights
#' are fitted to \eqn{(y_1, \dots, y_K, y_{K+})} by penalised least squares
#' (Eq. 12-15). The penalty is chosen by resampling species: 1500 resamples
#' are ordered by their distance to the left-out species and split into
#' deciles, and the first decile whose out-of-sample error is smallest at a
#' non-zero penalty sets it. `n_boot` count vectors are drawn from a normal
#' distribution with that covariance, rounded to non-negative integers and
#' re-estimated; the interval is their 2.5 and 97.5 percent quantiles and the
#' standard error their standard deviation. In Poisson-Gamma simulations the
#' 95 percent interval covered the truth in 93 to 96 percent of samples for
#' sparse data and about 89 percent for dense data, where the point estimate
#' is biased upwards (`vignette("richness-estimation")`).
#'
#' With `data = "incidence"` the observations are species-site occurrences
#' and \eqn{y_k} becomes \eqn{Q_k}, the number of species found at \eqn{k}
#' sites, the substitution that turns Chao1 into Chao2.
#'
#' @references
#' Shestopaloff, K. (2025). An accumulation rate curve estimator for total
#' species. Environmental and Ecological Statistics, 32, 311-345.
#' \doi{10.1007/s10651-025-00644-y}
#'
#' @seealso [chao1()], [chao2()], [ace()]
#'
#' @examples
#' set.seed(1)
#' species <- matrix(rpois(40 * 300, rgamma(300, 0.4, 2)), nrow = 40, byrow = TRUE)
#' arc(species, n_boot = 0)
#' \donttest{
#' arc(species, n_boot = 200)
#' }
#'
#' @export
arc <- function(x,
                data = c("abundance", "incidence"),
                K = 20L,
                adjust = FALSE,
                n_boot = 2000L,
                parallel = TRUE,
                n_cores = NULL) {
  data <- match.arg(data)
  stopifnot("K must be a single integer of at least 2" =
              is.numeric(K) && length(K) == 1 && K >= 2 && K == round(K))
  K <- as.integer(K)
  n_cores <- resolve_cores(n_cores, parallel)

  x <- as.matrix(x)
  counts <- if (data == "abundance") colSums(x) else colSums(x > 0)
  counts <- counts[counts > 0]
  S_obs <- length(counts)
  y <- tabulate(pmin(counts, K + 1L), K + 1L)

  point <- cpp_arc_point(y, S_obs, adjust, 0.5, .arc_max_points)
  if (!is.finite(point$estimate)) {
    stop("ARC could not fit the rate curve: too few species seen at most K times",
         call. = FALSE)
  }

  interval <- if (n_boot > 0) .arc_interval(y, S_obs, adjust, n_boot, n_cores) else NULL

  structure(
    list(
      estimator = "ARC",
      estimate = point$estimate,
      se = if (is.null(interval)) NA_real_ else interval$se,
      lower = if (is.null(interval)) NA_real_ else interval$lower,
      upper = if (is.null(interval)) NA_real_ else interval$upper,
      S_obs = S_obs,
      details = list(undetected = point$undetected, weights = point$weights,
                     K = K, data = data, adjust = adjust,
                     penalty = interval$penalty, boot = interval$boot)
    ),
    class = "spacc_estimate"
  )
}

# Points of the rate curve used in the fits.
.arc_max_points <- 1000L


# Parametric bootstrap interval (Shestopaloff 2025, Sect. 3.4) for the counts
# y = (y_1, ..., y_K, y_K+).
.arc_interval <- function(y, S_obs, adjust, n_boot, n_cores,
                          n_lambda = 30L, n_penalty = 40L, n_resample = 1500L,
                          n_groups = 10L) {
  R <- length(y)
  K <- R - 1L

  # Poisson mixture components on a log grid; rows are k = 1..K and K+
  top <- min(K, max(which(y > 0)) + 2)
  lambda <- exp(seq(-6, log(top), length.out = n_lambda))
  A <- rbind(vapply(lambda, function(l) stats::dpois(seq_len(K), l), numeric(K)),
             stats::ppois(K, lambda, lower.tail = FALSE))
  penalties <- c(0, exp(seq(-2.5, 4.5, length.out = n_penalty))) / S_obs

  # Resample species; the left-out species, rescaled to S_obs, are the test set
  category <- rep(seq_len(R), y)
  in_sample <- out_sample <- matrix(0, n_resample, R)
  for (b in seq_len(n_resample)) {
    drawn <- sample.int(S_obs, S_obs, replace = TRUE)
    left_out <- setdiff(seq_len(S_obs), drawn)
    in_sample[b, ] <- tabulate(category[drawn], R)
    if (length(left_out) > 0) {
      out_sample[b, ] <- tabulate(category[left_out], R) * S_obs / length(left_out)
    }
  }
  ord <- order(rowMeans((in_sample - out_sample)^2))
  in_sample <- in_sample[ord, , drop = FALSE]
  out_sample <- out_sample[ord, , drop = FALSE]

  fits <- cpp_arc_penalty(A, in_sample, out_sample, penalties, y, n_cores)

  # first group whose mean out-of-sample error is smallest at a non-zero
  # penalty; the smallest non-zero penalty otherwise
  group <- rep(seq_len(n_groups), each = ceiling(n_resample / n_groups))[seq_len(n_resample)]
  group_err <- rowsum(fits$err, group) / as.vector(table(group))
  best <- apply(group_err, 1, which.min)
  pick <- if (any(best > 1L)) best[best > 1L][1] else 2L

  v <- fits$v_obs[, pick]
  Sigma <- diag(drop(A %*% v), R) - A %*% (v * t(A))

  eig <- eigen(Sigma, symmetric = TRUE)
  root <- eig$vectors %*% diag(sqrt(pmax(eig$values, 0)), R)
  draws <- matrix(stats::rnorm(n_boot * R), n_boot, R) %*% t(root)
  draws <- round(pmax(sweep(draws, 2, y, "+"), 0))

  boot <- cpp_arc_draws(draws, adjust, 0.5, .arc_max_points, n_cores)
  ok <- boot[is.finite(boot)]
  q <- stats::quantile(ok, c(0.025, 0.975), names = FALSE)
  list(se = stats::sd(ok), lower = q[1], upper = q[2],
       penalty = penalties[pick], boot = boot)
}
