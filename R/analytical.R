#' Analytical Accumulation Methods
#'
#' These methods compute expected species accumulation without simulation.
#' They are faster but don't provide spatial information.
#'
#' @name analytical
#' @rdname analytical
NULL


#' Coleman Expected Accumulation
#'
#' Compute the expected species accumulation curve using the Coleman method
#' (Coleman et al. 1982). This is an analytical formula, no simulation needed.
#'
#' @param x A site-by-species matrix (presence/absence or abundance).
#' @return A data.frame with columns: sites, expected, sd
#'
#' @references
#' Coleman, B.D., Mares, M.A., Willig, M.R. & Hsieh, Y.H. (1982).
#' Randomness, area, and species richness. Ecology, 63, 1121-1133.
#'
#' @export
coleman <- function(x) {
  x <- as.matrix(x)
  x <- (x > 0) * 1  # presence/absence

  n_sites <- nrow(x)
  n_species <- ncol(x)

  # Frequency of each species (number of sites where present)
  freq <- colSums(x)

  # Expected richness at each sample size
  expected <- numeric(n_sites)
  variance <- numeric(n_sites)

  for (n in 1:n_sites) {
    # E[S_n] = sum_i (1 - (1 - f_i/N)^n) where f_i = frequency of species i
    # But Coleman uses: E[S_n] = sum_i (1 - choose(N-f_i, n)/choose(N, n))
    probs <- 1 - exp(lchoose(n_sites - freq, n) - lchoose(n_sites, n))
    probs[freq == 0] <- 0

    expected[n] <- sum(probs)

    # Variance (Coleman equation)
    var_terms <- probs * (1 - probs)
    variance[n] <- sum(var_terms)
  }

  data.frame(
    sites = 1:n_sites,
    expected = expected,
    sd = sqrt(variance)
  )
}


#' Exact (Mao Tau) Expected Accumulation
#'
#' Compute the expected species accumulation curve using sample-based
#' rarefaction (Mao Tau estimator). This is analytically identical to
#' the expected curve from random permutations.
#'
#' @param x A site-by-species matrix (presence/absence or abundance).
#' @return A data.frame with columns: sites, expected, sd, lower, upper
#'
#' @references
#' Colwell, R.K., Mao, C.X. & Chang, J. (2004). Interpolating, extrapolating,
#' and comparing incidence-based species accumulation curves. Ecology, 85, 2717-2727.
#'
#' @export
mao_tau <- function(x) {
  x <- as.matrix(x)
  x <- (x > 0) * 1  # presence/absence

  n_sites <- nrow(x)
  n_species <- ncol(x)

  # Frequency of each species
  freq <- colSums(x)
  freq <- freq[freq > 0]  # remove absent species
  S_obs <- length(freq)

  expected <- numeric(n_sites)
  variance <- numeric(n_sites)

  for (n in 1:n_sites) {
    # Mao Tau: E[S_n] = S_obs - sum_i choose(N-f_i, n)/choose(N, n)
    # Equivalent to: sum_i (1 - choose(N-f_i, n)/choose(N, n))
    probs <- 1 - exp(lchoose(n_sites - freq, n) - lchoose(n_sites, n))
    expected[n] <- sum(probs)

    # Variance (Colwell et al. 2004, eq. 6)
    if (n < n_sites) {
      cov_sum <- 0
      for (i in seq_along(freq)) {
        for (j in seq_along(freq)) {
          if (i != j) {
            # Covariance term
            pi <- 1 - exp(lchoose(n_sites - freq[i], n) - lchoose(n_sites, n))
            pj <- 1 - exp(lchoose(n_sites - freq[j], n) - lchoose(n_sites, n))
            pij <- 1 - exp(lchoose(n_sites - freq[i], n) - lchoose(n_sites, n)) -
                       exp(lchoose(n_sites - freq[j], n) - lchoose(n_sites, n)) +
                       exp(lchoose(n_sites - freq[i] - freq[j], n) - lchoose(n_sites, n))
            pij <- max(0, min(1, pij))  # bound to [0,1]
            cov_sum <- cov_sum + (pij - pi * pj)
          }
        }
      }
      var_terms <- probs * (1 - probs)
      variance[n] <- sum(var_terms) + cov_sum
    } else {
      variance[n] <- 0
    }
  }

  variance[variance < 0] <- 0
  sd_vals <- sqrt(variance)

  data.frame(
    sites = 1:n_sites,
    expected = expected,
    sd = sd_vals,
    lower = expected - 1.96 * sd_vals,
    upper = expected + 1.96 * sd_vals
  )
}


#' Collector's Curve
#'
#' Compute the species accumulation curve in the order sites appear
#' in the data (no randomization). Useful for understanding how data
#' was collected.
#'
#' @param x A site-by-species matrix.
#' @return A data.frame with columns: sites, species
#'
#' @export
collector <- function(x) {
  x <- as.matrix(x)
  x <- (x > 0) * 1

  n_sites <- nrow(x)
  species_seen <- logical(ncol(x))
  curve <- integer(n_sites)

  for (i in 1:n_sites) {
    species_seen <- species_seen | (x[i, ] > 0)
    curve[i] <- sum(species_seen)
  }

  data.frame(
    sites = 1:n_sites,
    species = curve
  )
}


#' Spatially-Constrained Rarefaction
#'
#' Compute expected species richness accounting for spatial autocorrelation
#' (Chiarucci et al. 2009). Uses distance-weighted sampling probabilities.
#'
#' @param x A site-by-species matrix.
#' @param coords A data.frame with x and y columns.
#' @param n_perm Number of permutations. Default 100.
#' @param bandwidth Distance bandwidth for spatial weighting.
#'
#' @return A data.frame with columns: sites, mean, sd, lower, upper
#'
#' @references
#' Chiarucci, A., Bacaro, G., Rocchini, D. & Fattorini, L. (2009).
#' Discovering and rediscovering the sample-based rarefaction formula in
#' the ecological literature. Community Ecology, 10, 195-199.
#'
#' @export
spatial_rarefaction <- function(x, coords, n_perm = 100, bandwidth = NULL) {
  x <- as.matrix(x)
  x <- (x > 0) * 1L

  n_sites <- nrow(x)

  # Compute distance matrix
  dist_mat <- as.matrix(stats::dist(cbind(coords$x, coords$y)))

  # Default bandwidth: median distance
 if (is.null(bandwidth)) {
    bandwidth <- stats::median(dist_mat[lower.tri(dist_mat)])
  }

  # Run permutations with spatially-weighted sampling
  curves <- matrix(0, n_perm, n_sites)

  for (p in 1:n_perm) {
    # Start from random site
    current <- sample(n_sites, 1)
    visited <- logical(n_sites)
    visited[current] <- TRUE
    species_seen <- x[current, ] > 0

    curves[p, 1] <- sum(species_seen)

    for (step in 2:n_sites) {
      # Weight unvisited sites by distance from current
      weights <- exp(-dist_mat[current, ] / bandwidth)
      weights[visited] <- 0
      weights <- weights / sum(weights)

      # Sample next site
      current <- sample(n_sites, 1, prob = weights)
      visited[current] <- TRUE
      species_seen <- species_seen | (x[current, ] > 0)

      curves[p, step] <- sum(species_seen)
    }
  }

  data.frame(
    sites = 1:n_sites,
    mean = colMeans(curves),
    sd = apply(curves, 2, stats::sd),
    lower = apply(curves, 2, stats::quantile, 0.025),
    upper = apply(curves, 2, stats::quantile, 0.975)
  )
}
