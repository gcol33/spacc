# Parameter-recovery and CI-coverage tests (issue #2).
#
# These go beyond the structural tests (shape / class / no-NaN) elsewhere: they
# check that the statistical machinery recovers a known truth and that the
# defining-formula endpoints of the accumulation curves are correct.

# ---- simulators ----------------------------------------------------------
sim_clustered <- function(S = 240, P = 120, ell = 0.12) {
  xy <- data.frame(x = runif(S), y = runif(S))
  cx <- runif(P); cy <- runif(P)
  m <- matrix(0L, S, P)
  for (k in seq_len(P)) {
    m[, k] <- rbinom(S, 1, exp(-((xy$x - cx[k])^2 + (xy$y - cy[k])^2) / (2 * ell^2)))
  }
  list(xy = xy, comm = m[, colSums(m) > 0, drop = FALSE])
}

sim_easy <- function(S = 200, P = 70) {
  xy <- data.frame(x = runif(S), y = runif(S))
  p <- runif(P, 0.08, 0.45)
  m <- sapply(p, function(pi) rbinom(S, 1, pi))
  list(xy = xy, comm = m[, colSums(m) > 0, drop = FALSE])
}

# Incidence pool: S_true species with occupancy probabilities over n sites;
# truth = S_true (including species missed by the sample).
pool_incidence <- function(S_true = 150, n = 80, pmin = 0.02, pmax = 0.30) {
  p <- runif(S_true, pmin, pmax)
  sapply(p, function(pi) rbinom(n, 1, pi))
}

# Abundance pool: S_true species with lognormal relative abundances, Nind
# individuals sampled multinomially into a single pooled sample.
pool_abundance <- function(S_true = 150, Nind = 2000, sdlog = 1.0) {
  ra <- stats::rlnorm(S_true, 0, sdlog); ra <- ra / sum(ra)
  matrix(as.integer(stats::rmultinom(1, Nind, ra)), nrow = 1)
}

recover_stats <- function(simfun, est, N, truth_fun = function(m) ncol(m)) {
  cover <- logical(N); bias <- numeric(N)
  for (i in seq_len(N)) {
    m <- simfun(); truth <- truth_fun(m)
    obs <- m[, colSums(m) > 0, drop = FALSE]
    e <- est(obs)
    bias[i] <- (e$estimate - truth) / truth
    cover[i] <- truth >= e$lower && truth <= e$upper
  }
  list(bias = median(bias), coverage = mean(cover))
}


# ---- incidence estimators recover truth ----------------------------------
test_that("incidence richness estimators recover the true pool size", {
  skip_on_cran()
  set.seed(101)
  N <- 40
  c2 <- recover_stats(pool_incidence, chao2, N)
  i2 <- recover_stats(pool_incidence, iChao2, N)
  j1 <- recover_stats(pool_incidence, function(x) jackknife(x, 1), N)
  j2 <- recover_stats(pool_incidence, function(x) jackknife(x, 2), N)

  # all near-unbiased
  expect_lt(abs(c2$bias), 0.06)
  expect_lt(abs(i2$bias), 0.06)
  expect_lt(abs(j1$bias), 0.06)
  expect_lt(abs(j2$bias), 0.06)

  # jackknife intervals achieve ~nominal coverage
  expect_gt(j1$coverage, 0.85)
  expect_gt(j2$coverage, 0.85)
  # chao2 / iChao2 are conservative lower-bound estimators: their log-CIs
  # under-cover (the issue observed ~68%), but clearly beat the broken 0/40.
  expect_gt(c2$coverage, 0.45)
  expect_gt(i2$coverage, 0.45)
})


# ---- abundance estimators recover truth ----------------------------------
test_that("abundance richness estimators recover the true pool size", {
  skip_on_cran()
  set.seed(106)
  N <- 40
  c1 <- recover_stats(pool_abundance, chao1, N)
  i1 <- recover_stats(pool_abundance, iChao1, N)
  ac <- recover_stats(pool_abundance, ace, N)

  expect_lt(abs(c1$bias), 0.06)
  expect_lt(abs(i1$bias), 0.06)
  expect_lt(abs(ac$bias), 0.07)

  expect_gt(c1$coverage, 0.60)
  expect_gt(i1$coverage, 0.75)
  expect_gt(ac$coverage, 0.75)
})


# ---- extrapolate(): recovery on a near-saturation regime -----------------
test_that("extrapolate() asymptote recovers truth on well-sampled data", {
  skip_on_cran()
  set.seed(202)
  N <- 20
  bias <- numeric(N)
  for (i in seq_len(N)) {
    d <- sim_easy(); truth <- ncol(d$comm)
    keep <- sample(nrow(d$comm), floor(nrow(d$comm) * 0.8))
    comm <- d$comm[keep, , drop = FALSE]; comm <- comm[, colSums(comm) > 0, drop = FALSE]
    sac <- spacc(comm, d$xy[keep, ], n_seeds = 40, progress = FALSE)
    fit <- suppressWarnings(extrapolate(sac, model = "asymptotic",
                                        interval = "none", compare = FALSE,
                                        warn_ratio = Inf))
    bias[i] <- (fit$asymptote - truth) / truth
  }
  expect_lt(abs(median(bias)), 0.15)
})


# ---- extrapolate(): bootstrap interval is wider than the profile one ------
# (issue #1: the profile interval is over-confident; the bootstrap interval
#  propagates across-seed variability and must not be narrower.)
test_that("bootstrap interval is wider than the nls profile interval", {
  skip_on_cran()
  set.seed(303)
  d <- sim_clustered()
  keep <- sample(nrow(d$comm), nrow(d$comm) %/% 2)
  sac <- spacc(d$comm[keep, ], d$xy[keep, ], method = "kncn",
               n_seeds = 60, progress = FALSE)
  fp <- suppressWarnings(extrapolate(sac, model = "michaelis-menten",
                                     interval = "profile", compare = FALSE,
                                     warn_ratio = Inf))
  fb <- suppressWarnings(extrapolate(sac, model = "michaelis-menten",
                                     interval = "bootstrap", R = 150,
                                     compare = FALSE, warn_ratio = Inf))
  w_prof <- diff(fp$asymptote_ci)
  w_boot <- diff(fb$asymptote_ci)
  expect_true(is.finite(w_boot))
  expect_gt(w_boot, w_prof)
})


# ---- extrapolate(): diagnostics fire on implausible extrapolation --------
test_that("extrapolate() warns on aggressive extrapolation and disagreement", {
  set.seed(404)
  d <- sim_clustered()
  keep <- sample(nrow(d$comm), nrow(d$comm) %/% 2)
  sac <- spacc(d$comm[keep, ], d$xy[keep, ], method = "kncn",
               n_seeds = 40, progress = FALSE)
  # lomolino strongly overshoots clustered data -> ratio + disagreement warnings
  w <- testthat::capture_warnings(
    extrapolate(sac, model = "lomolino", interval = "none", warn_ratio = 1.5)
  )
  expect_true(any(grepl("exceeds observed richness", w)))
  expect_true(any(grepl("disagrees with nonparametric", w)))
  fit <- suppressWarnings(extrapolate(sac, model = "michaelis-menten",
                                      interval = "bootstrap", R = 60))
  expect_warning(predict(fit, n = 100 * sac$n_sites), "beyond")
})


# ---- extrapolateArea(): area extrapolation is sane -----------------------
test_that("extrapolateArea() interpolates the observed area and extrapolates up", {
  skip_on_cran()
  set.seed(505)
  d <- sim_clustered(S = 400, P = 120)
  af <- extrapolateArea(d$comm, d$xy, n_subareas = 8,
                        model = "michaelis-menten", R = 120,
                        seed = 1, progress = FALSE)
  # predicting at the observed area reproduces observed richness within tolerance
  at_obs <- predict(af, area = af$observed_area)
  expect_lt(abs(at_obs - af$observed_richness) / af$observed_richness, 0.15)
  # extrapolating to a larger area predicts at least as many species
  expect_gte(af$predicted_richness, af$observed_richness - 1e-6)
  # the centred interval brackets the point estimate
  expect_lte(af$predicted_ci[1], af$predicted_richness + 1e-6)
  expect_gte(af$predicted_ci[2], af$predicted_richness - 1e-6)
  # T-S curve is monotone in area
  ts <- af$ts_curve
  expect_true(all(diff(ts$richness) >= -1e-6))
})


# ---- spaccCoverage(): endpoint matches the Chao coverage formula ---------
test_that("spaccCoverage endpoint richness and coverage match definitions", {
  set.seed(606)
  n <- 14; S <- 80
  coords <- data.frame(x = runif(n), y = runif(n))
  ab <- matrix(rpois(n * S, 0.15), n, S)         # few sites -> pooled singletons
  ab <- ab[, colSums(ab) > 0, drop = FALSE]
  S_obs <- ncol(ab)

  ab_tot <- colSums(ab)
  testthat::skip_if(sum(ab_tot == 1) == 0 || sum(ab_tot == 2) == 0,
                    "no singletons/doubletons in this draw")

  cv <- spaccCoverage(ab, coords, n_seeds = 6, progress = FALSE)
  end_rich <- cv$richness[, cv$n_sites]
  end_C <- cv$coverage[, cv$n_sites]

  expect_equal(unname(mean(end_rich)), S_obs)

  N <- sum(ab_tot); f1 <- sum(ab_tot == 1); f2 <- sum(ab_tot == 2)
  C_def <- 1 - (f1 / N) * ((N - 1) * f1 / ((N - 1) * f1 + 2 * f2))
  expect_equal(unname(mean(end_C)), C_def, tolerance = 1e-4)

  # interpolating to the observed coverage returns the observed richness
  interp <- interpolateCoverage(cv, target = mean(end_C))
  expect_equal(unname(mean(interp[[1]])), S_obs, tolerance = 0.5)
})


# ---- Hill endpoints equal the pooled Hill numbers ------------------------
test_that("spaccHill endpoints equal abundance-based Hill numbers of the pool", {
  set.seed(707)
  n <- 50; S <- 18
  coords <- data.frame(x = runif(n), y = runif(n))
  ab <- matrix(rpois(n * S, 1.5), n, S)
  ab <- ab[, colSums(ab) > 0, drop = FALSE]; S <- ncol(ab)

  h <- spaccHill(ab, coords, q = c(0, 1, 2), n_seeds = 4, progress = FALSE)
  end <- vapply(h$curves, function(m) mean(m[, ncol(m)]), numeric(1))

  p <- colSums(ab) / sum(ab)
  hill_q <- function(p, q) if (q == 1) exp(-sum(p * log(p))) else sum(p^q)^(1 / (1 - q))

  expect_equal(unname(end[1]), as.numeric(S))
  expect_equal(unname(end[2]), hill_q(p, 1), tolerance = 1e-6)
  expect_equal(unname(end[3]), hill_q(p, 2), tolerance = 1e-6)
})


# ---- Phylogenetic endpoints equal mpd / mntd of the pool -----------------
test_that("spaccPhylo endpoints equal mpd / mntd of the observed species", {
  skip_if_not_installed("ape")
  set.seed(808)
  n <- 40; S <- 14
  coords <- data.frame(x = runif(n), y = runif(n))
  ab <- matrix(rbinom(n * S, 1, 0.3), n, S)
  ab <- ab[, colSums(ab) > 0, drop = FALSE]; S <- ncol(ab)
  colnames(ab) <- paste0("sp", seq_len(S))

  tr <- ape::rtree(S); tr$tip.label <- colnames(ab)
  cph <- ape::cophenetic.phylo(tr)[colnames(ab), colnames(ab)]

  ph <- spaccPhylo(ab, coords, tr, metric = c("mpd", "mntd"),
                   n_seeds = 4, progress = FALSE)
  end <- vapply(ph$curves, function(m) mean(m[, ncol(m)]), numeric(1))

  off <- cph; diag(off) <- NA
  mpd_def <- mean(off, na.rm = TRUE)
  mntd_def <- mean(apply(off, 1, min, na.rm = TRUE))

  expect_equal(unname(end[1]), mpd_def, tolerance = 1e-6)
  expect_equal(unname(end[2]), mntd_def, tolerance = 1e-6)
})


# ---- Functional endpoints equal Rao Q / FDis of the pool -----------------
test_that("spaccFunc endpoints equal abundance-weighted Rao Q and FDis", {
  set.seed(909)
  n <- 40; S <- 14
  coords <- data.frame(x = runif(n), y = runif(n))
  ab <- matrix(rpois(n * S, 2), n, S)
  ab <- ab[, colSums(ab) > 0, drop = FALSE]; S <- ncol(ab)
  colnames(ab) <- paste0("sp", seq_len(S))
  traits <- matrix(rnorm(S * 3), S, 3,
                   dimnames = list(colnames(ab), c("t1", "t2", "t3")))

  fu <- spaccFunc(ab, coords, traits, metric = c("rao", "fdis"),
                  n_seeds = 4, progress = FALSE)
  end <- vapply(fu$curves, function(m) mean(m[, ncol(m)]), numeric(1))

  w <- colSums(ab) / sum(ab)
  dtr <- as.matrix(dist(traits))
  rao_def <- as.numeric(w %*% dtr %*% w)
  cen <- colSums(w * traits)
  fdis_def <- sum(w * sqrt(rowSums((traits - matrix(cen, S, ncol(traits), byrow = TRUE))^2)))

  expect_equal(unname(end[1]), rao_def, tolerance = 1e-6)
  expect_equal(unname(end[2]), fdis_def, tolerance = 1e-6)
})
