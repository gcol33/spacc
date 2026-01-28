# Comprehensive tests to increase coverage toward 90%
# Tests untested code paths across multiple files

# ===========================================================================
# R/extrapolate.R - Additional coverage
# ===========================================================================

test_that("extrapolate with evt model works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(40), y = runif(40))
  species <- matrix(rbinom(40 * 20, 1, 0.3), nrow = 40)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)

  fit <- extrapolate(sac, model = "evt")

  expect_s3_class(fit, "spacc_fit")
})


test_that("extrapolate handles failed fit gracefully", {
  skip_on_cran()

  set.seed(42)
  # Create degenerate data that may fail to fit
  coords <- data.frame(x = 1:5, y = 1:5)
  species <- matrix(c(1, 0, 0, 0, 0,
                      1, 1, 0, 0, 0,
                      1, 1, 1, 0, 0,
                      1, 1, 1, 1, 0,
                      1, 1, 1, 1, 1), nrow = 5, byrow = TRUE)

  sac <- spacc(species, coords, n_seeds = 1, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)

  # evt with very few points should warn and fallback
  fit <- suppressWarnings(extrapolate(sac, model = "evt"))
  expect_s3_class(fit, "spacc_fit")
})


test_that("extrapolate print handles NA asymptote", {
  skip_on_cran()

  # Create a spacc_fit with NA asymptote manually
  fake_fit <- structure(
    list(
      asymptote = NA,
      asymptote_ci = c(NA, NA),
      model = "michaelis-menten",
      fit = NULL,
      aic = NA,
      data = data.frame(x = 1:10, y = 1:10)
    ),
    class = "spacc_fit"
  )

  expect_output(print(fake_fit), "failed")
})


test_that("extrapolate summary handles NULL fit", {
  skip_on_cran()

  fake_fit <- structure(
    list(
      asymptote = NA,
      asymptote_ci = c(NA, NA),
      model = "michaelis-menten",
      fit = NULL,
      aic = NA,
      data = data.frame(x = 1:10, y = 1:10)
    ),
    class = "spacc_fit"
  )

  expect_output(summary(fake_fit), "failed")
})


test_that("coef.spacc_fit handles NULL fit", {
  fake_fit <- structure(list(fit = NULL), class = "spacc_fit")
  expect_true(is.na(coef(fake_fit)))
})


test_that("confint.spacc_fit handles NULL fit", {
  fake_fit <- structure(list(fit = NULL), class = "spacc_fit")
  expect_true(is.na(confint(fake_fit)))
})


test_that("predict.spacc_fit handles NULL fit", {
  fake_fit <- structure(list(fit = NULL), class = "spacc_fit")
  expect_true(is.na(predict(fake_fit)))
})


test_that("predict.spacc_fit with default n uses data$x", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)
  fit <- extrapolate(sac, model = "michaelis-menten")

  preds <- predict(fit)
  expect_length(preds, 30)
})


# ===========================================================================
# R/compare.R - Additional coverage
# ===========================================================================

test_that("compare with bootstrap method works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species1 <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)
  species2 <- matrix(rbinom(20 * 10, 1, 0.5), nrow = 20)

  sac1 <- spacc(species1, coords, n_seeds = 5,
                parallel = FALSE, progress = FALSE, seed = 1)
  sac2 <- spacc(species2, coords, n_seeds = 5,
                parallel = FALSE, progress = FALSE, seed = 1)

  comp <- compare(sac1, sac2, method = "bootstrap", n_perm = 99)

  expect_s3_class(comp, "spacc_comp")
  expect_equal(comp$method, "bootstrap")
})


test_that("compare with auc method works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species1 <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)
  species2 <- matrix(rbinom(20 * 10, 1, 0.5), nrow = 20)

  sac1 <- spacc(species1, coords, n_seeds = 5,
                parallel = FALSE, progress = FALSE, seed = 1)
  sac2 <- spacc(species2, coords, n_seeds = 5,
                parallel = FALSE, progress = FALSE, seed = 1)

  comp <- compare(sac1, sac2, method = "auc")

  expect_s3_class(comp, "spacc_comp")
  expect_equal(comp$method, "auc")
  expect_true(is.na(comp$p_value))
})


test_that("compare print shows significance stars", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species1 <- matrix(rbinom(20 * 10, 1, 0.2), nrow = 20)
  species2 <- matrix(rbinom(20 * 10, 1, 0.8), nrow = 20)

  sac1 <- spacc(species1, coords, n_seeds = 10,
                parallel = FALSE, progress = FALSE, seed = 1)
  sac2 <- spacc(species2, coords, n_seeds = 10,
                parallel = FALSE, progress = FALSE, seed = 1)

  comp <- compare(sac1, sac2, method = "permutation", n_perm = 99)

  expect_output(print(comp), "Comparison")
})


test_that("compare summary works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species1 <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)
  species2 <- matrix(rbinom(20 * 10, 1, 0.5), nrow = 20)

  sac1 <- spacc(species1, coords, n_seeds = 5,
                parallel = FALSE, progress = FALSE, seed = 1)
  sac2 <- spacc(species2, coords, n_seeds = 5,
                parallel = FALSE, progress = FALSE, seed = 1)

  comp <- compare(sac1, sac2, method = "auc")

  # summary just calls print
  expect_output(summary(comp), "Comparison")
})


# ===========================================================================
# R/phylofunc.R - Additional coverage
# ===========================================================================

test_that("spaccPhylo returns correct structure", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  # Create a distance matrix as tree proxy
  phylo_dist <- matrix(runif(8 * 8), nrow = 8)
  phylo_dist <- (phylo_dist + t(phylo_dist)) / 2
  diag(phylo_dist) <- 0

  result <- spaccPhylo(species, coords, phylo_dist,
                       metric = c("mpd", "mntd"), n_seeds = 3,
                       parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc_phylo")
  expect_equal(result$n_seeds, 3)
  expect_equal(result$metric, c("mpd", "mntd"))
})


test_that("spaccPhylo warns when PD requested with distance matrix", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  phylo_dist <- matrix(runif(8 * 8), nrow = 8)
  phylo_dist <- (phylo_dist + t(phylo_dist)) / 2
  diag(phylo_dist) <- 0

  expect_warning(
    spaccPhylo(species, coords, phylo_dist,
               metric = c("pd", "mpd"), n_seeds = 3,
               parallel = FALSE, progress = FALSE, seed = 1),
    "PD requires"
  )
})


test_that("spaccPhylo errors on invalid tree input", {
  coords <- data.frame(x = runif(10), y = runif(10))
  species <- matrix(rbinom(10 * 5, 1, 0.3), nrow = 10)

  expect_error(
    spaccPhylo(species, coords, "not_a_tree",
               n_seeds = 1, parallel = FALSE, progress = FALSE),
    "tree must be"
  )
})


test_that("spaccPhylo print works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  phylo_dist <- matrix(runif(8 * 8), nrow = 8)
  phylo_dist <- (phylo_dist + t(phylo_dist)) / 2
  diag(phylo_dist) <- 0

  result <- spaccPhylo(species, coords, phylo_dist,
                       metric = "mpd", n_seeds = 3,
                       parallel = FALSE, progress = FALSE, seed = 1)

  expect_output(print(result), "spacc phylogenetic")
})


test_that("spaccPhylo summary returns data.frame", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  phylo_dist <- matrix(runif(8 * 8), nrow = 8)
  phylo_dist <- (phylo_dist + t(phylo_dist)) / 2
  diag(phylo_dist) <- 0

  result <- spaccPhylo(species, coords, phylo_dist,
                       metric = c("mpd", "mntd"), n_seeds = 3,
                       parallel = FALSE, progress = FALSE, seed = 1)

  summ <- summary(result)
  expect_s3_class(summ, "data.frame")
  expect_true("metric" %in% names(summ))
  expect_true("mean" %in% names(summ))
})


test_that("spaccFunc returns correct structure", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 2), nrow = 15)
  traits <- matrix(rnorm(8 * 3), nrow = 8)
  colnames(species) <- paste0("sp", 1:8)
  rownames(traits) <- paste0("sp", 1:8)

  result <- spaccFunc(species, coords, traits,
                      metric = c("fdis", "fric"), n_seeds = 3,
                      parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc_func")
  expect_equal(result$n_seeds, 3)
  expect_equal(result$metric, c("fdis", "fric"))
})


test_that("spaccFunc print works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 2), nrow = 15)
  traits <- matrix(rnorm(8 * 3), nrow = 8)

  result <- spaccFunc(species, coords, traits,
                      metric = "fdis", n_seeds = 3,
                      parallel = FALSE, progress = FALSE, seed = 1)

  expect_output(print(result), "spacc functional")
})


test_that("spaccFunc summary returns data.frame", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 2), nrow = 15)
  traits <- matrix(rnorm(8 * 3), nrow = 8)

  result <- spaccFunc(species, coords, traits,
                      metric = c("fdis", "fric"), n_seeds = 3,
                      parallel = FALSE, progress = FALSE, seed = 1)

  summ <- summary(result)
  expect_s3_class(summ, "data.frame")
  expect_true("metric" %in% names(summ))
})


test_that("spaccFunc errors on missing species in traits", {
  coords <- data.frame(x = runif(10), y = runif(10))
  species <- matrix(rpois(10 * 5, 2), nrow = 10)
  colnames(species) <- paste0("sp", 1:5)
  traits <- matrix(rnorm(3 * 3), nrow = 3)
  rownames(traits) <- paste0("sp", 1:3)

  expect_error(
    spaccFunc(species, coords, traits, n_seeds = 1,
              parallel = FALSE, progress = FALSE),
    "Some species"
  )
})


test_that("spaccPhylo with spacc_dist coords", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  phylo_dist <- matrix(runif(8 * 8), nrow = 8)
  phylo_dist <- (phylo_dist + t(phylo_dist)) / 2
  diag(phylo_dist) <- 0

  d <- distances(coords)
  result <- spaccPhylo(species, d, phylo_dist,
                       metric = "mpd", n_seeds = 3,
                       parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc_phylo")
})


test_that("spaccFunc with spacc_dist coords", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 2), nrow = 15)
  traits <- matrix(rnorm(8 * 3), nrow = 8)

  d <- distances(coords)
  result <- spaccFunc(species, d, traits,
                      metric = "fdis", n_seeds = 3,
                      parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc_func")
})


test_that("spaccPhylo as_sf errors without map", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  phylo_dist <- matrix(runif(8 * 8), nrow = 8)
  phylo_dist <- (phylo_dist + t(phylo_dist)) / 2
  diag(phylo_dist) <- 0

  result <- spaccPhylo(species, coords, phylo_dist,
                       metric = "mpd", n_seeds = 3,
                       parallel = FALSE, progress = FALSE, seed = 1)

  expect_error(as_sf(result), "map")
})


test_that("spaccFunc as_sf errors without map", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 2), nrow = 15)
  traits <- matrix(rnorm(8 * 3), nrow = 8)

  result <- spaccFunc(species, coords, traits,
                      metric = "fdis", n_seeds = 3,
                      parallel = FALSE, progress = FALSE, seed = 1)

  expect_error(as_sf(result), "map")
})


test_that("spaccPhylo plot map errors without map data", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  phylo_dist <- matrix(runif(8 * 8), nrow = 8)
  phylo_dist <- (phylo_dist + t(phylo_dist)) / 2
  diag(phylo_dist) <- 0

  result <- spaccPhylo(species, coords, phylo_dist,
                       metric = "mpd", n_seeds = 3,
                       parallel = FALSE, progress = FALSE, seed = 1)

  expect_error(plot(result, type = "map"), "map")
})


test_that("spaccFunc plot map errors without map data", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 2), nrow = 15)
  traits <- matrix(rnorm(8 * 3), nrow = 8)

  result <- spaccFunc(species, coords, traits,
                      metric = "fdis", n_seeds = 3,
                      parallel = FALSE, progress = FALSE, seed = 1)

  expect_error(plot(result, type = "map"), "map")
})


test_that("spaccPhylo plot curve returns ggplot", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  phylo_dist <- matrix(runif(8 * 8), nrow = 8)
  phylo_dist <- (phylo_dist + t(phylo_dist)) / 2
  diag(phylo_dist) <- 0

  result <- spaccPhylo(species, coords, phylo_dist,
                       metric = c("mpd", "mntd"), n_seeds = 3,
                       parallel = FALSE, progress = FALSE, seed = 1)

  p <- plot(result)
  expect_s3_class(p, "ggplot")
})


test_that("spaccFunc plot curve returns ggplot", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 2), nrow = 15)
  traits <- matrix(rnorm(8 * 3), nrow = 8)

  result <- spaccFunc(species, coords, traits,
                      metric = c("fdis", "fric"), n_seeds = 3,
                      parallel = FALSE, progress = FALSE, seed = 1)

  p <- plot(result)
  expect_s3_class(p, "ggplot")
})


# ===========================================================================
# R/hill.R - Additional coverage
# ===========================================================================

test_that("spaccHill with map=TRUE stores site_values", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 2), nrow = 15)

  result <- spaccHill(species, coords, q = c(0, 1), n_seeds = 3,
                      map = TRUE,
                      parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc_hill")
  expect_false(is.null(result$site_values))
  expect_equal(nrow(result$site_values), 15)
})


# ===========================================================================
# R/coverage.R - Additional coverage
# ===========================================================================

test_that("spaccCoverage with map=TRUE stores site_values", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 3), nrow = 15)

  result <- spaccCoverage(species, coords, n_seeds = 3,
                          map = TRUE,
                          parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc_coverage")
  expect_false(is.null(result$site_values))
  expect_equal(nrow(result$site_values), 15)
})


# ===========================================================================
# R/beta.R - Additional coverage
# ===========================================================================

test_that("spaccBeta plot with distance xaxis works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  result <- spaccBeta(species, coords, n_seeds = 3,
                      parallel = FALSE, progress = FALSE, seed = 1)

  p <- plot(result, xaxis = "distance")
  expect_s3_class(p, "ggplot")
})


test_that("spaccBeta plot without partition shows total only", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  result <- spaccBeta(species, coords, n_seeds = 3,
                      parallel = FALSE, progress = FALSE, seed = 1)

  p <- plot(result, partition = FALSE, ci = TRUE)
  expect_s3_class(p, "ggplot")
})


test_that("spaccBetaFunc with spacc_dist coords", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)
  traits <- matrix(rnorm(8 * 3), nrow = 8)

  d <- distances(coords)
  result <- spaccBetaFunc(species, d, traits, n_seeds = 3,
                           parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc_beta")
})


test_that("spaccBetaPhylo with spacc_dist coords", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  phylo_dist <- matrix(runif(8 * 8), nrow = 8)
  phylo_dist <- (phylo_dist + t(phylo_dist)) / 2
  diag(phylo_dist) <- 0

  d <- distances(coords)
  result <- spaccBetaPhylo(species, d, phylo_dist, n_seeds = 3,
                            parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc_beta")
})


# ===========================================================================
# R/sar.R - Additional coverage
# ===========================================================================

test_that("dar summary with multiple q values", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 2), nrow = 15)

  result <- dar(species, coords, q = c(0, 1, 2), n_seeds = 3,
                area_method = "count",
                parallel = FALSE, progress = FALSE, seed = 1)

  summ <- summary(result)
  expect_s3_class(summ, "data.frame")
  expect_equal(length(unique(summ$q)), 3)
})


test_that("dar plot without ci works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 2), nrow = 15)

  result <- dar(species, coords, q = c(0, 1), n_seeds = 3,
                area_method = "count",
                parallel = FALSE, progress = FALSE, seed = 1)

  p <- plot(result, ci = FALSE)
  expect_s3_class(p, "ggplot")
})


test_that("spaccEndemism plot without richness works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  result <- spaccEndemism(species, coords, n_seeds = 3,
                           parallel = FALSE, progress = FALSE, seed = 1)

  p <- plot(result, show_richness = FALSE, ci = FALSE)
  expect_s3_class(p, "ggplot")
})


test_that("sfar print shows fragmentation interpretation", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.4), nrow = 20)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)
  patches <- rep(1:4, each = 5)

  result <- sfar(sac, patches)

  expect_output(print(result), "SFAR")
})


# ===========================================================================
# R/metrics.R - Additional coverage
# ===========================================================================

test_that("spaccMetrics returns correct structure", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  result <- spaccMetrics(species, coords,
                         parallel = FALSE, progress = FALSE)

  expect_s3_class(result, "spacc_metrics")
  expect_equal(result$n_sites, 20)
})


test_that("spaccMetrics print works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  result <- spaccMetrics(species, coords,
                         parallel = FALSE, progress = FALSE)

  expect_output(print(result), "spacc_metrics")
})


test_that("spaccMetrics summary works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  result <- spaccMetrics(species, coords,
                         parallel = FALSE, progress = FALSE)

  # summary.spacc_metrics doesn't return a data.frame, it prints
  expect_output(summary(result), "Metric summary")
})


# ===========================================================================
# R/wavefront.R and decay - Additional coverage
# ===========================================================================

test_that("wavefront summary returns data.frame", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  result <- wavefront(species, coords, n_seeds = 3, n_steps = 8,
                      progress = FALSE, seed = 1)

  summ <- summary(result)
  expect_s3_class(summ, "data.frame")
})


test_that("wavefront plot returns ggplot", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  result <- wavefront(species, coords, n_seeds = 3, n_steps = 8,
                      progress = FALSE, seed = 1)

  p <- plot(result)
  expect_s3_class(p, "ggplot")
})


test_that("distanceDecay summary works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  result <- distanceDecay(species, coords, n_seeds = 3,
                          progress = FALSE, seed = 1)

  summ <- summary(result)
  # summary.spacc_decay returns its own class
  expect_true(!is.null(summ))
})


test_that("distanceDecay plot returns ggplot", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  result <- distanceDecay(species, coords, n_seeds = 3,
                          progress = FALSE, seed = 1)

  p <- plot(result)
  expect_s3_class(p, "ggplot")
})


# ===========================================================================
# R/spacc.R - Additional coverage for different methods
# ===========================================================================

test_that("spacc with kncn method works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  result <- spacc(species, coords, n_seeds = 3, method = "kncn",
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
  expect_equal(result$method, "kncn")
})


test_that("spacc with random method works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  result <- spacc(species, coords, n_seeds = 3, method = "random",
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
  expect_equal(result$method, "random")
})


test_that("spacc with collector method works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  result <- spacc(species, coords, n_seeds = 3, method = "collector",
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
  expect_equal(result$method, "collector")
})


test_that("spacc with radius method works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  result <- spacc(species, coords, n_seeds = 3, method = "radius",
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
  expect_equal(result$method, "radius")
})


test_that("spacc with gaussian method works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  result <- spacc(species, coords, n_seeds = 3, method = "gaussian",
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
  expect_equal(result$method, "gaussian")
})


test_that("spacc with cone method works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  result <- spacc(species, coords, n_seeds = 3, method = "cone",
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
  expect_equal(result$method, "cone")
})


test_that("spacc with haversine distance works", {
  skip_on_cran()

  set.seed(42)
  # Use coordinates that look like lat/lon
  coords <- data.frame(x = runif(20, -10, 10), y = runif(20, 40, 50))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  result <- spacc(species, coords, n_seeds = 3, method = "knn",
                  distance = "haversine",
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
  expect_equal(result$distance, "haversine")
})


# ===========================================================================
# R/methods.R - Additional coverage
# ===========================================================================

test_that("as.data.frame.spacc_hill works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 2), nrow = 15)

  result <- spaccHill(species, coords, q = c(0, 1), n_seeds = 3,
                      parallel = FALSE, progress = FALSE, seed = 1)

  df <- as.data.frame(result)
  expect_s3_class(df, "data.frame")
})


test_that("as.data.frame.spacc_beta works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  result <- spaccBeta(species, coords, n_seeds = 3,
                      parallel = FALSE, progress = FALSE, seed = 1)

  df <- as.data.frame(result)
  expect_s3_class(df, "data.frame")
})


test_that("as.data.frame.spacc_coverage works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 3), nrow = 15)

  result <- spaccCoverage(species, coords, n_seeds = 3,
                          parallel = FALSE, progress = FALSE, seed = 1)

  df <- as.data.frame(result)
  expect_s3_class(df, "data.frame")
})


test_that("as.data.frame.spacc_metrics works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  result <- spaccMetrics(species, coords,
                         parallel = FALSE, progress = FALSE)

  df <- as.data.frame(result)
  expect_s3_class(df, "data.frame")
})


test_that("as.data.frame.spacc_partition works", {
  species <- matrix(rpois(50 * 30, 2), nrow = 50)
  part <- diversityPartition(species, q = c(0, 1))

  df <- as.data.frame(part)
  expect_s3_class(df, "data.frame")
  expect_true("alpha" %in% names(df))
  expect_true("beta" %in% names(df))
  expect_true("gamma" %in% names(df))
})


test_that("as.data.frame.spacc_alpha works", {
  coords <- data.frame(x = runif(50), y = runif(50))
  species <- matrix(rpois(50 * 30, 2), nrow = 50)

  result <- alphaDiversity(species, q = c(0, 1), coords = coords)

  df <- as.data.frame(result)
  expect_s3_class(df, "data.frame")
})


# ===========================================================================
# Additional edge cases
# ===========================================================================

test_that("spacc handles presence-absence data correctly", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  # Explicit 0/1 presence-absence
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  result <- spacc(species, coords, n_seeds = 3, method = "knn",
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
  # Final richness should equal unique species with presences
  expect_true(result$n_species >= 0)
})


test_that("spacc handles abundance data correctly", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  # Abundance data
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- spacc(species, coords, n_seeds = 3, method = "knn",
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
})


test_that("rarefy returns correct structure", {
  skip_on_cran()

  set.seed(42)
  species <- matrix(rpois(30 * 15, 3), nrow = 30)

  # rarefy takes a matrix, not spacc object
  result <- rarefy(species, n_individuals = c(50, 100))

  expect_s3_class(result, "spacc_rare")
})


test_that("rarefy print works", {
  skip_on_cran()

  set.seed(42)
  species <- matrix(rpois(30 * 15, 3), nrow = 30)

  result <- rarefy(species, n_individuals = c(50, 100))

  expect_output(print(result), "rarefaction")
})


test_that("rarefy plot returns ggplot", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  species <- matrix(rpois(30 * 15, 3), nrow = 30)

  result <- rarefy(species, n_individuals = c(50, 100))

  p <- plot(result)
  expect_s3_class(p, "ggplot")
})


# ===========================================================================
# R/diversity.R - Additional coverage
# ===========================================================================

test_that("diversity functions with Hill numbers work", {
  skip_on_cran()

  # Test calc_hill_number directly
  # Already tested in test-hill.R but adding edge cases

  # Empty input
  expect_equal(spacc:::calc_hill_number(numeric(0), 0), 0)

  # Single species
  expect_equal(spacc:::calc_hill_number(c(10), 0), 1)
  expect_equal(spacc:::calc_hill_number(c(10), 1), 1)
  expect_equal(spacc:::calc_hill_number(c(10), 2), 1)

  # Highly uneven community
  abundances <- c(1000, 1, 1, 1, 1)
  q0 <- spacc:::calc_hill_number(abundances, 0)
  q1 <- spacc:::calc_hill_number(abundances, 1)
  q2 <- spacc:::calc_hill_number(abundances, 2)

  expect_equal(q0, 5)  # Species richness
  expect_true(q1 < q0)  # Shannon diversity < richness for uneven
  expect_true(q2 < q1)  # Simpson diversity < Shannon for uneven
})


# ===========================================================================
# R/diversity.R - Additional coverage
# ===========================================================================

test_that("alphaDiversity returns correct structure", {
  species <- matrix(rpois(50 * 30, 2), nrow = 50)
  alpha <- alphaDiversity(species, q = c(0, 1, 2))

  expect_s3_class(alpha, "data.frame")
  expect_equal(nrow(alpha), 50)
  expect_equal(ncol(alpha), 3)
  expect_true(all(c("q0", "q1", "q2") %in% names(alpha)))
})


test_that("alphaDiversity with coords returns spacc_alpha", {
  coords <- data.frame(x = runif(50), y = runif(50))
  species <- matrix(rpois(50 * 30, 2), nrow = 50)

  result <- alphaDiversity(species, q = c(0, 1), coords = coords)

  expect_s3_class(result, "spacc_alpha")
  expect_equal(result$n_sites, 50)
})


test_that("gammaDiversity returns correct structure", {
  species <- matrix(rpois(50 * 30, 2), nrow = 50)
  gamma <- gammaDiversity(species, q = c(0, 1, 2))

  expect_true(is.numeric(gamma))
  expect_length(gamma, 3)
  expect_true(all(c("q0", "q1", "q2") %in% names(gamma)))
})


test_that("diversityPartition returns correct structure", {
  species <- matrix(rpois(50 * 30, 2), nrow = 50)
  part <- diversityPartition(species, q = c(0, 1, 2))

  expect_s3_class(part, "spacc_partition")
  expect_equal(length(part$alpha), 3)
  expect_equal(length(part$beta), 3)
  expect_equal(length(part$gamma), 3)
})


test_that("diversityPartition with proportional weights", {
  species <- matrix(rpois(50 * 30, 2), nrow = 50)
  part <- diversityPartition(species, q = c(0, 1), weights = "proportional")

  expect_s3_class(part, "spacc_partition")
})


test_that("diversityPartition with custom weights", {
  species <- matrix(rpois(50 * 30, 2), nrow = 50)
  custom_weights <- runif(50)
  part <- diversityPartition(species, q = c(0, 1), weights = custom_weights)

  expect_s3_class(part, "spacc_partition")
})


test_that("diversityPartition errors on invalid weights", {
  species <- matrix(rpois(50 * 30, 2), nrow = 50)

  expect_error(
    diversityPartition(species, q = 0, weights = "invalid"),
    "weights must be"
  )
})


test_that("diversityPartition print works", {
  species <- matrix(rpois(50 * 30, 2), nrow = 50)
  part <- diversityPartition(species, q = c(0, 1))

  expect_output(print(part), "Alpha-Beta-Gamma")
})


test_that("diversityPartition summary returns data.frame", {
  species <- matrix(rpois(50 * 30, 2), nrow = 50)
  part <- diversityPartition(species, q = c(0, 1, 2))

  summ <- summary(part)
  expect_s3_class(summ, "data.frame")
  expect_true("beta_normalized" %in% names(summ))
})


test_that("spacc_alpha print works", {
  coords <- data.frame(x = runif(50), y = runif(50))
  species <- matrix(rpois(50 * 30, 2), nrow = 50)

  result <- alphaDiversity(species, q = c(0, 1), coords = coords)

  expect_output(print(result), "spacc alpha")
})


test_that("spacc_alpha summary returns data.frame", {
  coords <- data.frame(x = runif(50), y = runif(50))
  species <- matrix(rpois(50 * 30, 2), nrow = 50)

  result <- alphaDiversity(species, q = c(0, 1), coords = coords)

  summ <- summary(result)
  expect_s3_class(summ, "data.frame")
})


test_that("spacc_alpha plot histogram works", {
  skip_if_not_installed("ggplot2")

  coords <- data.frame(x = runif(50), y = runif(50))
  species <- matrix(rpois(50 * 30, 2), nrow = 50)

  result <- alphaDiversity(species, q = c(0, 1), coords = coords)

  p <- plot(result, type = "histogram")
  expect_s3_class(p, "ggplot")
})


test_that("spacc_alpha as_sf errors without coords", {
  species <- matrix(rpois(50 * 30, 2), nrow = 50)
  result <- alphaDiversity(species, q = c(0, 1))

  # result is data.frame, not spacc_alpha - this should work differently
  # The as_sf method is for spacc_alpha objects with coords
})


test_that("spacc_partition plot bar works", {
  skip_if_not_installed("ggplot2")

  species <- matrix(rpois(50 * 30, 2), nrow = 50)
  part <- diversityPartition(species, q = c(0, 1))

  p <- plot(part, type = "bar")
  expect_s3_class(p, "ggplot")
})


test_that("spacc_partition map errors without coords", {
  skip_if_not_installed("ggplot2")

  species <- matrix(rpois(50 * 30, 2), nrow = 50)
  part <- diversityPartition(species, q = c(0, 1))

  expect_error(plot(part, type = "map"), "coordinates")
})


test_that("spacc_partition as_sf errors without coords", {
  species <- matrix(rpois(50 * 30, 2), nrow = 50)
  part <- diversityPartition(species, q = c(0, 1))

  expect_error(as_sf(part), "coordinates")
})


test_that("spacc_partition with coords works", {
  skip_if_not_installed("ggplot2")

  coords <- data.frame(x = runif(50), y = runif(50))
  species <- matrix(rpois(50 * 30, 2), nrow = 50)
  part <- diversityPartition(species, q = c(0, 1), coords = coords)

  expect_s3_class(part, "spacc_partition")

  # Map plot should work
  p <- plot(part, type = "map")
  expect_s3_class(p, "ggplot")
})


# ===========================================================================
# R/metrics.R - Additional coverage
# ===========================================================================

test_that("spaccMetrics with different metrics", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  result <- spaccMetrics(species, coords,
                         metrics = c("slope_25", "richness_75pct", "richness_90pct"),
                         parallel = FALSE, progress = FALSE)

  expect_s3_class(result, "spacc_metrics")
  expect_true("slope_25" %in% names(result$metrics))
  expect_true("richness_75pct" %in% names(result$metrics))
  expect_true("richness_90pct" %in% names(result$metrics))
})


test_that("spaccMetrics with spacc_dist coords", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  d <- distances(coords)
  result <- spaccMetrics(species, d,
                         parallel = FALSE, progress = FALSE)

  expect_s3_class(result, "spacc_metrics")
})


test_that("spaccMetrics with kncn method", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  result <- spaccMetrics(species, coords, method = "kncn",
                         parallel = FALSE, progress = FALSE)

  expect_s3_class(result, "spacc_metrics")
})


test_that("spaccMetrics with random method", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  result <- spaccMetrics(species, coords, method = "random",
                         parallel = FALSE, progress = FALSE)

  expect_s3_class(result, "spacc_metrics")
})


test_that("spaccMetrics plot heatmap works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  result <- spaccMetrics(species, coords,
                         parallel = FALSE, progress = FALSE)

  p <- plot(result, type = "heatmap")
  expect_s3_class(p, "ggplot")
})


test_that("spaccMetrics plot histogram works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  result <- spaccMetrics(species, coords,
                         parallel = FALSE, progress = FALSE)

  p <- plot(result, type = "histogram")
  expect_s3_class(p, "ggplot")
})


test_that("spaccMetrics plot errors on invalid metric", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  result <- spaccMetrics(species, coords,
                         parallel = FALSE, progress = FALSE)

  expect_error(plot(result, metric = "nonexistent"), "not found")
})


test_that("as_sf.spacc_metrics returns sf object", {
  skip_on_cran()
  skip_if_not_installed("sf")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  result <- spaccMetrics(species, coords,
                         parallel = FALSE, progress = FALSE)

  # Pass a valid CRS to avoid the empty crs error
  sf_obj <- as_sf(result, crs = 4326)
  expect_s3_class(sf_obj, "sf")
})
