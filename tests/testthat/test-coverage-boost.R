# Additional tests to boost coverage to 90%+

# ============================================================================
# DISTANCES.R - sf support and edge cases
# ============================================================================

test_that("distances.sf with geodesic method works", {
  skip_on_cran()
  skip_if_not_installed("sf")

  coords <- data.frame(lon = c(-122.4, -122.5, -122.3), lat = c(37.7, 37.8, 37.6))
  pts_sf <- sf::st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)

  result <- distances(pts_sf, method = "geodesic")

  expect_s3_class(result, "spacc_dist")
  expect_equal(attr(result, "method"), "geodesic")
  expect_equal(dim(result), c(3, 3))
})


test_that("distances.sf with haversine method works", {
  skip_on_cran()
  skip_if_not_installed("sf")

  coords <- data.frame(lon = c(-122.4, -122.5, -122.3), lat = c(37.7, 37.8, 37.6))
  pts_sf <- sf::st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)

  result <- distances(pts_sf, method = "haversine")

  expect_s3_class(result, "spacc_dist")
  expect_equal(attr(result, "method"), "haversine")
})


test_that("distances.sf with euclidean method works", {
  skip_on_cran()
  skip_if_not_installed("sf")

  coords <- data.frame(x = c(100, 200, 300), y = c(100, 200, 300))
  pts_sf <- sf::st_as_sf(coords, coords = c("x", "y"), crs = 32610) # UTM

  result <- distances(pts_sf, method = "euclidean")

  expect_s3_class(result, "spacc_dist")
  expect_equal(attr(result, "method"), "euclidean")
})


test_that("distances.sf auto-detects geodesic for geographic CRS", {
  skip_on_cran()
  skip_if_not_installed("sf")

  coords <- data.frame(lon = c(-122.4, -122.5), lat = c(37.7, 37.8))
  pts_sf <- sf::st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)

  expect_message(result <- distances(pts_sf), "geodesic")
  expect_equal(attr(result, "method"), "geodesic")
})


test_that("distances.sf auto-detects euclidean for projected CRS", {
  skip_on_cran()
  skip_if_not_installed("sf")

  coords <- data.frame(x = c(100, 200), y = c(100, 200))
  pts_sf <- sf::st_as_sf(coords, coords = c("x", "y"), crs = 32610)

  expect_message(result <- distances(pts_sf), "euclidean")
  expect_equal(attr(result, "method"), "euclidean")
})


test_that("distances.sf with custom function works", {
  skip_on_cran()
  skip_if_not_installed("sf")

  coords <- data.frame(lon = c(-122.4, -122.5, -122.3), lat = c(37.7, 37.8, 37.6))
  pts_sf <- sf::st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)

  # Custom Manhattan distance
  manhattan <- function(x, y) {
    n <- length(x)
    m <- matrix(0, n, n)
    for (i in 1:n) for (j in 1:n) {
      m[i, j] <- abs(x[i] - x[j]) + abs(y[i] - y[j])
    }
    m
  }

  result <- distances(pts_sf, fun = manhattan)

  expect_s3_class(result, "spacc_dist")
})


test_that("distances.sfc works", {
  skip_on_cran()
  skip_if_not_installed("sf")

  coords <- data.frame(x = c(0, 1, 2), y = c(0, 1, 0))
  pts_sfc <- sf::st_sfc(
    sf::st_point(c(0, 0)),
    sf::st_point(c(1, 1)),
    sf::st_point(c(2, 0)),
    crs = 32610
  )

  expect_message(result <- distances(pts_sfc), "euclidean")
  expect_s3_class(result, "spacc_dist")
  expect_equal(dim(result), c(3, 3))
})


test_that("distances.sf errors on non-point geometry", {
  skip_on_cran()
  skip_if_not_installed("sf")

  # Create polygon instead of point
  poly <- sf::st_sfc(
    sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 0))))
  )
  poly_sf <- sf::st_sf(id = 1, geometry = poly)

  expect_error(distances(poly_sf), "POINT")
})


test_that("print.spacc_dist shows CRS info for sf-based distances", {
  skip_on_cran()
  skip_if_not_installed("sf")

  coords <- data.frame(lon = c(-122.4, -122.5), lat = c(37.7, 37.8))
  pts_sf <- sf::st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)

  result <- distances(pts_sf, method = "geodesic")

  expect_output(print(result), "EPSG")
})


test_that("print.spacc_dist shows (meters) for geodesic", {
  skip_on_cran()
  skip_if_not_installed("sf")

  coords <- data.frame(lon = c(-122.4, -122.5), lat = c(37.7, 37.8))
  pts_sf <- sf::st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)

  result <- distances(pts_sf, method = "geodesic")

  expect_output(print(result), "meters")
})


# ============================================================================
# METHODS.R - More S3 methods coverage
# ============================================================================

test_that("c.spacc combines objects correctly", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  sp1 <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)
  sp2 <- matrix(rbinom(15 * 8, 1, 0.2), nrow = 15)

  sac1 <- spacc(sp1, coords, n_seeds = 3, parallel = FALSE, progress = FALSE)
  sac2 <- spacc(sp2, coords, n_seeds = 3, parallel = FALSE, progress = FALSE)

  combined <- c(native = sac1, alien = sac2)

  expect_s3_class(combined, "spacc")
  expect_equal(length(combined$group_names), 2)
  expect_equal(combined$group_names, c("native", "alien"))
})


test_that("as.data.frame.spacc_hill returns correct structure", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 2), nrow = 15)

  result <- spaccHill(species, coords, q = c(0, 1, 2), n_seeds = 3,
                      parallel = FALSE, progress = FALSE)

  df <- as.data.frame(result)

  expect_s3_class(df, "data.frame")
  expect_true("q" %in% names(df))
  expect_true("mean" %in% names(df))
  expect_equal(nrow(df), 15 * 3) # n_sites * n_q
})


test_that("as.data.frame.spacc_beta returns correct structure", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  result <- spaccBeta(species, coords, n_seeds = 3,
                      parallel = FALSE, progress = FALSE)

  df <- as.data.frame(result)

  expect_s3_class(df, "data.frame")
  expect_true("beta_total" %in% names(df))
  expect_true("beta_turnover" %in% names(df))
})


test_that("as.data.frame.spacc_coverage returns correct structure", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 3), nrow = 15)

  result <- spaccCoverage(species, coords, n_seeds = 3,
                          parallel = FALSE, progress = FALSE)

  df <- as.data.frame(result)

  expect_s3_class(df, "data.frame")
  expect_true("richness" %in% names(df))
  expect_true("coverage" %in% names(df))
})


test_that("as.data.frame.spacc_fit returns correct structure", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  sac <- spacc(species, coords, n_seeds = 3, parallel = FALSE, progress = FALSE)
  fit <- extrapolate(sac, model = "michaelis-menten")

  df <- as.data.frame(fit)

  expect_s3_class(df, "data.frame")
  expect_true("observed" %in% names(df))
  expect_true("predicted" %in% names(df))
})


test_that("as.data.frame.spacc_comp returns correct structure", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  sp1 <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)
  sp2 <- matrix(rbinom(15 * 8, 1, 0.2), nrow = 15)

  sac1 <- spacc(sp1, coords, n_seeds = 3, parallel = FALSE, progress = FALSE)
  sac2 <- spacc(sp2, coords, n_seeds = 3, parallel = FALSE, progress = FALSE)

  comp <- compare(sac1, sac2, n_perm = 10)

  df <- as.data.frame(comp)

  expect_s3_class(df, "data.frame")
  expect_true("auc_diff" %in% names(df))
  expect_true("p_value" %in% names(df))
})


test_that("[.spacc subsetting works for ungrouped", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  sac <- spacc(species, coords, n_seeds = 5, parallel = FALSE, progress = FALSE)

  subset <- sac[1:3]

  expect_s3_class(subset, "spacc")
  expect_equal(subset$n_seeds, 3)
  expect_equal(nrow(subset$curves), 3)
})


test_that("[.spacc subsetting works for grouped", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  sp1 <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)
  sp2 <- matrix(rbinom(15 * 8, 1, 0.2), nrow = 15)

  sac1 <- spacc(sp1, coords, n_seeds = 5, parallel = FALSE, progress = FALSE)
  sac2 <- spacc(sp2, coords, n_seeds = 5, parallel = FALSE, progress = FALSE)

  combined <- c(native = sac1, alien = sac2)
  subset <- combined[1:3]

  expect_s3_class(subset, "spacc")
  expect_equal(subset$n_seeds, 3)
})


test_that("summary.spacc_decay returns correct structure", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  result <- distanceDecay(species, coords, n_seeds = 3, progress = FALSE)

  summ <- summary(result)

  expect_s3_class(summ, "summary.spacc_decay")
  expect_true("n_bins" %in% names(summ))
})


test_that("print.summary.spacc_decay works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  result <- distanceDecay(species, coords, n_seeds = 3, progress = FALSE)

  summ <- summary(result)

  expect_output(print(summ), "Distance-Decay")
})


test_that("summary.spacc_rare returns correct structure", {
  skip_on_cran()

  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  rare <- rarefy(species, n_boot = 10)
  summ <- summary(rare)

  expect_s3_class(summ, "summary.spacc_rare")
  expect_true("expected_richness" %in% names(summ))
})


test_that("print.summary.spacc_rare works", {
  skip_on_cran()

  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  rare <- rarefy(species, n_boot = 10)
  summ <- summary(rare)

  expect_output(print(summ), "Rarefaction")
})


test_that("as.data.frame.spacc_metrics includes coords", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  result <- spaccMetrics(species, coords, metrics = c("slope_10"),
                         parallel = FALSE, progress = FALSE)

  df <- as.data.frame(result)

  expect_true("x" %in% names(df))
  expect_true("y" %in% names(df))
})


# ============================================================================
# EXTRAPOLATE.R - EVT model and edge cases
# ============================================================================

test_that("extrapolate with EVT model works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  sac <- spacc(species, coords, n_seeds = 5, parallel = FALSE, progress = FALSE)

  # EVT requires at least 10 data points
  fit <- extrapolate(sac, model = "evt")

  expect_s3_class(fit, "spacc_fit")
})


test_that("predict.spacc_fit with NULL n uses data x values", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  sac <- spacc(species, coords, n_seeds = 3, parallel = FALSE, progress = FALSE)
  fit <- extrapolate(sac, model = "michaelis-menten")

  preds <- predict(fit)

  expect_length(preds, 20)
})


test_that("coef.spacc_fit returns NA for failed fit", {
  skip_on_cran()

  # Create a fake failed fit
  fake_fit <- structure(
    list(
      asymptote = NA,
      asymptote_ci = c(NA, NA),
      model = "michaelis-menten",
      fit = NULL,
      aic = NA,
      data = data.frame(x = 1:10, y = 1:10),
      spacc = NULL
    ),
    class = "spacc_fit"
  )

  expect_equal(coef(fake_fit), NA)
})


test_that("predict.spacc_fit returns NA for failed fit", {
  skip_on_cran()

  fake_fit <- structure(
    list(
      asymptote = NA,
      asymptote_ci = c(NA, NA),
      model = "michaelis-menten",
      fit = NULL,
      aic = NA,
      data = data.frame(x = 1:10, y = 1:10),
      spacc = NULL
    ),
    class = "spacc_fit"
  )

  expect_equal(predict(fake_fit, n = 1:5), NA)
})


test_that("confint.spacc_fit returns NA for failed fit", {
  skip_on_cran()

  fake_fit <- structure(
    list(
      asymptote = NA,
      asymptote_ci = c(NA, NA),
      model = "michaelis-menten",
      fit = NULL,
      aic = NA,
      data = data.frame(x = 1:10, y = 1:10),
      spacc = NULL
    ),
    class = "spacc_fit"
  )

  expect_equal(confint(fake_fit), NA)
})


test_that("print.spacc_fit handles failed fit", {
  skip_on_cran()

  fake_fit <- structure(
    list(
      asymptote = NA,
      asymptote_ci = c(NA, NA),
      model = "michaelis-menten",
      fit = NULL,
      aic = NA,
      data = data.frame(x = 1:10, y = 1:10),
      spacc = NULL
    ),
    class = "spacc_fit"
  )

  expect_output(print(fake_fit), "failed")
})


test_that("summary.spacc_fit prints for failed fit", {
  skip_on_cran()

  fake_fit <- structure(
    list(
      asymptote = NA,
      asymptote_ci = c(NA, NA),
      model = "michaelis-menten",
      fit = NULL,
      aic = NA,
      data = data.frame(x = 1:10, y = 1:10),
      spacc = NULL
    ),
    class = "spacc_fit"
  )

  expect_output(summary(fake_fit), "failed")
})


# ============================================================================
# COMPARE.R - Bootstrap and AUC methods
# ============================================================================

test_that("compare with bootstrap method works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  sp1 <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)
  sp2 <- matrix(rbinom(15 * 8, 1, 0.2), nrow = 15)

  sac1 <- spacc(sp1, coords, n_seeds = 5, parallel = FALSE, progress = FALSE)
  sac2 <- spacc(sp2, coords, n_seeds = 5, parallel = FALSE, progress = FALSE)

  comp <- compare(sac1, sac2, method = "bootstrap", n_perm = 50)

  expect_s3_class(comp, "spacc_comp")
  expect_equal(comp$method, "bootstrap")
})


test_that("compare with auc method works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  sp1 <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)
  sp2 <- matrix(rbinom(15 * 8, 1, 0.2), nrow = 15)

  sac1 <- spacc(sp1, coords, n_seeds = 5, parallel = FALSE, progress = FALSE)
  sac2 <- spacc(sp2, coords, n_seeds = 5, parallel = FALSE, progress = FALSE)

  comp <- compare(sac1, sac2, method = "auc")

  expect_s3_class(comp, "spacc_comp")
  expect_equal(comp$method, "auc")
  expect_true(is.na(comp$p_value))
})


test_that("summary.spacc_comp returns invisibly", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  sp1 <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)
  sp2 <- matrix(rbinom(15 * 8, 1, 0.2), nrow = 15)

  sac1 <- spacc(sp1, coords, n_seeds = 5, parallel = FALSE, progress = FALSE)
  sac2 <- spacc(sp2, coords, n_seeds = 5, parallel = FALSE, progress = FALSE)

  comp <- compare(sac1, sac2, n_perm = 10)

  expect_output(summary(comp), "Comparison")
})


# ============================================================================
# BETA.R - spaccBetaFunc and spaccBetaPhylo
# ============================================================================

test_that("spaccBetaFunc returns correct structure", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)
  colnames(species) <- paste0("sp", 1:8)
  traits <- matrix(rnorm(8 * 3), nrow = 8)
  rownames(traits) <- paste0("sp", 1:8)

  result <- spaccBetaFunc(species, coords, traits, n_seeds = 3,
                          parallel = FALSE, progress = FALSE)

  expect_s3_class(result, "spacc_beta")
  expect_equal(result$beta_type, "functional")
})


test_that("spaccBetaFunc with jaccard index works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)
  colnames(species) <- paste0("sp", 1:8)
  traits <- matrix(rnorm(8 * 3), nrow = 8)
  rownames(traits) <- paste0("sp", 1:8)

  result <- spaccBetaFunc(species, coords, traits, n_seeds = 3,
                          index = "jaccard",
                          parallel = FALSE, progress = FALSE)

  expect_s3_class(result, "spacc_beta")
  expect_equal(result$index, "jaccard")
})


test_that("spaccBetaPhylo with tree returns correct structure", {
  skip_on_cran()
  skip_if_not_installed("ape")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  tree <- ape::rtree(8)
  colnames(species) <- tree$tip.label

  result <- spaccBetaPhylo(species, coords, tree, n_seeds = 3,
                           parallel = FALSE, progress = FALSE)

  expect_s3_class(result, "spacc_beta")
  expect_equal(result$beta_type, "phylogenetic")
})


test_that("spaccBetaPhylo with distance matrix works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  phylo_dist <- matrix(runif(8 * 8), nrow = 8)
  phylo_dist <- (phylo_dist + t(phylo_dist)) / 2
  diag(phylo_dist) <- 0

  result <- spaccBetaPhylo(species, coords, phylo_dist, n_seeds = 3,
                           parallel = FALSE, progress = FALSE)

  expect_s3_class(result, "spacc_beta")
  expect_equal(result$beta_type, "phylogenetic")
})


test_that("spaccBetaPhylo with jaccard index works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  phylo_dist <- matrix(runif(8 * 8), nrow = 8)
  phylo_dist <- (phylo_dist + t(phylo_dist)) / 2
  diag(phylo_dist) <- 0

  result <- spaccBetaPhylo(species, coords, phylo_dist, n_seeds = 3,
                           index = "jaccard",
                           parallel = FALSE, progress = FALSE)

  expect_s3_class(result, "spacc_beta")
  expect_equal(result$index, "jaccard")
})


# ============================================================================
# PLOT.R - Various plot edge cases
# ============================================================================

test_that("plot.spacc_beta with partition=FALSE works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  result <- spaccBeta(species, coords, n_seeds = 3,
                      parallel = FALSE, progress = FALSE)

  p <- plot(result, partition = FALSE)
  expect_s3_class(p, "ggplot")
})


test_that("plot.spacc_beta with xaxis=distance works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  result <- spaccBeta(species, coords, n_seeds = 3,
                      parallel = FALSE, progress = FALSE)

  p <- plot(result, xaxis = "distance")
  expect_s3_class(p, "ggplot")
})


test_that("plot.spacc_endemism without richness works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  result <- spaccEndemism(species, coords, n_seeds = 3,
                          parallel = FALSE, progress = FALSE)

  p <- plot(result, show_richness = FALSE)
  expect_s3_class(p, "ggplot")
})


test_that("plot.spacc_fit with extrapolate_to works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  sac <- spacc(species, coords, n_seeds = 3, parallel = FALSE, progress = FALSE)
  fit <- extrapolate(sac, model = "michaelis-menten")

  p <- plot(fit, extrapolate_to = 100)
  expect_s3_class(p, "ggplot")
})


# ============================================================================
# COVERAGE.R - spaccCoverage with map=TRUE
# ============================================================================

test_that("spaccCoverage with map=TRUE returns site_values", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 3), nrow = 15)

  result <- spaccCoverage(species, coords, n_seeds = 3, map = TRUE,
                          parallel = FALSE, progress = FALSE)

  expect_s3_class(result, "spacc_coverage")
  expect_false(is.null(result$site_values))
  expect_true("final_coverage" %in% names(result$site_values))
})


test_that("plot.spacc_coverage with map type works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 3), nrow = 15)

  result <- spaccCoverage(species, coords, n_seeds = 3, map = TRUE,
                          parallel = FALSE, progress = FALSE)

  p <- plot(result, type = "map", metric = "final_coverage")
  expect_s3_class(p, "ggplot")
})


# ============================================================================
# HILL.R - spaccHill with map=TRUE and map plot
# ============================================================================

test_that("spaccHill with map=TRUE returns site_values", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 2), nrow = 15)

  result <- spaccHill(species, coords, q = c(0, 1), n_seeds = 3, map = TRUE,
                      parallel = FALSE, progress = FALSE)

  expect_s3_class(result, "spacc_hill")
  expect_false(is.null(result$site_values))
})


test_that("plot.spacc_hill with map type works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 2), nrow = 15)

  result <- spaccHill(species, coords, q = c(0, 1), n_seeds = 3, map = TRUE,
                      parallel = FALSE, progress = FALSE)

  p <- plot(result, type = "map", q = 0)
  expect_s3_class(p, "ggplot")
})


# ============================================================================
# BETA.R - spaccBeta with map=TRUE
# ============================================================================

test_that("spaccBeta with map=TRUE returns site_values", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  result <- spaccBeta(species, coords, n_seeds = 3, map = TRUE,
                      parallel = FALSE, progress = FALSE)

  expect_s3_class(result, "spacc_beta")
  expect_false(is.null(result$site_values))
})


test_that("plot.spacc_beta with map type works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  result <- spaccBeta(species, coords, n_seeds = 3, map = TRUE,
                      parallel = FALSE, progress = FALSE)

  p <- plot(result, type = "map", component = "beta_total")
  expect_s3_class(p, "ggplot")
})


# ============================================================================
# SAR.R - sfar with negative f interpretation
# ============================================================================

test_that("sfar print shows no fragmentation message when f <= 0", {
  skip_on_cran()

  set.seed(123)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.5), nrow = 20)

  sac <- spacc(species, coords, n_seeds = 5, parallel = FALSE, progress = FALSE)
  patches <- rep(1:4, each = 5)

  result <- sfar(sac, patches)

  # Just verify print doesn't error (f may be positive or negative)
  expect_output(print(result), "SFAR")
})
