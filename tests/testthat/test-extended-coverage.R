# Additional tests for extended coverage
# Targets: extrapolateCoverage, as.data.frame methods, as_sf methods, wavefront, dar, endemism

# ==============================================================================
# Test extrapolateCoverage with different q values
# ==============================================================================

test_that("extrapolateCoverage with q=0 works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rpois(30 * 15, 3), nrow = 30)

  cov <- spaccCoverage(species, coords, n_seeds = 3, parallel = FALSE,
                       progress = FALSE, seed = 1)
  ext <- extrapolateCoverage(cov, target_coverage = c(0.80, 0.95, 0.99), q = 0)

  expect_s3_class(ext, "spacc_coverage_ext")
  expect_equal(ext$q, 0)
  expect_equal(ncol(ext$richness), 3)
})


test_that("extrapolateCoverage with q=1 works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rpois(30 * 15, 3), nrow = 30)

  cov <- spaccCoverage(species, coords, n_seeds = 3, parallel = FALSE,
                       progress = FALSE, seed = 1)
  ext <- extrapolateCoverage(cov, target_coverage = c(0.80, 0.95, 0.99), q = 1)

  expect_s3_class(ext, "spacc_coverage_ext")
  expect_equal(ext$q, 1)
})


test_that("extrapolateCoverage with q=2 works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rpois(30 * 15, 3), nrow = 30)

  cov <- spaccCoverage(species, coords, n_seeds = 3, parallel = FALSE,
                       progress = FALSE, seed = 1)
  ext <- extrapolateCoverage(cov, target_coverage = c(0.80, 0.95, 0.99), q = 2)

  expect_s3_class(ext, "spacc_coverage_ext")
  expect_equal(ext$q, 2)
})


test_that("plot.spacc_coverage_ext works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  cov <- spaccCoverage(species, coords, n_seeds = 3, parallel = FALSE,
                       progress = FALSE, seed = 1)
  ext <- extrapolateCoverage(cov, target_coverage = c(0.90, 0.95), q = 0)

  expect_silent(plot(ext))
})


# ==============================================================================
# Test wavefront function and as.data.frame.spacc_wavefront
# ==============================================================================

test_that("wavefront works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(25), y = runif(25))
  species <- matrix(rbinom(25 * 10, 1, 0.3), nrow = 25)

  result <- wavefront(species, coords, n_seeds = 3,
                      progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc_wavefront")
  expect_true("radius" %in% names(result))
  expect_true("curves" %in% names(result))
})


test_that("as.data.frame.spacc_wavefront works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(25), y = runif(25))
  species <- matrix(rbinom(25 * 10, 1, 0.3), nrow = 25)

  result <- wavefront(species, coords, n_seeds = 3,
                      progress = FALSE, seed = 1)
  df <- as.data.frame(result)

  expect_s3_class(df, "data.frame")
  expect_true(all(c("step", "radius", "sites", "mean", "lower", "upper", "sd") %in% names(df)))
})


# ==============================================================================
# Test dar (diversity-area relationship)
# ==============================================================================

test_that("dar works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rpois(30 * 15, 2), nrow = 30)

  result <- dar(species, coords, q = c(0, 1), n_seeds = 3,
                parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc_dar")
  expect_equal(result$q, c(0, 1))
})


test_that("dar has expected structure", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rpois(30 * 15, 2), nrow = 30)

  result <- dar(species, coords, q = c(0, 1), n_seeds = 3,
                parallel = FALSE, progress = FALSE, seed = 1)

  expect_true("area" %in% names(result))
  expect_true("q" %in% names(result))
  expect_true("hill" %in% names(result))
})


# ==============================================================================
# Test spaccEndemism
# ==============================================================================

test_that("spaccEndemism works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(25), y = runif(25))
  species <- matrix(rbinom(25 * 12, 1, 0.3), nrow = 25)

  result <- spaccEndemism(species, coords, n_seeds = 3,
                          parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc_endemism")
})


test_that("spaccEndemism has expected structure", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(25), y = runif(25))
  species <- matrix(rbinom(25 * 12, 1, 0.3), nrow = 25)

  result <- spaccEndemism(species, coords, n_seeds = 3,
                          parallel = FALSE, progress = FALSE, seed = 1)

  expect_true("richness" %in% names(result))
  expect_true("endemism" %in% names(result))
})


# ==============================================================================
# Test spaccPhylo and as.data.frame.spacc_phylo
# ==============================================================================

test_that("spaccPhylo works", {
  skip_on_cran()
  skip_if_not_installed("ape")

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.4), nrow = 20)
  colnames(species) <- paste0("sp", 1:10)

  # Create simple tree
  tree <- ape::rtree(10, tip.label = paste0("sp", 1:10))

  result <- spaccPhylo(species, coords, tree = tree, metric = "pd",
                       n_seeds = 3, parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc_phylo")
})


test_that("spaccPhylo has expected structure", {
  skip_on_cran()
  skip_if_not_installed("ape")

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.4), nrow = 20)
  colnames(species) <- paste0("sp", 1:10)

  tree <- ape::rtree(10, tip.label = paste0("sp", 1:10))

  result <- spaccPhylo(species, coords, tree = tree, metric = "pd",
                       n_seeds = 3, parallel = FALSE, progress = FALSE, seed = 1)

  expect_true("curves" %in% names(result))
  expect_true("metric" %in% names(result))
  expect_equal(result$metric, "pd")
})


# ==============================================================================
# Test spaccFunc and as.data.frame.spacc_func
# ==============================================================================

test_that("spaccFunc works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 8, 1, 0.4), nrow = 20)
  colnames(species) <- paste0("sp", 1:8)

  # Create simple trait matrix
  traits <- data.frame(
    trait1 = rnorm(8),
    trait2 = rnorm(8),
    row.names = paste0("sp", 1:8)
  )

  result <- spaccFunc(species, coords, traits = traits, metric = "fdis",
                      n_seeds = 3, parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc_func")
})


test_that("spaccFunc has expected structure", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 8, 1, 0.4), nrow = 20)
  colnames(species) <- paste0("sp", 1:8)

  traits <- data.frame(
    trait1 = rnorm(8),
    trait2 = rnorm(8),
    row.names = paste0("sp", 1:8)
  )

  result <- spaccFunc(species, coords, traits = traits, metric = "fdis",
                      n_seeds = 3, parallel = FALSE, progress = FALSE, seed = 1)

  expect_true("curves" %in% names(result))
  expect_true("metric" %in% names(result))
  expect_equal(result$metric, "fdis")
})


# ==============================================================================
# Test spaccBetaFunc and spaccBetaPhylo
# ==============================================================================

test_that("spaccBetaFunc works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 8, 1, 0.4), nrow = 20)
  colnames(species) <- paste0("sp", 1:8)

  traits <- data.frame(
    trait1 = rnorm(8),
    trait2 = rnorm(8),
    row.names = paste0("sp", 1:8)
  )

  result <- spaccBetaFunc(species, coords, traits = traits,
                          n_seeds = 3, parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc_beta")
})


test_that("spaccBetaPhylo works", {
  skip_on_cran()
  skip_if_not_installed("ape")

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 8, 1, 0.4), nrow = 20)
  colnames(species) <- paste0("sp", 1:8)

  tree <- ape::rtree(8, tip.label = paste0("sp", 1:8))

  result <- spaccBetaPhylo(species, coords, tree = tree,
                           n_seeds = 3, parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc_beta")
})


# ==============================================================================
# Test distanceDecay
# ==============================================================================

test_that("distanceDecay works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(25), y = runif(25))
  species <- matrix(rbinom(25 * 12, 1, 0.3), nrow = 25)

  result <- distanceDecay(species, coords, n_seeds = 3, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc_decay")
})


test_that("as.data.frame.spacc_decay works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(25), y = runif(25))
  species <- matrix(rbinom(25 * 12, 1, 0.3), nrow = 25)

  result <- distanceDecay(species, coords, n_seeds = 3, progress = FALSE, seed = 1)
  df <- as.data.frame(result)

  expect_s3_class(df, "data.frame")
  expect_true(all(c("distance", "mean", "lower", "upper", "sd") %in% names(df)))
})


# ==============================================================================
# Test spaccCoverage with progress = TRUE (cli paths)
# ==============================================================================

test_that("spaccCoverage with progress shows messages", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 2), nrow = 15)

  expect_message(
    spaccCoverage(species, coords, n_seeds = 2, parallel = FALSE,
                  progress = TRUE, seed = 1),
    regexp = NULL
  )
})


# ==============================================================================
# Test spaccHill with progress
# ==============================================================================

test_that("spaccHill with progress shows messages", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 2), nrow = 15)

  expect_message(
    spaccHill(species, coords, q = 0, n_seeds = 2, parallel = FALSE,
              progress = TRUE, seed = 1),
    regexp = NULL
  )
})


# ==============================================================================
# Test spaccBeta with progress
# ==============================================================================

test_that("spaccBeta with progress shows messages", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  expect_message(
    spaccBeta(species, coords, n_seeds = 2, parallel = FALSE,
              progress = TRUE, seed = 1),
    regexp = NULL
  )
})


# ==============================================================================
# Test spaccMetrics with progress
# ==============================================================================

test_that("spaccMetrics with progress shows messages", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  expect_message(
    spaccMetrics(species, coords, parallel = FALSE, progress = TRUE),
    regexp = NULL
  )
})


# ==============================================================================
# Test print.spacc_dist
# ==============================================================================

test_that("print.spacc_dist works", {
  skip_on_cran()

  coords <- data.frame(x = c(0, 1, 2), y = c(0, 0, 0))
  d <- distances(coords)

  expect_output(print(d), "spacc")
})


# ==============================================================================
# Test as_sf methods
# ==============================================================================

test_that("as_sf.spacc_coverage works when sf is available", {
  skip_on_cran()
  skip_if_not_installed("sf")

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 2), nrow = 20)

  result <- spaccCoverage(species, coords, n_seeds = 20, map = TRUE,
                          parallel = FALSE, progress = FALSE, seed = 1)
  sf_obj <- as_sf(result, crs = 4326)

  expect_s3_class(sf_obj, "sf")
})


test_that("as_sf.spacc_hill works when sf is available", {
  skip_on_cran()
  skip_if_not_installed("sf")

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 2), nrow = 20)

  result <- spaccHill(species, coords, q = 0, n_seeds = 20, map = TRUE,
                      parallel = FALSE, progress = FALSE, seed = 1)
  sf_obj <- as_sf(result, crs = 4326)

  expect_s3_class(sf_obj, "sf")
})


test_that("as_sf.spacc_beta works when sf is available", {
  skip_on_cran()
  skip_if_not_installed("sf")

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  result <- spaccBeta(species, coords, n_seeds = 20, map = TRUE,
                      parallel = FALSE, progress = FALSE, seed = 1)
  sf_obj <- as_sf(result, crs = 4326)

  expect_s3_class(sf_obj, "sf")
})


# ==============================================================================
# Test subsample function
# ==============================================================================

test_that("subsample with random method works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(50), y = runif(50))

  result <- subsample(coords, n = 25, method = "random", seed = 1)

  expect_equal(length(result), 25)
})


test_that("subsample with grid method works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(50), y = runif(50))

  result <- subsample(coords, n = 25, method = "grid", seed = 1)

  expect_true(length(result) > 0)
})


# ==============================================================================
# Test spacc with radius method
# ==============================================================================

test_that("spacc with radius method works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(25), y = runif(25))
  species <- matrix(rbinom(25 * 10, 1, 0.3), nrow = 25)

  result <- spacc(species, coords, method = "radius", n_seeds = 3,
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
  expect_equal(result$method, "radius")
})


test_that("spacc with gaussian method works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(25), y = runif(25))
  species <- matrix(rbinom(25 * 10, 1, 0.3), nrow = 25)

  result <- spacc(species, coords, method = "gaussian", n_seeds = 3,
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
  expect_equal(result$method, "gaussian")
})


test_that("spacc with cone method works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(25), y = runif(25))
  species <- matrix(rbinom(25 * 10, 1, 0.3), nrow = 25)

  result <- spacc(species, coords, method = "cone", n_seeds = 3,
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
  expect_equal(result$method, "cone")
})


test_that("spacc with collector method works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(25), y = runif(25))
  species <- matrix(rbinom(25 * 10, 1, 0.3), nrow = 25)

  result <- spacc(species, coords, method = "collector", n_seeds = 3,
                  parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc")
  expect_equal(result$method, "collector")
})


# ==============================================================================
# Test plot.spacc with different types
# ==============================================================================

test_that("plot.spacc with type='raw' works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  result <- spacc(species, coords, n_seeds = 5, parallel = FALSE,
                  progress = FALSE, seed = 1)

  expect_silent(plot(result, type = "raw"))
})


# ==============================================================================
# Test spaccPhylo with different metrics
# ==============================================================================

test_that("spaccPhylo with mpd metric works", {
  skip_on_cran()
  skip_if_not_installed("ape")

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.4), nrow = 20)
  colnames(species) <- paste0("sp", 1:10)

  tree <- ape::rtree(10, tip.label = paste0("sp", 1:10))

  result <- spaccPhylo(species, coords, tree = tree, metric = "mpd",
                       n_seeds = 3, parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc_phylo")
  expect_equal(result$metric, "mpd")
})


test_that("spaccPhylo with mntd metric works", {
  skip_on_cran()
  skip_if_not_installed("ape")

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.4), nrow = 20)
  colnames(species) <- paste0("sp", 1:10)

  tree <- ape::rtree(10, tip.label = paste0("sp", 1:10))

  result <- spaccPhylo(species, coords, tree = tree, metric = "mntd",
                       n_seeds = 3, parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc_phylo")
  expect_equal(result$metric, "mntd")
})


# ==============================================================================
# Test spaccFunc with fric metric
# ==============================================================================

test_that("spaccFunc with fric metric works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 8, 1, 0.5), nrow = 20)
  colnames(species) <- paste0("sp", 1:8)

  traits <- data.frame(
    trait1 = rnorm(8),
    trait2 = rnorm(8),
    row.names = paste0("sp", 1:8)
  )

  result <- spaccFunc(species, coords, traits = traits, metric = "fric",
                      n_seeds = 3, parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc_func")
  expect_equal(result$metric, "fric")
})
