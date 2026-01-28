test_that("plot.spacc returns ggplot object", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  sac <- spacc(species, coords, n_seeds = 3, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)

  p <- plot(sac)
  expect_s3_class(p, "ggplot")
})


test_that("plot.spacc with show_seeds works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  sac <- spacc(species, coords, n_seeds = 3, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)

  p <- plot(sac, show_seeds = TRUE)
  expect_s3_class(p, "ggplot")
})


test_that("plot.spacc with saturation marker works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  sac <- spacc(species, coords, n_seeds = 3, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)

  p <- plot(sac, saturation = TRUE)
  expect_s3_class(p, "ggplot")
})


test_that("plot.spacc with no CI works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  sac <- spacc(species, coords, n_seeds = 3, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)

  p <- plot(sac, ci = FALSE)
  expect_s3_class(p, "ggplot")
})


test_that("plot.spacc for grouped object works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  sp1 <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)
  sp2 <- matrix(rbinom(15 * 8, 1, 0.2), nrow = 15)

  sac1 <- spacc(sp1, coords, n_seeds = 3, method = "knn",
                parallel = FALSE, progress = FALSE, seed = 1)
  sac2 <- spacc(sp2, coords, n_seeds = 3, method = "knn",
                parallel = FALSE, progress = FALSE, seed = 1)

  combined <- c(native = sac1, alien = sac2)
  p <- plot(combined)
  expect_s3_class(p, "ggplot")
})


test_that("plot.spacc grouped with facet works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  sp1 <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)
  sp2 <- matrix(rbinom(15 * 8, 1, 0.2), nrow = 15)

  sac1 <- spacc(sp1, coords, n_seeds = 3, method = "knn",
                parallel = FALSE, progress = FALSE, seed = 1)
  sac2 <- spacc(sp2, coords, n_seeds = 3, method = "knn",
                parallel = FALSE, progress = FALSE, seed = 1)

  combined <- c(native = sac1, alien = sac2)
  p <- plot(combined, facet = TRUE)
  expect_s3_class(p, "ggplot")
})


test_that("plot.spacc_fit returns ggplot", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)
  fit <- extrapolate(sac, model = "michaelis-menten")

  p <- plot(fit)
  expect_s3_class(p, "ggplot")
})


test_that("plot.spacc_comp returns ggplot", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  sp1 <- matrix(rbinom(20 * 10, 1, 0.4), nrow = 20)
  sp2 <- matrix(rbinom(20 * 10, 1, 0.2), nrow = 20)

  sac1 <- spacc(sp1, coords, n_seeds = 3, method = "knn",
                parallel = FALSE, progress = FALSE, seed = 1)
  sac2 <- spacc(sp2, coords, n_seeds = 3, method = "knn",
                parallel = FALSE, progress = FALSE, seed = 1)

  comp <- compare(sac1, sac2, n_perm = 19)
  p <- plot(comp)
  expect_s3_class(p, "ggplot")
})


test_that("plot.spacc_rare returns ggplot", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  rare <- rarefy(species, n_boot = 10)
  p <- plot(rare)
  expect_s3_class(p, "ggplot")
})


test_that("plot.spacc_rare without CI works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  rare <- rarefy(species, n_boot = 10)
  p <- plot(rare, ci = FALSE)
  expect_s3_class(p, "ggplot")
})


test_that("plot.spacc_dist histogram works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  coords <- data.frame(x = runif(10), y = runif(10))
  d <- distances(coords)

  p <- plot(d, type = "histogram")
  expect_s3_class(p, "ggplot")
})


test_that("plot.spacc_dist heatmap works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  coords <- data.frame(x = runif(10), y = runif(10))
  d <- distances(coords)

  p <- plot(d, type = "heatmap")
  expect_s3_class(p, "ggplot")
})


test_that("plot.spacc_metrics heatmap works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  result <- spaccMetrics(species, coords,
                         metrics = c("slope_10", "auc"),
                         parallel = FALSE, progress = FALSE)

  p <- plot(result, metric = "slope_10", type = "heatmap")
  expect_s3_class(p, "ggplot")
})


test_that("plot.spacc_metrics histogram works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  result <- spaccMetrics(species, coords,
                         metrics = c("slope_10", "auc"),
                         parallel = FALSE, progress = FALSE)

  p <- plot(result, metric = "auc", type = "histogram")
  expect_s3_class(p, "ggplot")
})


test_that("plot.spacc_metrics errors on unknown metric", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  result <- spaccMetrics(species, coords,
                         metrics = c("slope_10"),
                         parallel = FALSE, progress = FALSE)

  expect_error(plot(result, metric = "nonexistent"), "not found")
})


test_that("plot.spacc_beta returns ggplot", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  result <- spaccBeta(species, coords, n_seeds = 3,
                      parallel = FALSE, progress = FALSE, seed = 1)

  p <- plot(result)
  expect_s3_class(p, "ggplot")
})


test_that("plot.spacc_hill returns ggplot", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 2), nrow = 20)

  result <- spaccHill(species, coords, q = c(0, 1, 2), n_seeds = 3,
                      parallel = FALSE, progress = FALSE)

  p <- plot(result)
  expect_s3_class(p, "ggplot")
})


test_that("plot.spacc_coverage returns ggplot", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- spaccCoverage(species, coords, n_seeds = 3,
                          parallel = FALSE, progress = FALSE, seed = 1)

  p <- plot(result)
  expect_s3_class(p, "ggplot")
})


test_that("plot.spacc_phylo returns ggplot", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  phylo_dist <- matrix(runif(8 * 8), nrow = 8)
  phylo_dist <- (phylo_dist + t(phylo_dist)) / 2
  diag(phylo_dist) <- 0

  result <- spaccPhylo(species, coords, phylo_dist,
                       metric = c("mpd", "mntd"), n_seeds = 3,
                       parallel = FALSE, progress = FALSE)

  p <- plot(result)
  expect_s3_class(p, "ggplot")
})


test_that("plot.spacc_func returns ggplot", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 2), nrow = 15)
  traits <- matrix(rnorm(8 * 3), nrow = 8)

  result <- spaccFunc(species, coords, traits,
                      metric = c("fdis", "fric"), n_seeds = 3,
                      parallel = FALSE, progress = FALSE)

  p <- plot(result)
  expect_s3_class(p, "ggplot")
})


test_that("plot.spacc_dar returns ggplot", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 2), nrow = 20)

  result <- dar(species, coords, q = c(0, 1), n_seeds = 3,
                area_method = "count",
                parallel = FALSE, progress = FALSE, seed = 1)

  p <- plot(result)
  expect_s3_class(p, "ggplot")
})


test_that("plot.spacc_endemism returns ggplot", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  result <- spaccEndemism(species, coords, n_seeds = 3,
                           parallel = FALSE, progress = FALSE, seed = 1)

  p <- plot(result)
  expect_s3_class(p, "ggplot")
})
