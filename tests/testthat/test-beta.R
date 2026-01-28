test_that("spaccBeta returns correct structure", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  result <- spaccBeta(species, coords, n_seeds = 5,
                      parallel = FALSE, progress = FALSE)

  expect_s3_class(result, "spacc_beta")
  expect_equal(result$n_sites, 20)
  expect_equal(result$n_seeds, 5)
  expect_equal(result$index, "sorensen")

  # Beta matrices should be n_seeds x (n_sites - 1)
  expect_equal(dim(result$beta_total), c(5, 19))
  expect_equal(dim(result$beta_turnover), c(5, 19))
  expect_equal(dim(result$beta_nestedness), c(5, 19))
})


test_that("Beta components sum to total", {
  skip_on_cran()

  set.seed(123)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  result <- spaccBeta(species, coords, n_seeds = 3,
                      parallel = FALSE, progress = FALSE)

  # turnover + nestedness should equal total (within tolerance)
  computed_total <- result$beta_turnover + result$beta_nestedness
  expect_equal(computed_total, result$beta_total, tolerance = 1e-10)
})


test_that("Beta diversity values are bounded 0-1", {
  skip_on_cran()

  set.seed(456)
  coords <- data.frame(x = runif(25), y = runif(25))
  species <- matrix(rbinom(25 * 12, 1, 0.35), nrow = 25)

  result <- spaccBeta(species, coords, n_seeds = 5,
                      parallel = FALSE, progress = FALSE)

  expect_true(all(result$beta_total >= 0 & result$beta_total <= 1))
  expect_true(all(result$beta_turnover >= 0 & result$beta_turnover <= 1))
  expect_true(all(result$beta_nestedness >= 0 & result$beta_nestedness <= 1))
})


test_that("Jaccard and Sorensen give different results", {
  skip_on_cran()

  set.seed(789)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  result_sor <- spaccBeta(species, coords, n_seeds = 3, index = "sorensen",
                          parallel = FALSE, progress = FALSE, seed = 1)
  result_jac <- spaccBeta(species, coords, n_seeds = 3, index = "jaccard",
                          parallel = FALSE, progress = FALSE, seed = 1)

  # Should be different
  expect_false(identical(result_sor$beta_total, result_jac$beta_total))
})


test_that("spaccBeta print works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  result <- spaccBeta(species, coords, n_seeds = 3,
                      parallel = FALSE, progress = FALSE, seed = 1)

  expect_output(print(result), "spacc beta diversity")
})


test_that("spaccBeta summary returns data.frame", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.3), nrow = 20)

  result <- spaccBeta(species, coords, n_seeds = 3,
                      parallel = FALSE, progress = FALSE, seed = 1)

  summ <- summary(result)
  expect_s3_class(summ, "data.frame")
  expect_true("step" %in% names(summ))
  expect_true("beta_total" %in% names(summ))
  expect_true("beta_turnover" %in% names(summ))
  expect_true("beta_nestedness" %in% names(summ))
  expect_true("mean_distance" %in% names(summ))
  expect_equal(nrow(summ), 19)
})


test_that("spaccBeta with map = TRUE stores site_values", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  result <- spaccBeta(species, coords, n_seeds = 3, map = TRUE,
                      parallel = FALSE, progress = FALSE, seed = 1)

  expect_false(is.null(result$site_values))
  expect_equal(nrow(result$site_values), 15)
  expect_true("beta_total" %in% names(result$site_values))
})


test_that("spaccBeta map error without map data", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  result <- spaccBeta(species, coords, n_seeds = 3,
                      parallel = FALSE, progress = FALSE, seed = 1)

  expect_error(as_sf(result), "No map data")
})


test_that("spaccBeta with spacc_dist coords", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  d <- distances(coords)
  result <- spaccBeta(species, d, n_seeds = 3,
                      parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc_beta")
})


test_that("spaccBetaFunc returns correct structure", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)
  traits <- matrix(rnorm(8 * 3), nrow = 8)
  colnames(species) <- paste0("sp", 1:8)
  rownames(traits) <- paste0("sp", 1:8)

  result <- spaccBetaFunc(species, coords, traits, n_seeds = 3,
                           parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc_beta")
  expect_equal(result$beta_type, "functional")
  expect_equal(result$n_seeds, 3)
})


test_that("spaccBetaFunc jaccard works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)
  traits <- matrix(rnorm(8 * 3), nrow = 8)

  result <- spaccBetaFunc(species, coords, traits, n_seeds = 3,
                           index = "jaccard",
                           parallel = FALSE, progress = FALSE, seed = 1)

  expect_equal(result$index, "jaccard")
})


test_that("spaccBetaPhylo returns correct structure", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  phylo_dist <- matrix(runif(8 * 8), nrow = 8)
  phylo_dist <- (phylo_dist + t(phylo_dist)) / 2
  diag(phylo_dist) <- 0

  result <- spaccBetaPhylo(species, coords, phylo_dist, n_seeds = 3,
                            parallel = FALSE, progress = FALSE, seed = 1)

  expect_s3_class(result, "spacc_beta")
  expect_equal(result$beta_type, "phylogenetic")
})


test_that("spaccBetaPhylo jaccard works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)

  phylo_dist <- matrix(runif(8 * 8), nrow = 8)
  phylo_dist <- (phylo_dist + t(phylo_dist)) / 2
  diag(phylo_dist) <- 0

  result <- spaccBetaPhylo(species, coords, phylo_dist, n_seeds = 3,
                            index = "jaccard",
                            parallel = FALSE, progress = FALSE, seed = 1)

  expect_equal(result$index, "jaccard")
})


test_that("spaccBetaPhylo errors on bad tree input", {
  coords <- data.frame(x = runif(10), y = runif(10))
  species <- matrix(rbinom(10 * 5, 1, 0.3), nrow = 10)

  expect_error(
    spaccBetaPhylo(species, coords, "not_a_tree", n_seeds = 1,
                    parallel = FALSE, progress = FALSE),
    "tree must be"
  )
})


test_that("spaccBetaFunc errors on missing species in traits", {
  coords <- data.frame(x = runif(10), y = runif(10))
  species <- matrix(rbinom(10 * 5, 1, 0.3), nrow = 10)
  colnames(species) <- paste0("sp", 1:5)
  traits <- matrix(rnorm(3 * 3), nrow = 3)
  rownames(traits) <- paste0("sp", 1:3)

  expect_error(
    spaccBetaFunc(species, coords, traits, n_seeds = 1,
                   parallel = FALSE, progress = FALSE),
    "Some species"
  )
})
