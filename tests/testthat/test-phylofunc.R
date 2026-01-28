test_that("MPD calculation is correct", {
  # Simple 4-species case with known distances
  dist_mat <- matrix(c(
    0, 1, 2, 3,
    1, 0, 1, 2,
    2, 1, 0, 1,
    3, 2, 1, 0
  ), nrow = 4)

  # All species present
  present <- c(TRUE, TRUE, TRUE, TRUE)

  # MPD = mean of all pairwise distances
  # Pairs: (1,2)=1, (1,3)=2, (1,4)=3, (2,3)=1, (2,4)=2, (3,4)=1
  # Mean = (1+2+3+1+2+1)/6 = 10/6 = 1.667
  expected_mpd <- 10/6

  result <- spacc:::calc_mpd(dist_mat, present)
  expect_equal(result, expected_mpd, tolerance = 1e-10)
})


test_that("MNTD calculation is correct", {
  dist_mat <- matrix(c(
    0, 1, 2, 3,
    1, 0, 1, 2,
    2, 1, 0, 1,
    3, 2, 1, 0
  ), nrow = 4)

  present <- c(TRUE, TRUE, TRUE, TRUE)

  # MNTD = mean of nearest neighbor distances
  # Species 1: nearest = species 2 (dist 1)
  # Species 2: nearest = species 1 or 3 (dist 1)
  # Species 3: nearest = species 2 or 4 (dist 1)
  # Species 4: nearest = species 3 (dist 1)
  # Mean = 1
  expected_mntd <- 1

  result <- spacc:::calc_mntd(dist_mat, present)
  expect_equal(result, expected_mntd, tolerance = 1e-10)
})


test_that("spaccPhylo returns correct structure", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  # Create simple distance matrix
  phylo_dist <- matrix(runif(8 * 8), nrow = 8)
  phylo_dist <- (phylo_dist + t(phylo_dist)) / 2  # Make symmetric
  diag(phylo_dist) <- 0

  result <- spaccPhylo(species, coords, phylo_dist,
                       metric = c("mpd", "mntd"), n_seeds = 3,
                       parallel = FALSE, progress = FALSE)

  expect_s3_class(result, "spacc_phylo")
  expect_equal(result$n_sites, 15)
  expect_equal(result$n_seeds, 3)
  expect_equal(result$metric, c("mpd", "mntd"))

  expect_equal(length(result$curves), 2)
})


test_that("spaccFunc returns correct structure", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 2), nrow = 15)

  traits <- matrix(rnorm(8 * 3), nrow = 8)

  result <- spaccFunc(species, coords, traits,
                      metric = c("fdis", "fric"), n_seeds = 3,
                      parallel = FALSE, progress = FALSE)

  expect_s3_class(result, "spacc_func")
  expect_equal(result$n_sites, 15)
  expect_equal(result$n_seeds, 3)
  expect_equal(result$n_traits, 3)
  expect_equal(result$metric, c("fdis", "fric"))

  expect_equal(length(result$curves), 2)
})


test_that("FDis increases with more species",
{
  skip_on_cran()

  set.seed(123)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 3), nrow = 20)
  traits <- matrix(rnorm(10 * 2), nrow = 10)

  result <- spaccFunc(species, coords, traits,
                      metric = "fdis", n_seeds = 3,
                      parallel = FALSE, progress = FALSE)

  # FDis should generally increase or stay stable as more species added
  # Just check final > initial for most seeds
  increases <- 0
  for (seed in 1:3) {
    curve <- result$curves$fdis[seed, ]
    # Find first non-zero value
    first_nonzero <- which(curve > 0)[1]
    if (!is.na(first_nonzero) && curve[length(curve)] >= curve[first_nonzero]) {
      increases <- increases + 1
    }
  }
  expect_true(increases >= 2)
})


test_that("spaccPhylo print works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  phylo_dist <- matrix(runif(8 * 8), nrow = 8)
  phylo_dist <- (phylo_dist + t(phylo_dist)) / 2
  diag(phylo_dist) <- 0

  result <- spaccPhylo(species, coords, phylo_dist,
                       metric = c("mpd", "mntd"), n_seeds = 3,
                       parallel = FALSE, progress = FALSE)

  expect_output(print(result), "spacc phylogenetic diversity")
})


test_that("spaccPhylo summary returns data.frame", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  phylo_dist <- matrix(runif(8 * 8), nrow = 8)
  phylo_dist <- (phylo_dist + t(phylo_dist)) / 2
  diag(phylo_dist) <- 0

  result <- spaccPhylo(species, coords, phylo_dist,
                       metric = c("mpd", "mntd"), n_seeds = 3,
                       parallel = FALSE, progress = FALSE)

  summ <- summary(result)
  expect_s3_class(summ, "data.frame")
  expect_true("metric" %in% names(summ))
  expect_true("mean" %in% names(summ))
})


test_that("spaccPhylo with spacc_dist coords", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  phylo_dist <- matrix(runif(8 * 8), nrow = 8)
  phylo_dist <- (phylo_dist + t(phylo_dist)) / 2
  diag(phylo_dist) <- 0

  d <- distances(coords)
  result <- spaccPhylo(species, d, phylo_dist,
                       metric = "mpd", n_seeds = 3,
                       parallel = FALSE, progress = FALSE)

  expect_s3_class(result, "spacc_phylo")
})


test_that("spaccPhylo errors on bad tree", {
  coords <- data.frame(x = runif(10), y = runif(10))
  species <- matrix(rbinom(10 * 5, 1, 0.3), nrow = 10)

  expect_error(
    spaccPhylo(species, coords, "not_a_tree",
               metric = "mpd", n_seeds = 1,
               parallel = FALSE, progress = FALSE),
    "tree must be"
  )
})


test_that("spaccPhylo PD warning for distance matrix", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  phylo_dist <- matrix(runif(8 * 8), nrow = 8)
  phylo_dist <- (phylo_dist + t(phylo_dist)) / 2
  diag(phylo_dist) <- 0

  expect_warning(
    spaccPhylo(species, coords, phylo_dist,
               metric = c("mpd", "pd"), n_seeds = 3,
               parallel = FALSE, progress = FALSE),
    "PD requires full tree"
  )
})


test_that("spaccFunc print works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 2), nrow = 15)
  traits <- matrix(rnorm(8 * 3), nrow = 8)

  result <- spaccFunc(species, coords, traits,
                      metric = c("fdis", "fric"), n_seeds = 3,
                      parallel = FALSE, progress = FALSE)

  expect_output(print(result), "spacc functional diversity")
})


test_that("spaccFunc summary returns data.frame", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 2), nrow = 15)
  traits <- matrix(rnorm(8 * 3), nrow = 8)

  result <- spaccFunc(species, coords, traits,
                      metric = c("fdis", "fric"), n_seeds = 3,
                      parallel = FALSE, progress = FALSE)

  summ <- summary(result)
  expect_s3_class(summ, "data.frame")
  expect_true("metric" %in% names(summ))
})


test_that("spaccFunc errors on species not in traits", {
  coords <- data.frame(x = runif(10), y = runif(10))
  species <- matrix(rpois(10 * 5, 2), nrow = 10)
  colnames(species) <- paste0("sp", 1:5)
  traits <- matrix(rnorm(3 * 2), nrow = 3)
  rownames(traits) <- paste0("sp", 1:3)

  expect_error(
    spaccFunc(species, coords, traits, metric = "fdis", n_seeds = 1,
              parallel = FALSE, progress = FALSE),
    "Some species"
  )
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
                      parallel = FALSE, progress = FALSE)

  expect_s3_class(result, "spacc_func")
})


test_that("spaccPhylo as_sf errors without map data", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  phylo_dist <- matrix(runif(8 * 8), nrow = 8)
  phylo_dist <- (phylo_dist + t(phylo_dist)) / 2
  diag(phylo_dist) <- 0

  result <- spaccPhylo(species, coords, phylo_dist,
                       metric = "mpd", n_seeds = 3,
                       parallel = FALSE, progress = FALSE)

  expect_error(as_sf(result), "No map data")
})


test_that("spaccFunc as_sf errors without map data", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 2), nrow = 15)
  traits <- matrix(rnorm(8 * 3), nrow = 8)

  result <- spaccFunc(species, coords, traits,
                      metric = "fdis", n_seeds = 3,
                      parallel = FALSE, progress = FALSE)

  expect_error(as_sf(result), "No map data")
})


test_that("spaccPhylo with map=TRUE works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  phylo_dist <- matrix(runif(8 * 8), nrow = 8)
  phylo_dist <- (phylo_dist + t(phylo_dist)) / 2
  diag(phylo_dist) <- 0

  result <- spaccPhylo(species, coords, phylo_dist,
                       metric = "mpd", n_seeds = 3, map = TRUE,
                       parallel = FALSE, progress = FALSE)

  expect_s3_class(result, "spacc_phylo")
  expect_true(!is.null(result$site_values))
})


test_that("spaccFunc with map=TRUE works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 2), nrow = 15)
  traits <- matrix(rnorm(8 * 3), nrow = 8)

  result <- spaccFunc(species, coords, traits,
                      metric = "fdis", n_seeds = 3, map = TRUE,
                      parallel = FALSE, progress = FALSE)

  expect_s3_class(result, "spacc_func")
  expect_true(!is.null(result$site_values))
})


test_that("plot.spacc_phylo map errors without map data", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  phylo_dist <- matrix(runif(8 * 8), nrow = 8)
  phylo_dist <- (phylo_dist + t(phylo_dist)) / 2
  diag(phylo_dist) <- 0

  result <- spaccPhylo(species, coords, phylo_dist,
                       metric = "mpd", n_seeds = 3,
                       parallel = FALSE, progress = FALSE)

  expect_error(plot(result, type = "map"), "map")
})


test_that("plot.spacc_func map errors without map data", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 2), nrow = 15)
  traits <- matrix(rnorm(8 * 3), nrow = 8)

  result <- spaccFunc(species, coords, traits,
                      metric = "fdis", n_seeds = 3,
                      parallel = FALSE, progress = FALSE)

  expect_error(plot(result, type = "map"), "map")
})
