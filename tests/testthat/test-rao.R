test_that("calc_rao matches manual Rao's Q", {
  d <- matrix(c(0, 1, 2,
                1, 0, 3,
                2, 3, 0), 3, 3, byrow = TRUE)
  ab <- c(2, 3, 5)
  p <- ab / sum(ab)
  manual <- sum(outer(p, p) * d)          # sum_i sum_j p_i p_j d_ij
  expect_equal(spacc:::calc_rao(d, ab), manual)
})


test_that("calc_rao is zero for one species or zero distances", {
  d <- matrix(c(0, 1, 1, 0), 2, 2)
  expect_equal(spacc:::calc_rao(d, c(5, 0)), 0)               # only one present
  expect_equal(spacc:::calc_rao(matrix(0, 2, 2), c(2, 3)), 0) # all distances 0
})


test_that("calc_rao_traits matches manual trait-based Rao", {
  traits <- matrix(c(0, 0,
                     1, 0,
                     0, 2), 3, 2, byrow = TRUE)
  ab <- c(1, 2, 4)
  p <- ab / sum(ab)
  d <- as.matrix(stats::dist(traits))     # Euclidean trait distance
  manual <- sum(outer(p, p) * d)
  expect_equal(spacc:::calc_rao_traits(traits, ab), manual)
})


test_that("spaccPhylo rao final step matches calc_rao on the full community", {
  skip_if_not_installed("ape")
  set.seed(11)
  tree <- ape::rtree(8)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 2), nrow = 15)
  colnames(species) <- tree$tip.label

  res <- spaccPhylo(species, coords, tree, metric = "rao",
                    n_seeds = 4, progress = FALSE)
  expect_true("rao" %in% res$metric)
  rao_mat <- res$curves[[which(res$metric == "rao")]]
  expect_equal(dim(rao_mat), c(4L, 15L))

  phylo_d <- as.matrix(ape::cophenetic.phylo(tree))[colnames(species), colnames(species)]
  full <- colSums(species)
  expect_equal(round(unique(rao_mat[, 15]), 8),
               round(spacc:::calc_rao(phylo_d, full), 8))
})


test_that("spaccPhylo rao does not alter mpd/mntd", {
  skip_if_not_installed("ape")
  set.seed(13)
  tree <- ape::rtree(8)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)
  colnames(species) <- tree$tip.label

  a <- spaccPhylo(species, coords, tree, metric = c("mpd", "mntd"),
                  n_seeds = 3, seed = 7, progress = FALSE)
  b <- spaccPhylo(species, coords, tree, metric = c("mpd", "mntd", "rao"),
                  n_seeds = 3, seed = 7, progress = FALSE)
  expect_equal(a$curves[["mpd"]], b$curves[["mpd"]])
  expect_equal(a$curves[["mntd"]], b$curves[["mntd"]])
})


test_that("spaccFunc rao final step matches calc_rao_traits", {
  set.seed(12)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 6, 2), nrow = 15)
  traits <- matrix(rnorm(6 * 3), nrow = 6)
  rownames(traits) <- paste0("sp", 1:6)
  colnames(species) <- rownames(traits)

  res <- spaccFunc(species, coords, traits, metric = "rao",
                   n_seeds = 4, progress = FALSE)
  rao_mat <- res$curves[[which(res$metric == "rao")]]
  full <- colSums(species)
  expect_equal(round(unique(rao_mat[, 15]), 8),
               round(spacc:::calc_rao_traits(traits, full), 8))
})
