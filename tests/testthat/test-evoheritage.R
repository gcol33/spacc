evo_fixture <- function(n_sites = 30, n_species = 20, seed = 1) {
  set.seed(seed)
  tree <- ape::rtree(n_species)
  species <- matrix(rbinom(n_sites * n_species, 1, 0.3), nrow = n_sites)
  colnames(species) <- tree$tip.label
  coords <- data.frame(x = runif(n_sites), y = runif(n_sites))
  list(tree = tree, species = species, coords = coords)
}

test_that("calc_evo_heritage recovers PD at rho = 0 and richness at large rho", {
  skip_if_not_installed("ape")
  f <- evo_fixture()
  present <- seq_along(f$tree$tip.label) %in% c(2, 5, 7, 11)
  pd <- spacc:::calc_faith_pd(f$tree$edge, f$tree$edge.length,
                              length(f$tree$tip.label), present)
  expect_equal(spacc:::calc_evo_heritage(f$tree$edge, f$tree$edge.length,
                                         length(f$tree$tip.label), present, 0), pd)
  expect_equal(spacc:::calc_evo_heritage(f$tree$edge, f$tree$edge.length,
                                         length(f$tree$tip.label), present, 1e4),
               4, tolerance = 1e-6)
})

test_that("evoheritage curve at rho = 0 equals the pd curve", {
  skip_if_not_installed("ape")
  f <- evo_fixture()
  res <- spaccPhylo(f$species, f$coords, f$tree, metric = c("pd", "evoheritage"),
                    rho = 0, n_seeds = 5, parallel = FALSE, progress = FALSE, seed = 3)
  expect_equal(res$curves$evoheritage, res$curves$pd)
})

test_that("vector rho gives one curve per value, each equal to its own run", {
  skip_if_not_installed("ape")
  f <- evo_fixture()
  rho <- c(0.1, 1, 10)
  res <- spaccPhylo(f$species, f$coords, f$tree, metric = c("pd", "evoheritage"),
                    rho = rho, n_seeds = 4, parallel = FALSE, progress = FALSE, seed = 3)
  labels <- paste0("evoheritage_rho", rho)
  expect_equal(res$metric, c("pd", labels))
  expect_equal(res$rho, rho)
  for (k in seq_along(rho)) {
    single <- spaccPhylo(f$species, f$coords, f$tree, metric = "evoheritage",
                         rho = rho[k], n_seeds = 4, parallel = FALSE, progress = FALSE,
                         seed = 3)
    expect_equal(res$curves[[labels[k]]], single$curves$evoheritage)
  }
  expect_s3_class(summary(res), "data.frame")
  expect_s3_class(plot(res), "ggplot")
})

test_that("evoheritage curve is non-decreasing along the accumulation", {
  skip_if_not_installed("ape")
  f <- evo_fixture()
  res <- spaccPhylo(f$species, f$coords, f$tree, metric = "evoheritage", rho = 2,
                    n_seeds = 5, parallel = FALSE, progress = FALSE, seed = 4)
  steps <- t(apply(res$curves$evoheritage, 1, diff))
  expect_true(all(steps >= -1e-12))
})

test_that("tree metrics match species to tips by name, not by column position", {
  skip_if_not_installed("ape")
  f <- evo_fixture()
  perm <- rev(seq_len(ncol(f$species)))
  a <- spaccPhylo(f$species, f$coords, f$tree, metric = c("pd", "evoheritage"),
                  rho = 0.5, n_seeds = 3, parallel = FALSE, progress = FALSE, seed = 9)
  b <- spaccPhylo(f$species[, perm], f$coords, f$tree, metric = c("pd", "evoheritage"),
                  rho = 0.5, n_seeds = 3, parallel = FALSE, progress = FALSE, seed = 9)
  expect_equal(a$curves$pd, b$curves$pd)
  expect_equal(a$curves$evoheritage, b$curves$evoheritage)
})

test_that("PD ignores tree tips absent from the community", {
  skip_if_not_installed("ape")
  f <- evo_fixture()
  keep <- 1:12
  sub <- f$species[, keep]
  res <- spaccPhylo(sub, f$coords, f$tree, metric = "pd",
                    n_seeds = 2, parallel = FALSE, progress = FALSE, seed = 2)
  all_present <- f$tree$tip.label %in% colnames(sub)[colSums(sub) > 0]
  expect_equal(res$curves$pd[1, nrow(sub)],
               spacc:::calc_faith_pd(f$tree$edge, f$tree$edge.length,
                                     length(f$tree$tip.label), all_present))
})

test_that("evoheritage is dropped with a warning for a distance matrix", {
  skip_if_not_installed("ape")
  f <- evo_fixture()
  d <- ape::cophenetic.phylo(f$tree)
  expect_warning(
    res <- spaccPhylo(f$species, f$coords, d, metric = c("mpd", "evoheritage"),
                      n_seeds = 2, parallel = FALSE, progress = FALSE),
    "requires a full tree")
  expect_equal(res$metric, "mpd")
})
