red_fixture <- function(n_sp = 12, n_sites = 20, seed = 6) {
  set.seed(seed)
  tree <- ape::rcoal(n_sp)
  species <- matrix(rpois(n_sites * n_sp, 1.2), nrow = n_sites,
                    dimnames = list(NULL, tree$tip.label))
  list(tree = tree, species = species,
       coords = data.frame(x = runif(n_sites), y = runif(n_sites)))
}

test_that("phylogenetic Hill numbers match iNEXT.3D mean PD", {
  skip_if_not_installed("ape")
  skip_if_not_installed("iNEXT.3D")
  f <- red_fixture()
  ab <- colSums(f$species)
  ours <- spacc:::calc_hill_phylo(f$tree$edge, f$tree$edge.length,
                                  length(f$tree$tip.label), ab, c(0, 1, 2))
  ref <- iNEXT.3D::ObsAsy3D(ab, diversity = "PD", q = c(0, 1, 2), datatype = "abundance",
                            method = "Observed", PDtree = f$tree, PDtype = "meanPD",
                            nboot = 0)
  expect_equal(ours, ref$qPD[order(ref$Order.q)], tolerance = 1e-8)
})

test_that("similarity Hill numbers follow the Leinster-Cobbold formula", {
  set.seed(2)
  ab <- c(5, 0, 3, 1, 8)
  Z <- matrix(runif(25, 0, 0.6), 5)
  Z <- (Z + t(Z)) / 2
  diag(Z) <- 1
  p <- ab / sum(ab)
  keep <- p > 0
  zp <- as.vector(Z[keep, keep] %*% p[keep])
  ref <- c((sum(p[keep] * zp^(0 - 1)))^(1 / (1 - 0)),
           exp(-sum(p[keep] * log(zp))),
           (sum(p[keep] * zp^(2 - 1)))^(1 / (1 - 2)))
  expect_equal(spacc:::calc_hill_similarity(Z, ab, c(0, 1, 2)), ref)
  expect_equal(spacc:::calc_hill_similarity(diag(5), ab, c(0, 1, 2)),
               vapply(c(0, 1, 2), function(q) spacc:::calc_hill_number(ab, q), numeric(1)))
})

test_that("redundancy is 1 on a star tree at every step and order", {
  skip_if_not_installed("ape")
  f <- red_fixture()
  star <- ape::stree(12, "star", tip.label = f$tree$tip.label)
  star$edge.length <- rep(1, nrow(star$edge))
  red <- spaccRedundancy(f$species, f$coords, tree = star, n_seeds = 3,
                         parallel = FALSE, progress = FALSE, seed = 1)
  for (m in red$curves) expect_equal(m, matrix(1, 3, 20), ignore_attr = TRUE)
})

test_that("redundancy is 1 for equidistant traits", {
  f <- red_fixture()
  traits <- diag(12)
  rownames(traits) <- colnames(f$species)
  red <- spaccRedundancy(f$species, f$coords, traits = traits, n_seeds = 3,
                         parallel = FALSE, progress = FALSE, seed = 1)
  for (m in red$curves) expect_equal(m, matrix(1, 3, 20), ignore_attr = TRUE, tolerance = 1e-10)
})

test_that("redundancy is at least 1 and its parts match spaccHill and the kernels", {
  skip_if_not_installed("ape")
  f <- red_fixture()
  red <- spaccRedundancy(f$species, f$coords, tree = f$tree, n_seeds = 4,
                         parallel = FALSE, progress = FALSE, seed = 2)
  for (m in red$curves) expect_true(all(m >= 1 - 1e-10))
  hill <- spaccHill(f$species, f$coords, q = c(0, 1, 2), n_seeds = 4,
                    parallel = FALSE, progress = FALSE, seed = 2)
  expect_equal(unname(red$taxonomic), unname(hill$curves), ignore_attr = TRUE)
  final <- spacc:::calc_hill_phylo(f$tree$edge, f$tree$edge.length, 12,
                                   colSums(f$species), c(0, 1, 2))
  expect_equal(vapply(red$dimension, function(m) m[1, 20], numeric(1)), final,
               ignore_attr = TRUE)
  expect_s3_class(summary(red), "data.frame")
  expect_s3_class(plot(red), "ggplot")
  expect_output(print(red), "phylogenetic redundancy")
})

test_that("spaccRedundancy needs exactly one of tree and traits", {
  f <- red_fixture()
  expect_error(spaccRedundancy(f$species, f$coords, progress = FALSE), "either")
})
