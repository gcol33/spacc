# Deprecated function names forward to their current front doors.
# spaccBetaFunc/spaccBetaPhylo -> spaccBeta(traits=/tree=)
# diversityProfileFunc/diversityProfilePhylo -> diversityProfile(traits=/tree=)
# wavefront -> spaccWavefront

test_that("spaccBetaFunc is deprecated and forwards to spaccBeta(traits=)", {
  skip_on_cran()
  set.seed(1)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)
  colnames(species) <- paste0("sp", 1:8)
  traits <- matrix(rnorm(8 * 3), nrow = 8, dimnames = list(paste0("sp", 1:8), NULL))

  expect_warning(
    spaccBetaFunc(species, coords, traits, n_seeds = 3,
                  parallel = FALSE, progress = FALSE, seed = 1),
    "deprecated"
  )

  set.seed(7)
  new <- spaccBeta(species, coords, traits = traits, n_seeds = 3,
                   parallel = FALSE, progress = FALSE)
  set.seed(7)
  old <- suppressWarnings(
    spaccBetaFunc(species, coords, traits, n_seeds = 3,
                  parallel = FALSE, progress = FALSE)
  )
  expect_equal(old$beta_total, new$beta_total)
  expect_equal(old$beta_type, "functional")
})


test_that("spaccBetaPhylo is deprecated and forwards to spaccBeta(tree=)", {
  skip_on_cran()
  set.seed(1)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.3), nrow = 15)
  phylo_dist <- matrix(runif(8 * 8), nrow = 8)
  phylo_dist <- (phylo_dist + t(phylo_dist)) / 2
  diag(phylo_dist) <- 0

  expect_warning(
    spaccBetaPhylo(species, coords, phylo_dist, n_seeds = 3,
                   parallel = FALSE, progress = FALSE, seed = 1),
    "deprecated"
  )

  set.seed(7)
  new <- spaccBeta(species, coords, tree = phylo_dist, n_seeds = 3,
                   parallel = FALSE, progress = FALSE)
  set.seed(7)
  old <- suppressWarnings(
    spaccBetaPhylo(species, coords, phylo_dist, n_seeds = 3,
                   parallel = FALSE, progress = FALSE)
  )
  expect_equal(old$beta_total, new$beta_total)
  expect_equal(old$beta_type, "phylogenetic")
})


test_that("diversityProfileFunc is deprecated and forwards to diversityProfile(traits=)", {
  species <- matrix(rpois(20 * 10, 2), nrow = 20,
                    dimnames = list(NULL, paste0("sp", 1:10)))
  traits <- data.frame(body_size = rnorm(10), row.names = paste0("sp", 1:10))

  expect_warning(diversityProfileFunc(species, traits, q = c(0, 1, 2)), "deprecated")

  new <- diversityProfile(species, traits = traits, q = c(0, 1, 2))
  old <- suppressWarnings(diversityProfileFunc(species, traits, q = c(0, 1, 2)))
  expect_equal(old$per_site, new$per_site)
  expect_equal(old$profile_type, "functional")
})


test_that("diversityProfilePhylo is deprecated and forwards to diversityProfile(tree=)", {
  skip_if_not_installed("ape")
  set.seed(1)
  species <- matrix(rpois(20 * 10, 2), nrow = 20,
                    dimnames = list(NULL, paste0("sp", 1:10)))
  tree <- ape::rcoal(10, tip.label = paste0("sp", 1:10))

  expect_warning(diversityProfilePhylo(species, tree, q = c(0, 1, 2)), "deprecated")

  new <- diversityProfile(species, tree = tree, q = c(0, 1, 2))
  old <- suppressWarnings(diversityProfilePhylo(species, tree, q = c(0, 1, 2)))
  expect_equal(old$per_site, new$per_site)
  expect_equal(old$profile_type, "phylogenetic")
})


# Renamed accumulation front door (wavefront -> spaccWavefront) ---------------

test_that("wavefront is deprecated and forwards to spaccWavefront", {
  skip_on_cran()
  set.seed(1)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 12, 1, 0.35), nrow = 20)

  expect_warning(
    wavefront(species, coords, n_seeds = 3, n_steps = 8, progress = FALSE),
    "deprecated"
  )

  set.seed(5)
  new <- spaccWavefront(species, coords, n_seeds = 4, n_steps = 10, progress = FALSE)
  set.seed(5)
  old <- suppressWarnings(
    wavefront(species, coords, n_seeds = 4, n_steps = 10, progress = FALSE)
  )
  expect_equal(old$curves, new$curves)
  expect_s3_class(new, "spacc_wavefront")
})
