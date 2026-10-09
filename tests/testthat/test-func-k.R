k_fixture <- function(seed = 8) {
  set.seed(seed)
  n_sp <- 15
  traits <- cbind(size = runif(n_sp), beak = rnorm(n_sp))
  rownames(traits) <- paste0("sp", seq_len(n_sp))
  species <- matrix(rpois(25 * n_sp, 1.5), nrow = 25, dimnames = list(NULL, rownames(traits)))
  list(species = species, traits = traits,
       coords = data.frame(x = runif(25), y = runif(25)))
}

test_that("metric = 'K' returns the four facets and K, all in [0, 1]", {
  f <- k_fixture()
  res <- spaccFunc(f$species, f$coords, f$traits, metric = c("fdis", "K"),
                   n_seeds = 4, parallel = FALSE, progress = FALSE, seed = 1)
  expect_equal(res$metric, c("fdis", "k_richness", "k_biomass_evenness",
                             "k_trait_evenness", "k_dispersion", "K"))
  for (m in res$metric[-1]) {
    v <- res$curves[[m]]
    expect_true(all(is.na(v) | (v >= 0 & v <= 1 + 1e-12)), label = m)
  }
  expect_s3_class(summary(res), "data.frame")
  expect_s3_class(plot(res), "ggplot")
})

test_that("final-step K equals calc_k_facets on the pooled community", {
  f <- k_fixture()
  res <- spaccFunc(f$species, f$coords, f$traits, metric = "K",
                   n_seeds = 3, parallel = FALSE, progress = FALSE, seed = 2)
  L <- apply(f$traits, 2, function(v) diff(range(v)))
  ref <- spacc:::calc_k_facets(f$traits, colSums(f$species), L)
  for (m in names(ref)) {
    expect_equal(res$curves[[m]][, nrow(f$species)], rep(unname(ref[m]), 3), label = m)
  }
})

test_that("K is the geometric mean of the per-trait K values", {
  f <- k_fixture()
  ab <- colSums(f$species)
  L <- apply(f$traits, 2, function(v) diff(range(v)))
  both <- spacc:::calc_k_facets(f$traits, ab, L)
  per_axis <- vapply(1:2, function(j) {
    spacc:::calc_k_facets(f$traits[, j, drop = FALSE], ab, L[j])[["K"]]
  }, numeric(1))
  expect_equal(unname(both["K"]), sqrt(prod(per_axis)))
})

test_that("axis_length_k scales the richness facet only", {
  f <- k_fixture()
  a <- spaccFunc(f$species, f$coords, f$traits[, 1, drop = FALSE], metric = "K",
                 n_seeds = 2, parallel = FALSE, progress = FALSE, seed = 3)
  b <- spaccFunc(f$species, f$coords, f$traits[, 1, drop = FALSE], metric = "K",
                 n_seeds = 2, parallel = FALSE, progress = FALSE, seed = 3,
                 axis_length_k = 2 * diff(range(f$traits[, 1])))
  expect_equal(b$curves$k_richness, a$curves$k_richness / 2)
  expect_equal(b$curves$k_dispersion, a$curves$k_dispersion)
})

test_that("a constant trait is rejected for K", {
  f <- k_fixture()
  tr <- cbind(f$traits, flat = 1)
  expect_error(spaccFunc(f$species, f$coords, tr, metric = "K", n_seeds = 2,
                         progress = FALSE), "positive length")
})
