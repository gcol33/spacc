sr_fixture <- function(seed = 10, n = 120) {
  set.seed(seed)
  coords <- data.frame(x = runif(n), y = runif(n))
  design <- data.frame(x = c(0.3, 0.4, 0.5, 0.35, 0.45), y = c(0.3, 0.35, 0.3, 0.45, 0.42))
  centres <- runif(50)
  gradient <- sapply(centres, function(cx) rbinom(n, 1, 0.8 * (abs(coords$x - cx) < 0.12)))
  homog <- sapply(runif(50, 0.02, 0.3), function(pi) rbinom(n, 1, pi))
  list(coords = coords, design = design, gradient = gradient, homog = homog)
}

test_that("shift_rotate SES is small where the permutation SES is large", {
  skip_on_cran()
  f <- sr_fixture()
  sac <- spacc(f$gradient, f$coords, method = "knn", focal_points = f$design, progress = FALSE)
  perm <- ses(sac, f$gradient, null_model = "frequency", n_perm = 49, progress = FALSE, seed = 1)
  sr <- ses(sac, f$gradient, null_model = "shift_rotate", n_perm = 49, progress = FALSE, seed = 1)
  steps <- 10:60
  expect_gt(max(abs(perm$ses[steps])), 10)
  expect_lt(max(abs(sr$ses[steps])), 3)
})

test_that("shift_rotate SES is centred near zero without spatial structure", {
  skip_on_cran()
  f <- sr_fixture()
  sac <- spacc(f$homog, f$coords, method = "knn", focal_points = f$design, progress = FALSE)
  sr <- ses(sac, f$homog, null_model = "shift_rotate", n_perm = 49, progress = FALSE, seed = 3)
  expect_lt(abs(mean(sr$ses[5:100])), 1.5)
})

test_that("moved designs stay inside the study area", {
  f <- sr_fixture()
  set.seed(4)
  hull <- grDevices::chull(f$coords$x, f$coords$y)
  for (i in 1:20) {
    moved <- spacc:::.shift_rotate_design(f$design, f$coords)
    expect_gte(nrow(moved), 3)
    expect_true(all(moved$x >= min(f$coords$x) & moved$x <= max(f$coords$x)))
    expect_true(all(moved$y >= min(f$coords$y) & moved$y <= max(f$coords$y)))
  }
})

test_that("shift_rotate needs a fixed knn design", {
  f <- sr_fixture()
  random_seeds <- spacc(f$homog, f$coords, method = "knn", n_seeds = 5, progress = FALSE)
  walk <- spacc(f$homog, f$coords, method = "nn_walk", n_seeds = 5, progress = FALSE)
  expect_error(ses(random_seeds, f$homog, null_model = "shift_rotate", n_perm = 5,
                   progress = FALSE), "user-supplied focal_points")
  expect_error(ses(walk, f$homog, null_model = "shift_rotate", n_perm = 5,
                   progress = FALSE), "user-supplied focal_points")
})

test_that("ses() runs on spacc_phylo and spacc_func objects", {
  skip_on_cran()
  skip_if_not_installed("ape")
  set.seed(5)
  coords <- data.frame(x = runif(20), y = runif(20))
  tree <- ape::rcoal(12)
  sp <- matrix(rbinom(20 * 12, 1, 0.3), 20, dimnames = list(NULL, tree$tip.label))
  ph <- spaccPhylo(sp, coords, tree, metric = c("pd", "evoheritage"), rho = c(0.5, 2),
                   n_seeds = 3, parallel = FALSE, progress = FALSE)
  s_ph <- ses(ph, sp, metric = "evoheritage_rho2", n_perm = 9, progress = FALSE, seed = 1)
  expect_equal(s_ph$n_valid, 9)
  expect_equal(s_ph$metric, "evoheritage_rho2")

  traits <- matrix(rnorm(24), 12, dimnames = list(tree$tip.label, NULL))
  fn <- spaccFunc(sp, coords, traits, metric = "K", n_seeds = 3, parallel = FALSE,
                  progress = FALSE)
  s_fn <- ses(fn, sp, metric = "K", n_perm = 9, progress = FALSE, seed = 1)
  expect_equal(s_fn$n_valid, 9)
})
