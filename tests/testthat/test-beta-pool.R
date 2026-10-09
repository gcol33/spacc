pool_fixture <- function(seed = 2) {
  set.seed(seed)
  coords <- data.frame(x = runif(24), y = runif(24))
  species <- matrix(rbinom(24 * 30, 1, 0.25), nrow = 24,
                    dimnames = list(NULL, paste0("sp", 1:30)))
  list(coords = coords, species = species)
}

test_that("one site per pool reproduces the site-level partition", {
  f <- pool_fixture()
  b <- spaccBeta(f$species, f$coords, n_seeds = 4, parallel = FALSE, progress = FALSE,
                 seed = 1, pool = seq_len(24))
  expect_equal(b$pooled$beta_total, b$beta_total)
  expect_equal(b$pooled$beta_turnover, b$beta_turnover)
  expect_equal(b$pooled$beta_nestedness, b$beta_nestedness)
  expect_equal(b$pooled$distance, b$distance)
})

test_that("pooled partition has one step per additional pool", {
  f <- pool_fixture()
  habitat <- rep(c("forest", "edge", "field"), each = 8)
  b <- spaccBeta(f$species, f$coords, n_seeds = 5, parallel = FALSE, progress = FALSE,
                 seed = 1, pool = habitat)
  expect_equal(b$pooled$n_pools, 3)
  expect_equal(dim(b$pooled$beta_total), c(5, 2))
  expect_true(all(b$pooled$beta_total >= 0 & b$pooled$beta_total <= 1))
  expect_equal(b$pooled$beta_total, b$pooled$beta_turnover + b$pooled$beta_nestedness)
  expect_equal(nrow(summary(b, resolution = "pooled")), 2)
  expect_output(print(b), "3 pools")
  expect_s3_class(plot(b, resolution = "both"), "ggplot")
  expect_s3_class(plot(b, resolution = "pooled", partition = FALSE), "ggplot")
})

test_that("pools of identical composition give zero pooled beta", {
  f <- pool_fixture()
  # every pool holds the same set of species
  sp <- f$species
  sp[, ] <- 0L
  sp[cbind(1:24, rep(1:8, 3))] <- 1L
  habitat <- rep(1:3, each = 8)
  b <- spaccBeta(sp, f$coords, n_seeds = 3, parallel = FALSE, progress = FALSE,
                 seed = 1, pool = habitat)
  expect_equal(b$pooled$beta_total, matrix(0, 3, 2), ignore_attr = TRUE)
  expect_true(any(b$beta_total > 0))
})

test_that("pool works with the functional partition and is validated", {
  f <- pool_fixture()
  traits <- matrix(rnorm(60), nrow = 30, dimnames = list(colnames(f$species), NULL))
  b <- spaccBeta(f$species, f$coords, traits = traits, n_seeds = 2, parallel = FALSE,
                 progress = FALSE, pool = rep(1:4, each = 6))
  expect_equal(dim(b$pooled$beta_total), c(2, 3))
  expect_error(spaccBeta(f$species, f$coords, pool = rep(1, 24), progress = FALSE),
               "at least two")
  expect_error(spaccBeta(f$species, f$coords, pool = 1:5, progress = FALSE),
               "one value per site")
  expect_error(summary(spaccBeta(f$species, f$coords, n_seeds = 2, progress = FALSE),
                       resolution = "pooled"), "pool")
})
