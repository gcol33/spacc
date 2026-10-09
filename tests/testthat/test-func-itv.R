itv_fixture <- function(seed = 4) {
  set.seed(seed)
  n_sites <- 30
  coords <- data.frame(x = seq(0, 1, length.out = n_sites), y = runif(n_sites))
  species <- matrix(rpois(n_sites * 10, 1.5), nrow = n_sites,
                    dimnames = list(NULL, paste0("sp", 1:10)))
  means <- data.frame(species = colnames(species), size = rnorm(10), beak = rnorm(10))
  list(coords = coords, species = species, means = means)
}

test_that("one row per species reproduces the species-mean curves", {
  f <- itv_fixture()
  traits <- as.matrix(f$means[, c("size", "beak")])
  rownames(traits) <- f$means$species
  a <- spaccFunc(f$species, f$coords, traits, metric = c("fdis", "fric", "rao", "K"),
                 n_seeds = 3, parallel = FALSE, progress = FALSE, seed = 1)
  b <- spaccFunc(f$species, f$coords, f$means, metric = c("fdis", "fric", "rao", "K"),
                 n_seeds = 3, parallel = FALSE, progress = FALSE, seed = 1, itv = TRUE)
  expect_equal(b$curves, a$curves)
  expect_true(b$itv)
})

test_that("site rows split their species' abundance and site-less rows fill the rest", {
  f <- itv_fixture()
  x <- f$species
  long <- rbind(
    data.frame(species = "sp1", site = c(3, 3), size = c(1, 2)),
    data.frame(species = "sp1", site = NA, size = 0),
    data.frame(species = "sp2", site = NA, size = 5)
  )
  rows <- spacc:::.itv_rows(x[, 1:2], long, f$coords)
  expect_equal(rows$weights[3, 1:2], rep(x[3, "sp1"] / 2, 2), ignore_attr = TRUE)
  expect_equal(rows$weights[-3, 3], x[-3, "sp1"], ignore_attr = TRUE)
  expect_equal(rows$weights[3, 3], 0)
  expect_equal(rows$weights[, 4], x[, "sp2"], ignore_attr = TRUE)
})

test_that("a widespread species varying along a gradient adds trait range", {
  set.seed(7)
  n_sites <- 40
  coords <- data.frame(x = seq(0, 1, length.out = n_sites), y = 0)
  # sp1 occurs everywhere, its size tracks x; four local species near the centre
  x <- cbind(sp1 = rep(3, n_sites),
             matrix(rpois(n_sites * 4, 1), n_sites, dimnames = list(NULL, paste0("sp", 2:5))))
  local_size <- c(sp2 = 0.45, sp3 = 0.5, sp4 = 0.55, sp5 = 0.5)
  long <- rbind(data.frame(species = "sp1", site = seq_len(n_sites), size = coords$x),
                data.frame(species = names(local_size), site = NA, size = local_size))
  means <- rbind(data.frame(species = "sp1", site = NA, size = mean(coords$x)),
                 long[long$species != "sp1", ])
  itv <- spaccFunc(x, coords, long, metric = c("fdis", "fric"), n_seeds = 5,
                   parallel = FALSE, progress = FALSE, seed = 2, itv = TRUE)
  mean_curve <- spaccFunc(x, coords, means, metric = c("fdis", "fric"), n_seeds = 5,
                          parallel = FALSE, progress = FALSE, seed = 2, itv = TRUE)
  final <- function(r, m) mean(r$curves[[m]][, n_sites])
  expect_gt(final(itv, "fric"), final(mean_curve, "fric"))
  expect_gt(final(itv, "fdis"), final(mean_curve, "fdis"))
})

test_that("ITV traits are validated", {
  f <- itv_fixture()
  expect_error(spaccFunc(f$species, f$coords, f$means[, -1], itv = TRUE, progress = FALSE),
               "species")
  bad <- f$means
  bad$species[1] <- "unknown"
  expect_error(spaccFunc(f$species, f$coords, bad, itv = TRUE, progress = FALSE),
               "Species not in x")
  expect_error(spaccFunc(f$species, f$coords, data.frame(species = "sp1", site = 99, size = 1),
                         itv = TRUE, progress = FALSE), "not found in coords")
})
