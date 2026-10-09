net_fixture <- function(seed = 1) {
  set.seed(seed)
  coords <- data.frame(x = runif(25), y = runif(25))
  n <- 300
  interactions <- data.frame(
    site = sample(25, n, replace = TRUE),
    lower = paste0("plant", sample(8, n, replace = TRUE, prob = 8:1)),
    upper = paste0("bee", sample(15, n, replace = TRUE)),
    count = rpois(n, 1) + 1L
  )
  list(coords = coords, interactions = interactions)
}

test_that("interaction richness is bounded by the product of the species curves", {
  f <- net_fixture()
  net <- spaccNetwork(f$interactions, f$coords, n_seeds = 4, parallel = FALSE,
                      progress = FALSE, seed = 2)
  bound <- net$lower[["q0.0"]] * net$upper[["q0.0"]]
  expect_true(all(net$interaction[["q0.0"]] <= bound))
  ok <- is.finite(net$connectance)
  expect_true(all(net$connectance[ok] > 0 & net$connectance[ok] <= 1))
  expect_equal(net$n_interactions,
               nrow(unique(f$interactions[, c("lower", "upper")])))
})

test_that("final-step values are the Hill numbers and coverage of the pooled network", {
  f <- net_fixture()
  net <- spaccNetwork(f$interactions, f$coords, q = c(0, 1, 2), n_seeds = 3,
                      parallel = FALSE, progress = FALSE, seed = 3)
  freq <- tapply(f$interactions$count,
                 paste(f$interactions$lower, f$interactions$upper), sum)
  for (i in seq_along(net$q)) {
    expect_equal(net$interaction[[i]][, 25],
                 rep(spacc:::calc_hill_number(as.numeric(freq), net$q[i]), 3))
  }
  expect_equal(net$coverage$interaction[, 25],
               rep(spacc:::calc_coverage(as.integer(freq)), 3))
})

test_that("interaction Hill numbers match iNEXT.3D observed diversity", {
  skip_if_not_installed("iNEXT.3D")
  f <- net_fixture()
  net <- spaccNetwork(f$interactions, f$coords, q = c(0, 1, 2), n_seeds = 2,
                      parallel = FALSE, progress = FALSE, seed = 4)
  freq <- as.numeric(tapply(f$interactions$count,
                            paste(f$interactions$lower, f$interactions$upper), sum))
  ref <- iNEXT.3D::ObsAsy3D(freq, diversity = "TD", q = c(0, 1, 2), datatype = "abundance",
                            method = "Observed", nboot = 0)
  expect_equal(vapply(net$interaction, function(m) m[1, 25], numeric(1)),
               ref$qTD[order(ref$Order.q)], ignore_attr = TRUE)
})

test_that("sites are matched by coords$site and S3 methods work", {
  f <- net_fixture()
  coords <- f$coords
  coords$site <- paste0("s", 25:1)
  inter <- f$interactions
  inter$site <- paste0("s", 26 - inter$site)
  a <- spaccNetwork(f$interactions, f$coords, n_seeds = 2, parallel = FALSE,
                    progress = FALSE, seed = 5)
  b <- spaccNetwork(inter, coords, n_seeds = 2, parallel = FALSE, progress = FALSE, seed = 5)
  expect_equal(a$interaction, b$interaction)
  expect_output(print(a), "interactions between 8 lower and 15 upper")
  expect_s3_class(summary(a), "data.frame")
  for (cmp in c("interaction", "species", "connectance", "coverage")) {
    expect_s3_class(plot(a, component = cmp), "ggplot")
  }
  bad <- f$interactions
  bad$site[1] <- 99
  expect_error(spaccNetwork(bad, f$coords, progress = FALSE), "not found in coords")
})
