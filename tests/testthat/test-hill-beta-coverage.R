beta3d_fixture <- function(n_sites = 12, n_species = 40, seed = 5) {
  set.seed(seed)
  coords <- data.frame(x = runif(n_sites), y = runif(n_sites))
  lambda <- rlnorm(n_species, 0, 1.2)
  species <- t(vapply(seq_len(n_sites), function(i) {
    rpois(n_species, lambda * runif(n_species, 0.2, 1.8))
  }, numeric(n_species)))
  storage.mode(species) <- "integer"
  list(species = species, coords = coords)
}

test_that("calc_hill_at_size matches iNEXT.3D rarefaction and extrapolation", {
  skip_if_not_installed("iNEXT.3D")
  i3 <- asNamespace("iNEXT.3D")
  set.seed(11)
  x <- as.integer(round(rlnorm(60, 1.5, 1.4)))
  x <- x[x > 0]
  n <- sum(x)
  m <- c(1, 5.5, floor(n / 3), n / 2 + 0.3, n - 1, n, n + 7, 2 * n)
  q <- c(0, 1, 2)
  ours <- spacc:::calc_hill_at_size(x, m, q)
  ref <- t(sapply(q, function(qq) i3$TD.m.est(x, m, qq)))
  expect_equal(ours$hill, ref, tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(ours$coverage, i3$Coverage(x, "abundance", m), tolerance = 1e-10)
})

test_that("calc_hill_at_coverage matches iNEXT.3D invChat", {
  skip_if_not_installed("iNEXT.3D")
  i3 <- asNamespace("iNEXT.3D")
  set.seed(12)
  x <- as.integer(round(rlnorm(50, 1.2, 1.3)))
  x <- x[x > 0]
  C <- c(0.4, 0.8, 0.97, 0.995)
  ref <- i3$invChat.Ind(x, c(0, 1, 2), C)
  ours <- spacc:::calc_hill_at_coverage(x, C, c(0, 1, 2))
  expect_equal(as.vector(t(ours$hill)), ref$qTD, tolerance = 1e-4)
})

test_that("standardize = 'none' is unchanged", {
  f <- beta3d_fixture()
  a <- spaccHillBeta(f$species, f$coords, n_seeds = 4, parallel = FALSE,
                     progress = FALSE, seed = 1)
  b <- spaccHillBeta(f$species, f$coords, n_seeds = 4, parallel = FALSE,
                     progress = FALSE, seed = 1, standardize = "none")
  expect_identical(a$beta, b$beta)
  expect_null(a$dissimilarity)
  expect_null(a$coverage)
})

final_step <- function(res, comp) {
  vapply(res[[comp]], function(m) m[1, ncol(m)], numeric(1))
}

beta3d_reference <- function(species, level) {
  out <- iNEXT.beta3D::iNEXTbeta3D(list(d = t(species)), diversity = "TD", q = c(0, 1, 2),
                                   datatype = "abundance", base = "coverage",
                                   level = level, nboot = 0)[[1]]
  meta <- c("Dataset", "Order.q", "SC", "Size", "Method", "s.e.", "LCL", "UCL", "Diversity")
  pick <- function(tab) tab[[setdiff(names(tab), meta)]][order(tab$Order.q)]
  list(gamma = pick(out$gamma), alpha = pick(out$alpha), beta = pick(out$beta),
       `1-CqN` = pick(out$`1-C`), `1-UqN` = pick(out$`1-U`),
       `1-VqN` = pick(out$`1-V`), `1-SqN` = pick(out$`1-S`))
}

test_that("final-step standardised beta matches iNEXT.beta3D at the same coverage", {
  skip_on_cran()
  skip_if_not_installed("iNEXT.beta3D")
  f <- beta3d_fixture()
  res <- spaccHillBeta(f$species, f$coords, n_seeds = 1, parallel = FALSE,
                       progress = FALSE, seed = 1, standardize = "coverage")
  C <- res$coverage[1, ncol(res$coverage)]
  ref <- beta3d_reference(f$species, C)
  for (comp in c("gamma", "alpha", "beta")) {
    expect_equal(unname(final_step(res, comp)), ref[[comp]], tolerance = 1e-4,
                 label = comp)
  }
  for (comp in names(res$dissimilarity)) {
    expect_equal(unname(final_step(res$dissimilarity, comp)), ref[[comp]],
                 tolerance = 1e-4, label = comp)
  }
})

test_that("fixed coverage, including extrapolation, matches iNEXT.beta3D", {
  skip_on_cran()
  skip_if_not_installed("iNEXT.beta3D")
  f <- beta3d_fixture()
  for (level in c(0.8, 0.99)) {
    res <- spaccHillBeta(f$species, f$coords, n_seeds = 1, parallel = FALSE,
                         progress = FALSE, seed = 1, standardize = "coverage",
                         coverage = level)
    ref <- beta3d_reference(f$species, level)
    expect_equal(unname(final_step(res, "beta")), ref$beta, tolerance = 1e-4,
                 label = sprintf("beta at C = %g", level))
    expect_true(all(res$coverage == level))
  }
})

test_that("standardised beta is 1 at the first step and dissimilarities are NA there", {
  f <- beta3d_fixture()
  res <- spaccHillBeta(f$species, f$coords, n_seeds = 3, parallel = FALSE,
                       progress = FALSE, seed = 2, standardize = "coverage")
  for (b in res$beta) expect_equal(b[, 1], rep(1, 3))
  for (d in res$dissimilarity) for (m in d) expect_true(all(is.na(m[, 1])))
  s <- summary(res)
  expect_setequal(unique(s$component),
                  c("gamma", "alpha", "beta", "1-CqN", "1-UqN", "1-VqN", "1-SqN"))
  expect_s3_class(plot(res, component = "1-CqN"), "ggplot")
  expect_s3_class(plot(res), "ggplot")
})

test_that("identical sites give standardised beta near 1 throughout", {
  set.seed(3)
  row <- rpois(30, 8)
  species <- matrix(rep(row, each = 10), nrow = 10)
  coords <- data.frame(x = runif(10), y = runif(10))
  res <- spaccHillBeta(species, coords, n_seeds = 2, parallel = FALSE, progress = FALSE,
                       standardize = "coverage", seed = 1)
  for (b in res$beta) expect_true(all(abs(b - 1) < 0.05))
})

test_that("coverage arguments are validated", {
  f <- beta3d_fixture()
  expect_error(spaccHillBeta(f$species, f$coords, coverage = 0.9, progress = FALSE),
               "requires standardize")
  expect_error(spaccHillBeta(f$species / 2, f$coords, standardize = "coverage",
                             progress = FALSE), "integer abundances")
  expect_error(spaccHillBeta(f$species, f$coords, q = 0.5, standardize = "coverage",
                             coverage = 0.9, progress = FALSE), "integer q")
})
