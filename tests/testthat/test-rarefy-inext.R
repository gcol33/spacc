rare_fixture <- function(seed = 4) {
  set.seed(seed)
  matrix(rpois(30 * 25, rep(rlnorm(25, 0, 1.2), each = 30)), nrow = 30)
}

test_that("rarefy() equals the iNEXT.3D rarefied Hill numbers at every order", {
  skip_if_not_installed("iNEXT.3D")
  i3 <- asNamespace("iNEXT.3D")
  x <- rare_fixture()
  totals <- colSums(x)
  totals <- totals[totals > 0]
  sizes <- c(1, 5, 20, 60, floor(sum(totals) / 2), sum(totals))
  for (q in c(0, 1, 2)) {
    r <- rarefy(x, n_individuals = sizes, q = q, n_boot = 2)
    expect_equal(r$expected, i3$TD.m.est(totals, sizes, q), tolerance = 1e-8,
                 label = sprintf("q = %g", q))
  }
})

test_that("rarefied Shannon and Simpson grow with sample size", {
  x <- rare_fixture()
  for (q in c(1, 2)) {
    r <- rarefy(x, n_individuals = c(5, 20, 80, 200), q = q, n_boot = 2)
    expect_true(all(diff(r$expected) > 0))
  }
})

test_that("extrapolateCoverage() reads the spatial curve within and the pooled sample beyond", {
  skip_if_not_installed("iNEXT.3D")
  i3 <- asNamespace("iNEXT.3D")
  # three singleton species keep the pooled coverage below 1
  singles <- matrix(0L, 30, 3)
  singles[cbind(c(3, 11, 20), 1:3)] <- 1L
  x <- cbind(rare_fixture(), singles)
  coords <- data.frame(x = runif(30), y = runif(30))
  hc <- spaccHillCoverage(x, coords, q = c(0, 1, 2), n_seeds = 4, progress = FALSE,
                          seed = 1)
  # every seed ends with the same pooled sample and coverage
  final_C <- hc$coverage[1, ncol(hc$coverage)]
  expect_true(final_C < 1)
  beyond <- (final_C + 1) / 2
  totals <- colSums(x)
  totals <- totals[totals > 0]
  for (q in c(0, 1, 2)) {
    ext <- extrapolateCoverage(hc, target_coverage = c(final_C * 0.9, beyond, 1), q = q)
    expect_equal(unname(ext$richness[, 2]), rep(i3$invChat.Ind(totals, q, beyond)$qTD, 4),
                 tolerance = 1e-4)
    expect_equal(unname(ext$richness[, 3]), rep(i3$Diversity_profile(totals, q), 4),
                 tolerance = 1e-8)
    within <- vapply(seq_len(4), function(s) {
      spacc:::interpolate_at_coverage(hc$hills[[q + 1]][s, ], hc$coverage[s, ], final_C * 0.9)
    }, numeric(1))
    expect_equal(unname(ext$richness[, 1]), within)
  }
})

test_that("extrapolateCoverage() refuses orders a spacc_coverage object does not carry", {
  x <- rare_fixture()
  coords <- data.frame(x = runif(30), y = runif(30))
  cv <- spaccCoverage(x, coords, n_seeds = 2, progress = FALSE, seed = 1)
  expect_error(extrapolateCoverage(cv, 0.99, q = 1), "richness only")
  set.seed(9)
  sparse <- matrix(rpois(30 * 60, 0.08), nrow = 30)
  cv_chiu <- spaccCoverage(sparse, coords, n_seeds = 2, progress = FALSE, seed = 1,
                           coverage = "chiu")
  expect_lt(max(cv_chiu$coverage[, 30]), 1)
  expect_error(extrapolateCoverage(cv_chiu, 1), "coverage = \"chao\"")
})
