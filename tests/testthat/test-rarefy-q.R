test_that("rarefy accepts arbitrary q and recovers Hill number at full sample", {
  set.seed(31)
  species <- matrix(rpois(20 * 12, 3), nrow = 20)
  totals <- colSums(species); totals <- totals[totals > 0]
  N <- sum(totals); p <- totals / N

  for (q in c(0.5, 3)) {
    res <- rarefy(species, q = q, n_boot = 10)
    expect_s3_class(res, "spacc_rare")
    hill_q <- sum(p^q)^(1 / (1 - q))
    expect_equal(unname(res$expected[length(res$expected)]), hill_q, tolerance = 1e-6)
    expect_true(all(is.finite(res$expected)))
  }
})


test_that("rarefy q=3 no longer falls back to richness", {
  set.seed(32)
  species <- matrix(rpois(20 * 12, 3), nrow = 20)
  r0 <- rarefy(species, q = 0, n_boot = 5)
  r3 <- rarefy(species, q = 3, n_boot = 5)
  expect_false(isTRUE(all.equal(r3$expected, r0$expected)))
})


test_that("rarefy rejects negative q", {
  species <- matrix(rpois(10 * 6, 3), nrow = 10)
  expect_error(rarefy(species, q = -1), "non-negative")
})
