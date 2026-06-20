test_that("as_spacc maps a permutation specaccum onto seed curves", {
  # Hand-built specaccum (method = "random"): sites x permutations
  perm <- matrix(c(2, 4, 5,
                   3, 4, 5,
                   1, 3, 5), nrow = 3, byrow = FALSE)
  fake <- structure(
    list(method = "random", sites = 1:3,
         richness = rowMeans(perm), perm = perm),
    class = "specaccum"
  )

  sp <- as_spacc(fake)

  expect_s3_class(sp, "spacc")
  # curves are permutations-by-sites: the transpose of perm
  expect_equal(sp$curves, t(perm), ignore_attr = TRUE)
  expect_equal(sp$n_sites, 3L)
  expect_equal(sp$n_seeds, 3L)
  # final richness across permutations is the species pool
  expect_equal(sp$n_species, 5L)
  expect_equal(sp$method, "random")
  expect_null(sp$coords)
})

test_that("as_spacc handles a single mean curve (no perm matrix)", {
  fake <- structure(
    list(method = "collector", sites = 1:3, richness = c(2, 4, 5)),
    class = "specaccum"
  )
  sp <- as_spacc(fake)
  expect_equal(nrow(sp$curves), 1L)
  expect_equal(sp$n_sites, 3L)
  expect_equal(sp$n_species, 5L)
})

test_that("converted object works with summary, print, and compare", {
  perm_a <- matrix(rpois(5 * 6, 4), nrow = 5)  # sites x perms
  perm_a <- t(apply(perm_a, 1, cumsum))[, seq_len(6)]  # monotone-ish per perm
  fake_a <- structure(list(method = "random", sites = 1:5,
                           richness = rowMeans(perm_a), perm = perm_a),
                      class = "specaccum")
  fake_b <- structure(list(method = "random", sites = 1:5,
                           richness = rowMeans(perm_a) + 1, perm = perm_a + 1),
                      class = "specaccum")

  a <- as_spacc(fake_a)
  b <- as_spacc(fake_b)

  expect_s3_class(summary(a), "summary.spacc")
  expect_output(print(a), "spacc")

  cmp <- compare(a, b, n_perm = 99L)
  expect_s3_class(cmp, "spacc_comp")
})

test_that("as_spacc errors clearly on unsupported input", {
  expect_error(as_spacc(1:5), "specaccum")
  expect_error(as_spacc(data.frame(a = 1)), "specaccum")
})

test_that("as_spacc round-trips a real vegan::specaccum object", {
  skip_if_not_installed("vegan")

  m <- matrix(rbinom(20 * 8, 1, 0.4), nrow = 20, ncol = 8)
  classic <- vegan::specaccum(m, method = "random", permutations = 15)

  sp <- as_spacc(classic)
  expect_s3_class(sp, "spacc")
  expect_equal(sp$n_sites, 20L)
  expect_equal(sp$n_seeds, 15L)
  expect_equal(sp$curves, t(classic$perm), ignore_attr = TRUE)
  expect_no_error(summary(sp))

  skip_if_not_installed("ggplot2")
  expect_s3_class(plot(sp), "ggplot")
})
