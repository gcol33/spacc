# The spatially constrained curve of Rarefy::directionalSAC() with a site
# distance matrix as gradient: for every plot, the collector curve of the
# sites sorted by distance from it, averaged over plots. Written out here so
# the equivalence with spacc's fixed-focus traversal is checked without the
# (archived) Rarefy package; dev_notes/rarefy_crosscheck.R runs Rarefy itself.
directional_scr <- function(community, coords) {
  d <- as.matrix(stats::dist(coords))
  curves <- vapply(seq_len(nrow(community)), function(i) {
    ord <- order(d[, i])
    seen <- apply(community[ord, , drop = FALSE] > 0, 2, cumsum) > 0
    rowSums(seen)
  }, numeric(nrow(community)))
  rowMeans(curves)
}

test_that("knn with every site as focal point is Rarefy's directional SAC", {
  set.seed(9)
  coords <- data.frame(x = runif(30), y = runif(30))
  community <- matrix(rbinom(30 * 25, 1, 0.2), nrow = 30)
  sac <- spacc(community, coords, method = "knn", focal_points = coords, progress = FALSE)
  expect_equal(unname(colMeans(sac$curves)), directional_scr(community, coords))
})

test_that("as_spacc converts numeric curves and data frame columns", {
  curve <- c(4, 7, 9, 10, 10.5)
  a <- as_spacc(curve)
  expect_s3_class(a, "spacc")
  expect_equal(a$curves, matrix(curve, nrow = 1))
  expect_equal(a$n_sites, 5)

  scr <- data.frame(N_SCR = curve, N_Exact = curve + 1)
  expect_equal(as_spacc(scr)$curves, matrix(curve, nrow = 1))
  expect_equal(as_spacc(scr, column = "N_Exact")$curves, matrix(curve + 1, nrow = 1))
  rare <- data.frame(Rarefaction = curve, IC_up = curve + 1, IC_low = curve - 1)
  expect_equal(as_spacc(rare)$curves, matrix(curve, nrow = 1))
  expect_error(as_spacc(data.frame(a = 1:3)), "column")
  expect_s3_class(plot(as_spacc(scr)), "ggplot")
})
