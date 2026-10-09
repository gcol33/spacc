dar_fixture <- function(seed = 3) {
  set.seed(seed)
  coords <- data.frame(x = runif(25), y = runif(25))
  species <- matrix(rpois(25 * 15, 2), nrow = 25, dimnames = list(NULL, paste0("sp", 1:15)))
  list(coords = coords, species = species)
}

shoelace_hull <- function(x, y) {
  if (length(x) < 3) return(0)
  h <- grDevices::chull(x, y)
  xs <- x[h]; ys <- y[h]
  abs(sum(xs * c(ys[-1], ys[1]) - c(xs[-1], xs[1]) * ys)) / 2
}

test_that("incremental hull area equals the convex hull of the first k sites", {
  f <- dar_fixture()
  orders <- rbind(sample(25), sample(25))
  area <- spacc:::cpp_area_along_orders(f$coords$x, f$coords$y, orders - 1L, "convex_hull")
  for (s in 1:2) {
    ref <- vapply(1:25, function(k) {
      idx <- orders[s, 1:k]
      shoelace_hull(f$coords$x[idx], f$coords$y[idx])
    }, numeric(1))
    expect_equal(area[s, ], ref, tolerance = 1e-12)
  }
})

test_that("Voronoi cells tile the bounding box", {
  f <- dar_fixture()
  orders <- matrix(sample(25), nrow = 1)
  area <- spacc:::cpp_area_along_orders(f$coords$x, f$coords$y, orders - 1L, "voronoi")
  bbox <- diff(range(f$coords$x)) * diff(range(f$coords$y))
  expect_equal(area[1, 25], bbox, tolerance = 1e-10)
  expect_true(all(diff(area[1, ]) > 0))
  # a regular grid gives equal interior cells
  g <- expand.grid(x = 0:4, y = 0:4)
  ga <- spacc:::cpp_area_along_orders(g$x, g$y, matrix(0:24, nrow = 1), "voronoi")
  cells <- diff(c(0, ga[1, ]))
  interior <- g$x %in% 1:3 & g$y %in% 1:3
  expect_equal(cells[interior], rep(1, 9))
  expect_equal(sum(cells), 16)
})

test_that("dar(matrix) equals dar(spaccHill(matrix))", {
  f <- dar_fixture()
  a <- dar(f$species, f$coords, q = c(0, 1), n_seeds = 3, area_method = "convex_hull",
           parallel = FALSE, progress = FALSE, seed = 4)
  h <- spaccHill(f$species, f$coords, q = c(0, 1), n_seeds = 3, parallel = FALSE,
                 progress = FALSE, seed = 4)
  b <- dar(h, area_method = "convex_hull")
  expect_identical(a$curves, b$curves)
  expect_identical(a$area, b$area)
  expect_equal(a$q, c(0, 1))
  expect_true("q" %in% names(summary(a)))
})

test_that("dar works on functional, phylogenetic, richness and redundancy curves", {
  skip_if_not_installed("ape")
  f <- dar_fixture()
  traits <- matrix(rnorm(30), nrow = 15, dimnames = list(colnames(f$species), NULL))
  func <- spaccFunc(f$species, f$coords, traits, metric = c("fdis", "fric"),
                    n_seeds = 3, parallel = FALSE, progress = FALSE, seed = 5)
  d <- dar(func, area_method = "voronoi")
  expect_equal(d$metric, c("fdis", "fric"))
  expect_equal(d$power_law$metric, c("fdis", "fric"))
  expect_true(all(is.finite(d$power_law$z)))
  expect_true("metric" %in% names(summary(d)))
  expect_s3_class(plot(d, log_scale = TRUE), "ggplot")
  expect_output(print(d), "fdis: z =")

  tree <- ape::rcoal(15, tip.label = colnames(f$species))
  expect_s3_class(dar(spaccPhylo(f$species, f$coords, tree, metric = "pd", n_seeds = 2,
                                 parallel = FALSE, progress = FALSE)), "spacc_dar")
  expect_s3_class(dar(spaccRedundancy(f$species, f$coords, tree = tree, n_seeds = 2,
                                      parallel = FALSE, progress = FALSE)), "spacc_dar")
  sac <- spacc((f$species > 0) * 1, f$coords, n_seeds = 2, progress = FALSE)
  rich <- dar(sac, area_method = "count")
  expect_equal(rich$metric, "richness")
  expect_equal(unname(rich$area[1, ]), 1:25)
})

test_that("the richness power-law exponent is recovered from a known SAR", {
  # richness grows as k^0.5 against an area proportional to k
  f <- dar_fixture()
  d <- dar(f$species, f$coords, q = 0, n_seeds = 2, area_method = "count",
           parallel = FALSE, progress = FALSE, seed = 1)
  d$curves <- list(q0.0 = matrix(rep(sqrt(1:25), each = 2), nrow = 2))
  fit <- spacc:::.power_law(colMeans(d$area), colMeans(d$curves[[1]]))
  expect_equal(unname(fit["z"]), 0.5)
})
