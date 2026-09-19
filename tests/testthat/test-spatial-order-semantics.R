test_that("knn ranks every site from one fixed focal point", {
  coords <- data.frame(
    x = c(0, 1, 2, 0),
    y = c(0, 0, 0, 1.5)
  )
  focus <- data.frame(x = 0, y = 0)

  observed <- spacc:::cpp_focal_order(
    coords$x, coords$y, focus$x, focus$y, "euclidean"
  )
  reference <- order(sqrt((coords$x - focus$x)^2 + (coords$y - focus$y)^2)) - 1L

  expect_equal(as.integer(observed[1, ]), reference)
  expect_equal(as.integer(observed[1, ]), c(0L, 1L, 3L, 2L))
})

test_that("nn_walk recursively uses the current site", {
  coords <- data.frame(
    x = c(0, 1, 2, 0),
    y = c(0, 0, 0, 1.5)
  )
  distances <- spacc:::cpp_distance_matrix(coords$x, coords$y, "euclidean")

  observed <- spacc:::cpp_nn_walk_order(distances, 0L)

  expect_equal(as.integer(observed[1, ]), c(0L, 1L, 2L, 3L))
})

test_that("public knn uses the supplied fixed focal point", {
  coords <- data.frame(
    x = c(0, 1, 2, 0),
    y = c(0, 0, 0, 1.5)
  )
  species <- rbind(
    c(1, 0, 0),
    c(0, 1, 0),
    c(0, 0, 1),
    c(0, 1, 0)
  )

  fixed_focus <- spacc(
    species, coords, method = "knn",
    focal_points = data.frame(x = 0, y = 0),
    parallel = FALSE, progress = FALSE
  )
  expect_equal(as.integer(fixed_focus$curves[1, ]), c(1L, 2L, 2L, 3L))
  expect_equal(fixed_focus$focal_points, data.frame(x = 0, y = 0))
})

test_that("knn agrees with an independent observed-focus sSBR definition", {
  coords <- data.frame(
    x = c(0, 1, 3, 4, 2),
    y = c(0, 2, 1, 4, 3)
  )
  species <- rbind(
    c(1, 0, 0, 0),
    c(1, 1, 0, 0),
    c(0, 1, 1, 0),
    c(0, 0, 1, 1),
    c(1, 0, 0, 1)
  )
  foci <- coords[c(1, 4), , drop = FALSE]
  reference <- t(vapply(seq_len(nrow(foci)), function(i) {
    ordering <- order((coords$x - foci$x[i])^2 + (coords$y - foci$y[i])^2)
    cumulative <- integer(nrow(coords))
    seen <- rep(FALSE, ncol(species))
    for (step in seq_along(ordering)) {
      seen <- seen | species[ordering[step], ] > 0
      cumulative[step] <- sum(seen)
    }
    cumulative
  }, integer(nrow(coords))))

  observed <- spacc(
    species, coords, method = "knn", focal_points = foci,
    parallel = FALSE, progress = FALSE
  )

  expect_equal(observed$curves, reference)
})

test_that("haversine focal points are sampled in an equal-area hull", {
  coords <- data.frame(
    x = c(179.5, -179.5, -179.5, 179.5),
    y = c(10, 10, 11, 11)
  )

  set.seed(42)
  foci <- spacc:::.sample_focal_points(coords, 20, distance = "haversine")

  expect_equal(nrow(foci), 20L)
  expect_true(all(is.finite(as.matrix(foci))))
  expect_true(all(foci$y >= 9.9 & foci$y <= 11.1))
  expect_true(all(abs(foci$x) >= 179.4))
})

test_that("default spatial method is canonical knn", {
  coords <- data.frame(x = c(0, 1, 0), y = c(0, 0, 1))
  species <- diag(3)

  default <- spacc(species, coords, n_seeds = 3, seed = 42,
                   parallel = FALSE, progress = FALSE)
  explicit <- spacc(species, coords, method = "knn", n_seeds = 3, seed = 42,
                    parallel = FALSE, progress = FALSE)

  expect_equal(default$curves, explicit$curves)
  expect_equal(default$focal_points, explicit$focal_points)
  expect_equal(default$method, "knn")
})

test_that("radius is consolidated into fixed-focus knn", {
  coords <- data.frame(x = 1:4, y = 1:4)
  species <- diag(4)

  expect_error(
    spacc(species, coords, method = "radius", progress = FALSE),
    "should be one of"
  )
})
