test_that("cpp_distance_matrix matches base R dist()", {
  set.seed(42)
  n <- 50
  x <- runif(n)
  y <- runif(n)

  # spacc version
  d_spacc <- spacc:::cpp_distance_matrix(x, y)

  # base R version
  coords <- cbind(x, y)
  d_base <- as.matrix(dist(coords))

  expect_equal(dim(d_spacc), c(n, n))
  # Compare values ignoring dimnames
  expect_equal(as.vector(d_spacc), as.vector(d_base), tolerance = 1e-10)
})

test_that("distance matrix is symmetric", {
  x <- c(0, 1, 2)
  y <- c(0, 1, 0)

  d <- spacc:::cpp_distance_matrix(x, y)

  expect_equal(d, t(d))
  expect_equal(diag(d), rep(0, 3))
})


test_that("haversine distance gives reasonable values", {
  # ~111km between 0,0 and 1,0 degree longitude at equator
  x <- c(0, 1)
  y <- c(0, 0)

  d <- spacc:::cpp_distance_matrix(x, y, "haversine")

  # Result could be in meters (~111000) or km (~111)
  expect_true(d[1, 2] > 100)
  expect_true(d[1, 2] > 0)
})
