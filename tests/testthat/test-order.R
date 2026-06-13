test_that("user-defined order reproduces collector when order is identity", {
  set.seed(1)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 12, 1, 0.3), nrow = 20)

  collector <- spacc(species, coords, method = "collector", progress = FALSE)
  user <- spacc(species, coords, order = seq_len(20), progress = FALSE)

  expect_equal(unname(user$curves[1, ]), unname(collector$curves[1, ]))
  expect_equal(user$method, "user")
})


test_that("a single ordering yields one curve ending at total richness", {
  set.seed(2)
  coords <- data.frame(x = runif(25), y = runif(25))
  species <- matrix(rbinom(25 * 15, 1, 0.3), nrow = 25)

  ord <- sample(25)
  res <- spacc(species, coords, order = ord, progress = FALSE)

  expect_s3_class(res, "spacc")
  expect_equal(res$n_seeds, 1L)
  expect_equal(nrow(res$curves), 1L)
  expect_equal(res$curves[1, 25], sum(colSums(species) > 0))
  # Curve is non-decreasing
  expect_true(all(diff(res$curves[1, ]) >= 0))
})


test_that("multiple orderings produce one curve each with shared endpoint", {
  set.seed(3)
  coords <- data.frame(x = runif(18), y = runif(18))
  species <- matrix(rbinom(18 * 10, 1, 0.4), nrow = 18)

  orders <- t(replicate(8, sample(18)))
  res <- spacc(species, coords, order = orders, progress = FALSE)

  expect_equal(res$n_seeds, 8L)
  expect_equal(dim(res$curves), c(8L, 18L))
  total <- sum(colSums(species) > 0)
  expect_true(all(res$curves[, 18] == total))
})


test_that("order accepts a list of vectors", {
  set.seed(4)
  coords <- data.frame(x = runif(12), y = runif(12))
  species <- matrix(rbinom(12 * 8, 1, 0.4), nrow = 12)

  res <- spacc(species, coords,
               order = list(seq_len(12), rev(seq_len(12))),
               progress = FALSE)

  expect_equal(res$n_seeds, 2L)
  expect_equal(dim(res$curves), c(2L, 12L))
})


test_that("invalid orderings are rejected", {
  coords <- data.frame(x = runif(10), y = runif(10))
  species <- matrix(rbinom(10 * 6, 1, 0.4), nrow = 10)

  # Wrong length
  expect_error(spacc(species, coords, order = 1:5, progress = FALSE),
               "site indices per ordering")
  # Not a permutation (duplicate, missing site)
  expect_error(spacc(species, coords, order = c(1:9, 9), progress = FALSE),
               "permutation")
})


test_that("order is incompatible with support and time", {
  coords <- data.frame(x = runif(10), y = runif(10))
  species <- matrix(rbinom(10 * 6, 1, 0.4), nrow = 10)

  expect_error(
    spacc(species, coords, order = sample(10), time = runif(10), progress = FALSE),
    "cannot be combined with `time`"
  )
})


test_that("summary works on a user-defined accumulation object", {
  set.seed(5)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 12, 1, 0.3), nrow = 20)

  res <- spacc(species, coords, order = t(replicate(10, sample(20))), progress = FALSE)
  s <- summary(res)

  expect_s3_class(s, "summary.spacc")
  expect_equal(length(s$mean), 20L)
  expect_true(all(s$lower <= s$mean + 1e-9))
  expect_true(all(s$upper >= s$mean - 1e-9))
})
