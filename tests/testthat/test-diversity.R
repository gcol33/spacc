test_that("alphaDiversity returns correct structure", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- alphaDiversity(species, q = c(0, 1, 2))

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 20)
  expect_equal(ncol(result), 3)
})


test_that("alphaDiversity q=0 equals site richness", {
  set.seed(42)
  species <- matrix(rpois(15 * 8, 3), nrow = 15)

  result <- alphaDiversity(species, q = 0)

  # q=0 should equal number of species present at each site
  expected <- rowSums(species > 0)
  expect_equal(as.numeric(result[, 1]), expected)
})


test_that("alphaDiversity q=0 >= q=1 >= q=2", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- alphaDiversity(species, q = c(0, 1, 2))

  # Hill number inequality: q0 >= q1 >= q2
  for (i in 1:nrow(result)) {
    expect_true(result[i, 1] >= result[i, 2] - 1e-10,
                info = paste("Row", i, "q0 >= q1"))
    expect_true(result[i, 2] >= result[i, 3] - 1e-10,
                info = paste("Row", i, "q1 >= q2"))
  }
})


test_that("alphaDiversity with coords returns spacc_alpha", {
  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 3), nrow = 15)

  result <- alphaDiversity(species, q = c(0, 1, 2), coords = coords)

  expect_s3_class(result, "spacc_alpha")
})


test_that("gammaDiversity returns named vector", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- gammaDiversity(species, q = c(0, 1, 2))

  expect_true(is.numeric(result))
  expect_length(result, 3)
})


test_that("gammaDiversity q=0 equals total species count", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- gammaDiversity(species, q = 0)

  expected <- sum(colSums(species) > 0)
  expect_equal(as.numeric(result), expected)
})


test_that("gammaDiversity q=0 >= q=1 >= q=2", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- gammaDiversity(species, q = c(0, 1, 2))

  expect_true(result[1] >= result[2] - 1e-10)
  expect_true(result[2] >= result[3] - 1e-10)
})


test_that("diversityPartition returns correct structure", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- diversityPartition(species, q = c(0, 1, 2))

  expect_s3_class(result, "spacc_partition")
})


test_that("diversityPartition gamma = alpha * beta", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- diversityPartition(species, q = c(0, 1, 2))

  # Multiplicative partitioning: gamma = alpha * beta
  for (i in seq_along(result$q)) {
    computed_gamma <- result$alpha[i] * result$beta[i]
    expect_equal(computed_gamma, result$gamma[i], tolerance = 1e-6,
                 info = paste("q =", result$q[i]))
  }
})


test_that("diversityPartition beta >= 1", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- diversityPartition(species, q = c(0, 1, 2))

  expect_true(all(result$beta >= 1 - 1e-10))
})


test_that("diversityPartition print works", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- diversityPartition(species, q = c(0, 1, 2))

  expect_output(print(result), "Diversity Partition")
})


test_that("diversityPartition summary returns data.frame", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- diversityPartition(species, q = c(0, 1, 2))

  summ <- summary(result)
  expect_s3_class(summ, "data.frame")
  expect_true("beta_normalized" %in% names(summ))
})


test_that("alphaDiversity spacc_alpha print works", {
  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 3), nrow = 15)

  result <- alphaDiversity(species, q = c(0, 1, 2), coords = coords)

  expect_output(print(result), "spacc alpha diversity")
})


test_that("alphaDiversity spacc_alpha summary returns data.frame", {
  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 3), nrow = 15)

  result <- alphaDiversity(species, q = c(0, 1, 2), coords = coords)

  summ <- summary(result)
  expect_s3_class(summ, "data.frame")
  expect_true("mean" %in% names(summ))
})


test_that("plot.spacc_alpha histogram works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rpois(15 * 8, 3), nrow = 15)

  result <- alphaDiversity(species, q = c(0, 1, 2), coords = coords)

  p <- plot(result, type = "histogram")
  expect_s3_class(p, "ggplot")
})


test_that("plot.spacc_alpha map errors without coords", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  species <- matrix(rpois(15 * 8, 3), nrow = 15)

  # No coords -> plain data.frame, not spacc_alpha
  result <- alphaDiversity(species, q = c(0, 1))
  expect_false(inherits(result, "spacc_alpha"))
})


test_that("diversityPartition with weights works", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- diversityPartition(species, q = c(0, 1, 2), weights = "proportional")

  expect_s3_class(result, "spacc_partition")
  expect_true(all(result$beta >= 1 - 1e-10))
})


test_that("plot.spacc_partition works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- diversityPartition(species, q = c(0, 1, 2), coords = coords)

  p <- plot(result)
  expect_s3_class(p, "ggplot")
})


# diversityProfile tests --------------------------------------------------

test_that("diversityProfile returns correct structure", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- diversityProfile(species, q = seq(0, 2, by = 0.5))

  expect_s3_class(result, "spacc_profile")
  expect_true(!is.null(result$per_site))
  expect_true(!is.null(result$regional))
  expect_equal(nrow(result$per_site), 20)
  expect_equal(ncol(result$per_site), length(seq(0, 2, by = 0.5)))
})


test_that("diversityProfile D_q is non-increasing in q", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- diversityProfile(species, q = seq(0, 3, by = 0.5))

  # Check per-site: each row should be non-increasing
  for (i in seq_len(nrow(result$per_site))) {
    vals <- result$per_site[i, ]
    for (j in 2:length(vals)) {
      expect_true(vals[j] <= vals[j - 1] + 1e-10,
                  info = sprintf("Site %d, q[%d] vs q[%d]", i, j - 1, j))
    }
  }

  # Check regional
  for (j in 2:length(result$regional)) {
    expect_true(result$regional[j] <= result$regional[j - 1] + 1e-10)
  }
})


test_that("diversityProfile D_0 equals richness", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- diversityProfile(species, q = c(0, 1, 2))

  # D_0 per site = observed richness
  expected <- rowSums(species > 0)
  expect_equal(as.numeric(result$per_site[, 1]), expected)
})


test_that("diversityProfile regional >= mean alpha", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- diversityProfile(species, q = seq(0, 2, by = 0.5))

  mean_alpha <- colMeans(result$per_site)
  # Gamma >= mean alpha is a fundamental property
  for (j in seq_along(result$q)) {
    expect_true(result$regional[j] >= mean_alpha[j] - 1e-10,
                info = sprintf("q = %.1f", result$q[j]))
  }
})


test_that("diversityProfile type = 'per_site' omits regional", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- diversityProfile(species, q = c(0, 1), type = "per_site")
  expect_true(!is.null(result$per_site))
  expect_null(result$regional)
})


test_that("diversityProfile type = 'regional' omits per_site", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- diversityProfile(species, q = c(0, 1), type = "regional")
  expect_null(result$per_site)
  expect_true(!is.null(result$regional))
})


test_that("diversityProfile print works", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)
  result <- diversityProfile(species, q = seq(0, 2, by = 0.5))

  expect_output(print(result), "Diversity profile")
})


test_that("diversityProfile summary returns data.frame", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)
  result <- diversityProfile(species, q = seq(0, 2, by = 0.5))

  summ <- summary(result)
  expect_s3_class(summ, "data.frame")
  expect_true("mean_alpha" %in% names(summ))
  expect_true("gamma" %in% names(summ))
})


test_that("plot.spacc_profile works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)
  result <- diversityProfile(species, q = seq(0, 2, by = 0.5))

  p <- plot(result)
  expect_s3_class(p, "ggplot")
})


test_that("as_sf.spacc_profile works", {
  skip_if_not_installed("sf")

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- diversityProfile(species, q = c(0, 1, 2), coords = coords)
  sf_obj <- as_sf(result)

  expect_s3_class(sf_obj, "sf")
  expect_equal(nrow(sf_obj), 20)
  expect_true("q0" %in% names(sf_obj))
  expect_true("q1" %in% names(sf_obj))
})


test_that("as_sf.spacc_profile errors without coords", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)
  result <- diversityProfile(species, q = c(0, 1))

  expect_error(as_sf(result), "No coordinates")
})


test_that("as_sf.spacc_profile errors without per_site", {
  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- diversityProfile(species, q = c(0, 1), type = "regional", coords = coords)
  expect_error(as_sf(result), "No per-site data")
})
