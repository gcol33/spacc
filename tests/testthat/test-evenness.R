test_that("evenness hill returns correct structure", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- evenness(species, q = seq(0.5, 2, by = 0.5))

  expect_s3_class(result, "spacc_evenness")
  expect_equal(result$type, "hill")
  expect_equal(nrow(result$per_site), 20)
  expect_equal(ncol(result$per_site), 4)
  expect_true(!is.null(result$regional))
})


test_that("hill evenness bounded [0, 1]", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- evenness(species, q = seq(0.5, 3, by = 0.5))

  expect_true(all(result$per_site >= -1e-10))
  expect_true(all(result$per_site <= 1 + 1e-10))
  expect_true(all(result$regional >= -1e-10))
  expect_true(all(result$regional <= 1 + 1e-10))
})


test_that("perfectly even community has evenness = 1", {
  # All species have equal abundance
  species <- matrix(5, nrow = 10, ncol = 10)

  result <- evenness(species, q = c(0.5, 1, 2))

  # All sites should have evenness = 1
  for (j in seq_len(ncol(result$per_site))) {
    expect_equal(result$per_site[, j], rep(1, 10), tolerance = 1e-8)
  }
})


test_that("S = 1 site returns evenness = 1", {
  # Only one species present at each site
  species <- matrix(0, nrow = 5, ncol = 10)
  for (i in 1:5) species[i, i] <- 10

  result <- evenness(species, q = c(0.5, 1, 2))

  expect_true(all(result$per_site == 1))
})


test_that("pielou evenness returns correct structure", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- evenness(species, type = "pielou")

  expect_s3_class(result, "spacc_evenness")
  expect_equal(result$type, "pielou")
  expect_true(is.numeric(result$per_site))
  expect_length(result$per_site, 20)
})


test_that("pielou bounded [0, 1]", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- evenness(species, type = "pielou")

  expect_true(all(result$per_site >= -1e-10))
  expect_true(all(result$per_site <= 1 + 1e-10))
})


test_that("simpson evenness returns correct structure", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- evenness(species, type = "simpson")

  expect_s3_class(result, "spacc_evenness")
  expect_equal(result$type, "simpson")
  expect_true(is.numeric(result$per_site))
  expect_length(result$per_site, 20)
})


test_that("simpson bounded [0, 1]", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- evenness(species, type = "simpson")

  expect_true(all(result$per_site >= -1e-10))
  expect_true(all(result$per_site <= 1 + 1e-10))
})


test_that("print.spacc_evenness works for all types", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  expect_output(print(evenness(species, q = c(0.5, 1, 2))), "Evenness")
  expect_output(print(evenness(species, type = "pielou")), "Pielou")
  expect_output(print(evenness(species, type = "simpson")), "Simpson")
})


test_that("summary.spacc_evenness returns data.frame", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- evenness(species, q = seq(0.5, 2, by = 0.5))
  summ <- summary(result)

  expect_s3_class(summ, "data.frame")
  expect_true("mean" %in% names(summ))
  expect_true("regional" %in% names(summ))
})


test_that("as.data.frame.spacc_evenness works", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- evenness(species, q = seq(0.5, 2, by = 0.5))
  df <- as.data.frame(result)

  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 4)  # 4 q values
})


test_that("plot.spacc_evenness profile works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)
  result <- evenness(species, q = seq(0.5, 2, by = 0.5))

  p <- plot(result)
  expect_s3_class(p, "ggplot")
})


test_that("plot.spacc_evenness histogram works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)
  result <- evenness(species, type = "pielou")

  p <- plot(result, type = "histogram")
  expect_s3_class(p, "ggplot")
})


test_that("as_sf.spacc_evenness works for hill type", {
  skip_if_not_installed("sf")

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- evenness(species, q = c(0.5, 1, 2), coords = coords)
  sf_obj <- as_sf(result)

  expect_s3_class(sf_obj, "sf")
  expect_equal(nrow(sf_obj), 20)
  expect_true("E_q0.5" %in% names(sf_obj))
  expect_true("E_q1" %in% names(sf_obj))
})


test_that("as_sf.spacc_evenness works for pielou type", {
  skip_if_not_installed("sf")

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 3), nrow = 20)

  result <- evenness(species, type = "pielou", coords = coords)
  sf_obj <- as_sf(result)

  expect_s3_class(sf_obj, "sf")
  expect_true("pielou_evenness" %in% names(sf_obj))
})


test_that("as_sf.spacc_evenness errors without coords", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 3), nrow = 20)
  result <- evenness(species)

  expect_error(as_sf(result), "No coordinates")
})
