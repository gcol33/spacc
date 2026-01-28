test_that("support = 'auto' splits by country", {
  skip_if_not_installed("sf")
  skip_if_not_installed("areaOfEffect")
  set.seed(1)
  n <- 40

  # Sites in France (~47N, 2E) and Germany (~51N, 10E)
  coords <- data.frame(
    x = c(runif(20, 0, 4), runif(20, 8, 12)),
    y = c(runif(20, 44, 48), runif(20, 48, 52))
  )

  # Skip if areaOfEffect country detection doesn't work (e.g. missing internal data)
  skip_if_not(
    tryCatch({
      pts <- sf::st_as_sf(coords, coords = c("x", "y"), crs = 4326)
      areaOfEffect::aoe(pts)
      TRUE
    }, error = function(e) FALSE),
    "areaOfEffect country detection not functional"
  )

  species <- matrix(rbinom(n * 15, 1, 0.4), nrow = n)

  result <- spacc(species, coords, n_seeds = 3, method = "knn",
                  support = "auto", distance = "haversine",
                  seed = 42, progress = FALSE)

  expect_s3_class(result, "spacc")
  expect_true(spacc:::is_grouped(result))
  expect_true(length(result$group_names) >= 2)
  expect_true(is.list(result$curves))
  expect_true(!is.null(result$support$auto))
  expect_true(result$support$auto)
})

test_that("support = 'auto' with groups produces compound names", {
  skip_if_not_installed("sf")
  skip_if_not_installed("areaOfEffect")

  set.seed(1)
  n <- 40

  coords <- data.frame(
    x = c(runif(20, 0, 4), runif(20, 8, 12)),
    y = c(runif(20, 44, 48), runif(20, 48, 52))
  )

  skip_if_not(
    tryCatch({
      pts <- sf::st_as_sf(coords, coords = c("x", "y"), crs = 4326)
      areaOfEffect::aoe(pts)
      TRUE
    }, error = function(e) FALSE),
    "areaOfEffect country detection not functional"
  )

  species <- matrix(rbinom(n * 10, 1, 0.4), nrow = n)
  groups <- rep(c("birds", "mammals"), each = 5)

  result <- spacc(species, coords, n_seeds = 3, method = "knn",
                  support = "auto", distance = "haversine",
                  groups = groups, seed = 42, progress = FALSE)

  expect_s3_class(result, "spacc")
  expect_true(spacc:::is_grouped(result))
  # Should have compound names like "France.birds", "France.mammals", etc.
  expect_true(all(grepl("\\.", result$group_names)))
})
