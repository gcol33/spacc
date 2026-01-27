test_that("groups parameter returns spacc_multi", {
  set.seed(1)
  species <- matrix(rbinom(50 * 20, 1, 0.4), nrow = 50)
  coords <- data.frame(x = runif(50), y = runif(50))
  groups <- rep(c("A", "B"), each = 10)

  result <- spacc(species, coords, n_seeds = 5, method = "knn",
                  groups = groups, seed = 42, progress = FALSE)

  expect_s3_class(result, "spacc_multi")
  expect_equal(length(result$objects), 2)
  expect_equal(result$names, c("A", "B"))
})


test_that("groups gives same result as manual split", {
  set.seed(1)
  species <- matrix(rbinom(50 * 20, 1, 0.4), nrow = 50)
  coords <- data.frame(x = runif(50), y = runif(50))
  groups <- rep(c("A", "B"), each = 10)

  # Grouped call
  grouped <- spacc(species, coords, n_seeds = 5, method = "knn",
                   groups = groups, seed = 42, progress = FALSE)

  # Manual split
  sac_a <- spacc(species[, groups == "A"], coords, n_seeds = 5, method = "knn",
                 seed = 42, progress = FALSE)
  sac_b <- spacc(species[, groups == "B"], coords, n_seeds = 5, method = "knn",
                 seed = 42, progress = FALSE)

  expect_equal(grouped$objects[[1]]$curves, sac_a$curves)
  expect_equal(grouped$objects[[2]]$curves, sac_b$curves)
})


test_that("groups validates length", {
  species <- matrix(1, nrow = 10, ncol = 5)
  coords <- data.frame(x = 1:10, y = 1:10)

  expect_error(
    spacc(species, coords, groups = c("A", "B"), progress = FALSE),
    "groups must have length equal to ncol"
  )
})


test_that("groups works with single group (identity)", {
  set.seed(1)
  species <- matrix(rbinom(50 * 10, 1, 0.4), nrow = 50)
  coords <- data.frame(x = runif(50), y = runif(50))
  groups <- rep("all", 10)

  result <- spacc(species, coords, n_seeds = 5, method = "knn",
                  groups = groups, seed = 42, progress = FALSE)

  direct <- spacc(species, coords, n_seeds = 5, method = "knn",
                  seed = 42, progress = FALSE)

  expect_s3_class(result, "spacc_multi")
  expect_equal(result$objects[[1]]$curves, direct$curves)
})
