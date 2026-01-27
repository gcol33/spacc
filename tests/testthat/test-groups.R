test_that("groups parameter returns grouped spacc", {
  set.seed(1)
  species <- matrix(rbinom(50 * 20, 1, 0.4), nrow = 50)
  coords <- data.frame(x = runif(50), y = runif(50))
  groups <- rep(c("A", "B"), each = 10)

  result <- spacc(species, coords, n_seeds = 5, method = "knn",
                  groups = groups, seed = 42, progress = FALSE)

  expect_s3_class(result, "spacc")
  expect_true(spacc:::is_grouped(result))
  expect_equal(length(result$group_names), 2)
  expect_equal(result$group_names, c("A", "B"))
  expect_true(is.list(result$curves))
  expect_equal(names(result$curves), c("A", "B"))
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

  expect_equal(grouped$curves[["A"]], sac_a$curves)
  expect_equal(grouped$curves[["B"]], sac_b$curves)
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

  expect_s3_class(result, "spacc")
  expect_true(spacc:::is_grouped(result))
  expect_equal(result$curves[["all"]], direct$curves)
})


test_that("c.spacc creates grouped spacc", {
  set.seed(1)
  species <- matrix(rbinom(50 * 10, 1, 0.4), nrow = 50)
  coords <- data.frame(x = runif(50), y = runif(50))

  sac1 <- spacc(species, coords, n_seeds = 5, seed = 42, progress = FALSE)
  sac2 <- spacc(species, coords, n_seeds = 5, seed = 43, progress = FALSE)

  combined <- c(a = sac1, b = sac2)
  expect_s3_class(combined, "spacc")
  expect_true(spacc:::is_grouped(combined))
  expect_equal(combined$group_names, c("a", "b"))
})


test_that("summary works on grouped spacc", {
  set.seed(1)
  species <- matrix(rbinom(50 * 20, 1, 0.4), nrow = 50)
  coords <- data.frame(x = runif(50), y = runif(50))
  groups <- rep(c("A", "B"), each = 10)

  result <- spacc(species, coords, n_seeds = 5, method = "knn",
                  groups = groups, seed = 42, progress = FALSE)

  s <- summary(result)
  expect_true(is.list(s))
  expect_equal(names(s), c("A", "B"))
  expect_s3_class(s[["A"]], "summary.spacc")
})


test_that("as.data.frame works on grouped spacc", {
  set.seed(1)
  species <- matrix(rbinom(50 * 20, 1, 0.4), nrow = 50)
  coords <- data.frame(x = runif(50), y = runif(50))
  groups <- rep(c("A", "B"), each = 10)

  result <- spacc(species, coords, n_seeds = 5, method = "knn",
                  groups = groups, seed = 42, progress = FALSE)

  df <- as.data.frame(result)
  expect_true("group" %in% names(df))
  expect_equal(sort(unique(df$group)), c("A", "B"))
})
