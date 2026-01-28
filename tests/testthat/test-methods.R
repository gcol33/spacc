test_that("print.spacc works for ungrouped", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  sac <- spacc(species, coords, n_seeds = 3, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)

  expect_output(print(sac), "spacc")
})


test_that("summary.spacc returns correct structure", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rbinom(20 * 10, 1, 0.4), nrow = 20)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)

  summ <- summary(sac)

  expect_s3_class(summ, "summary.spacc")
  expect_equal(length(summ$sites), 20)
  expect_equal(length(summ$mean), 20)
  expect_equal(length(summ$lower), 20)
  expect_equal(length(summ$upper), 20)
  expect_true(!is.na(summ$saturation_point))
  expect_true(summ$saturation_point >= 1 && summ$saturation_point <= 20)
})


test_that("print.summary.spacc works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  sac <- spacc(species, coords, n_seeds = 3, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)

  expect_output(print(summary(sac)), "Spatial Species Accumulation")
})


test_that("as.data.frame.spacc works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  sac <- spacc(species, coords, n_seeds = 3, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)

  df <- as.data.frame(sac)

  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 15)
  expect_true(all(c("sites", "mean", "lower", "upper", "sd") %in% names(df)))
})


test_that("[.spacc subsets seeds", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  species <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)

  sac <- spacc(species, coords, n_seeds = 5, method = "knn",
               parallel = FALSE, progress = FALSE, seed = 1)

  sub <- sac[1:3]

  expect_s3_class(sub, "spacc")
  expect_equal(sub$n_seeds, 3)
  expect_equal(nrow(sub$curves), 3)
})


test_that("c.spacc combines objects", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  sp1 <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)
  sp2 <- matrix(rbinom(15 * 8, 1, 0.2), nrow = 15)

  sac1 <- spacc(sp1, coords, n_seeds = 3, method = "knn",
                parallel = FALSE, progress = FALSE, seed = 1)
  sac2 <- spacc(sp2, coords, n_seeds = 3, method = "knn",
                parallel = FALSE, progress = FALSE, seed = 1)

  combined <- c(native = sac1, alien = sac2)

  expect_s3_class(combined, "spacc")
  expect_true(spacc:::is_grouped(combined))
  expect_equal(combined$group_names, c("native", "alien"))
  expect_equal(length(combined$curves), 2)
})


test_that("is_grouped detects grouped objects", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  sp1 <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)
  sp2 <- matrix(rbinom(15 * 8, 1, 0.2), nrow = 15)

  sac1 <- spacc(sp1, coords, n_seeds = 3, method = "knn",
                parallel = FALSE, progress = FALSE, seed = 1)
  sac2 <- spacc(sp2, coords, n_seeds = 3, method = "knn",
                parallel = FALSE, progress = FALSE, seed = 1)

  expect_false(spacc:::is_grouped(sac1))

  combined <- c(a = sac1, b = sac2)
  expect_true(spacc:::is_grouped(combined))
})


test_that("print.spacc works for grouped", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  sp1 <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)
  sp2 <- matrix(rbinom(15 * 8, 1, 0.2), nrow = 15)

  sac1 <- spacc(sp1, coords, n_seeds = 3, method = "knn",
                parallel = FALSE, progress = FALSE, seed = 1)
  sac2 <- spacc(sp2, coords, n_seeds = 3, method = "knn",
                parallel = FALSE, progress = FALSE, seed = 1)

  combined <- c(native = sac1, alien = sac2)

  expect_output(print(combined), "groups")
})


test_that("summary.spacc works for grouped", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  sp1 <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)
  sp2 <- matrix(rbinom(15 * 8, 1, 0.2), nrow = 15)

  sac1 <- spacc(sp1, coords, n_seeds = 3, method = "knn",
                parallel = FALSE, progress = FALSE, seed = 1)
  sac2 <- spacc(sp2, coords, n_seeds = 3, method = "knn",
                parallel = FALSE, progress = FALSE, seed = 1)

  combined <- c(native = sac1, alien = sac2)
  summ <- summary(combined)

  expect_true(is.list(summ))
  expect_equal(names(summ), c("native", "alien"))
})


test_that("as.data.frame.spacc works for grouped", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(15), y = runif(15))
  sp1 <- matrix(rbinom(15 * 8, 1, 0.4), nrow = 15)
  sp2 <- matrix(rbinom(15 * 8, 1, 0.2), nrow = 15)

  sac1 <- spacc(sp1, coords, n_seeds = 3, method = "knn",
                parallel = FALSE, progress = FALSE, seed = 1)
  sac2 <- spacc(sp2, coords, n_seeds = 3, method = "knn",
                parallel = FALSE, progress = FALSE, seed = 1)

  combined <- c(native = sac1, alien = sac2)
  df <- as.data.frame(combined)

  expect_s3_class(df, "data.frame")
  expect_true("group" %in% names(df))
  expect_equal(nrow(df), 30)  # 15 sites x 2 groups
})
