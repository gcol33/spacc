test_that("compare returns correct structure", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  sp1 <- matrix(rbinom(20 * 10, 1, 0.4), nrow = 20)
  sp2 <- matrix(rbinom(20 * 10, 1, 0.2), nrow = 20)

  sac1 <- spacc(sp1, coords, n_seeds = 5, method = "knn",
                parallel = FALSE, progress = FALSE, seed = 1)
  sac2 <- spacc(sp2, coords, n_seeds = 5, method = "knn",
                parallel = FALSE, progress = FALSE, seed = 1)

  comp <- compare(sac1, sac2, n_perm = 99)

  expect_s3_class(comp, "spacc_comp")
  expect_true(!is.na(comp$auc_diff))
  expect_true(!is.na(comp$p_value))
  expect_true(comp$p_value >= 0 && comp$p_value <= 1)
  expect_equal(comp$method, "permutation")
})


test_that("compare works with bootstrap method", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  sp1 <- matrix(rbinom(20 * 10, 1, 0.4), nrow = 20)
  sp2 <- matrix(rbinom(20 * 10, 1, 0.2), nrow = 20)

  sac1 <- spacc(sp1, coords, n_seeds = 5, method = "knn",
                parallel = FALSE, progress = FALSE, seed = 1)
  sac2 <- spacc(sp2, coords, n_seeds = 5, method = "knn",
                parallel = FALSE, progress = FALSE, seed = 1)

  comp <- compare(sac1, sac2, method = "bootstrap", n_perm = 99)

  expect_s3_class(comp, "spacc_comp")
  expect_equal(comp$method, "bootstrap")
})


test_that("compare works with auc method", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  sp1 <- matrix(rbinom(20 * 10, 1, 0.4), nrow = 20)
  sp2 <- matrix(rbinom(20 * 10, 1, 0.2), nrow = 20)

  sac1 <- spacc(sp1, coords, n_seeds = 5, method = "knn",
                parallel = FALSE, progress = FALSE, seed = 1)
  sac2 <- spacc(sp2, coords, n_seeds = 5, method = "knn",
                parallel = FALSE, progress = FALSE, seed = 1)

  comp <- compare(sac1, sac2, method = "auc")

  expect_s3_class(comp, "spacc_comp")
  expect_true(is.na(comp$p_value))
  expect_true(!is.na(comp$auc_diff))
})


test_that("compare errors on mismatched sites", {
  skip_on_cran()

  set.seed(42)
  coords1 <- data.frame(x = runif(20), y = runif(20))
  coords2 <- data.frame(x = runif(15), y = runif(15))
  sp1 <- matrix(rbinom(20 * 10, 1, 0.4), nrow = 20)
  sp2 <- matrix(rbinom(15 * 10, 1, 0.4), nrow = 15)

  sac1 <- spacc(sp1, coords1, n_seeds = 3, method = "knn",
                parallel = FALSE, progress = FALSE, seed = 1)
  sac2 <- spacc(sp2, coords2, n_seeds = 3, method = "knn",
                parallel = FALSE, progress = FALSE, seed = 1)

  expect_error(compare(sac1, sac2), "same number of sites")
})


test_that("print.spacc_comp works", {
  skip_on_cran()

  set.seed(42)
  coords <- data.frame(x = runif(20), y = runif(20))
  sp1 <- matrix(rbinom(20 * 10, 1, 0.4), nrow = 20)
  sp2 <- matrix(rbinom(20 * 10, 1, 0.2), nrow = 20)

  sac1 <- spacc(sp1, coords, n_seeds = 5, method = "knn",
                parallel = FALSE, progress = FALSE, seed = 1)
  sac2 <- spacc(sp2, coords, n_seeds = 5, method = "knn",
                parallel = FALSE, progress = FALSE, seed = 1)

  comp <- compare(sac1, sac2, n_perm = 99)

  expect_output(print(comp), "Comparison")
})
