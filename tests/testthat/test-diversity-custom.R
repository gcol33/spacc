test_that("spaccDiversity reproduces collector richness when fun counts species", {
  set.seed(21)
  coords <- data.frame(x = runif(25), y = runif(25))
  species <- matrix(rbinom(25 * 12, 1, 0.3), nrow = 25)

  custom <- spaccDiversity(species, coords, function(comm) sum(comm > 0),
                           method = "collector", progress = FALSE)
  builtin <- spacc(species, coords, method = "collector", progress = FALSE)

  expect_s3_class(custom, "spacc_diversity")
  expect_s3_class(custom, "spacc")
  expect_equal(as.numeric(custom$curves[1, ]), as.numeric(builtin$curves[1, ]))
})


test_that("spaccDiversity knn richness is monotone and ends at total richness", {
  set.seed(22)
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)

  res <- spaccDiversity(species, coords, function(comm) sum(comm > 0),
                        method = "knn", n_seeds = 6, progress = FALSE)
  expect_equal(dim(res$curves), c(6L, 30L))
  expect_true(all(res$curves[, 30] == sum(colSums(species) > 0)))
  expect_true(all(apply(res$curves, 1, function(r) all(diff(r) >= 0))))
})


test_that("spaccDiversity inherits spacc summary and data.frame methods", {
  set.seed(23)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 2), nrow = 20)

  res <- spaccDiversity(species, coords, function(comm) sum(comm > 0),
                        method = "random", n_seeds = 8, progress = FALSE)
  s <- summary(res)
  expect_s3_class(s, "summary.spacc")
  df <- as.data.frame(res)
  expect_true(all(c("sites", "mean", "lower", "upper") %in% names(df)))
})


test_that("spaccDiversity passes extra args and incidence flag", {
  set.seed(25)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 2), nrow = 20)

  # incidence = TRUE feeds 0/1; counting non-zero gives richness
  res <- spaccDiversity(species, coords, function(comm, thr) sum(comm >= thr),
                        thr = 1, method = "collector", incidence = TRUE,
                        progress = FALSE)
  expect_equal(res$curves[1, 20], sum(colSums(species) > 0))
})


test_that("plot.spacc_diversity returns a ggplot with a metric-neutral label", {
  skip_if_not_installed("ggplot2")
  set.seed(26)
  coords <- data.frame(x = runif(20), y = runif(20))
  species <- matrix(rpois(20 * 10, 2), nrow = 20)

  res <- spaccDiversity(species, coords, function(comm) sum(comm > 0),
                        method = "knn", n_seeds = 5, progress = FALSE)
  p <- plot(res)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$y, "Cumulative diversity")
  expect_equal(p$labels$title, "Custom Diversity Accumulation")

  p2 <- plot(res, ylab = "Shannon")
  expect_equal(p2$labels$y, "Shannon")
})


test_that("spaccDiversity rejects non-scalar fun output", {
  coords <- data.frame(x = runif(10), y = runif(10))
  species <- matrix(rbinom(10 * 5, 1, 0.4), nrow = 10)
  expect_error(
    spaccDiversity(species, coords, function(comm) comm,
                   method = "collector", progress = FALSE),
    "single numeric value"
  )
})


test_that("spaccDiversity works across orderings and ends at total richness", {
  set.seed(24)
  coords <- data.frame(x = runif(18), y = runif(18))
  species <- matrix(rbinom(18 * 8, 1, 0.4), nrow = 18)
  total <- sum(colSums(species) > 0)
  f <- function(comm) sum(comm > 0)

  for (m in c("kncn", "radius", "collector")) {
    res <- spaccDiversity(species, coords, f, method = m, n_seeds = 4, progress = FALSE)
    expect_true(all(res$curves[, 18] == total), info = m)
  }
})
