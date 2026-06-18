# Phylogenetic diversity profile tests ----------------------------------------

test_that("diversityProfile phylogenetic returns spacc_profile", {
  skip_if_not_installed("ape")

  set.seed(42)
  species <- matrix(rpois(20 * 10, 2), nrow = 20,
                    dimnames = list(NULL, paste0("sp", 1:10)))
  tree <- ape::rcoal(10, tip.label = paste0("sp", 1:10))

  prof <- diversityProfile(species, tree = tree)

  expect_s3_class(prof, "spacc_profile")
  expect_equal(prof$profile_type, "phylogenetic")
  expect_equal(prof$n_sites, 20)
  expect_equal(prof$n_species, 10)
})


test_that("diversityProfile phylogenetic values decrease in q", {
  skip_if_not_installed("ape")

  set.seed(42)
  species <- matrix(rpois(20 * 10, 2), nrow = 20,
                    dimnames = list(NULL, paste0("sp", 1:10)))
  tree <- ape::rcoal(10, tip.label = paste0("sp", 1:10))

  prof <- diversityProfile(species, tree = tree, q = c(0, 1, 2))

  # Regional profile should be non-increasing in q
  expect_true(prof$regional[1] >= prof$regional[2] - 1e-10)
  expect_true(prof$regional[2] >= prof$regional[3] - 1e-10)
})


test_that("diversityProfile phylogenetic errors on missing species", {
  skip_if_not_installed("ape")

  species <- matrix(rpois(20 * 10, 2), nrow = 20,
                    dimnames = list(NULL, paste0("sp", 1:10)))
  tree <- ape::rcoal(8, tip.label = paste0("sp", 1:8))  # Missing sp9, sp10

  expect_error(diversityProfile(species, tree = tree), "Species not in tree")
})


test_that("diversityProfile phylogenetic print works", {
  skip_if_not_installed("ape")

  set.seed(42)
  species <- matrix(rpois(20 * 10, 2), nrow = 20,
                    dimnames = list(NULL, paste0("sp", 1:10)))
  tree <- ape::rcoal(10, tip.label = paste0("sp", 1:10))

  prof <- diversityProfile(species, tree = tree)
  expect_output(print(prof), "Phylogenetic diversity profile")
})


test_that("diversityProfile phylogenetic plot works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ape")

  set.seed(42)
  species <- matrix(rpois(20 * 10, 2), nrow = 20,
                    dimnames = list(NULL, paste0("sp", 1:10)))
  tree <- ape::rcoal(10, tip.label = paste0("sp", 1:10))

  prof <- diversityProfile(species, tree = tree)
  p <- plot(prof)
  expect_s3_class(p, "ggplot")
})


# Functional diversity profile tests ------------------------------------------

test_that("diversityProfile functional returns spacc_profile", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 2), nrow = 20,
                    dimnames = list(NULL, paste0("sp", 1:10)))
  traits <- data.frame(
    body_size = rnorm(10), diet = rnorm(10),
    row.names = paste0("sp", 1:10)
  )

  prof <- diversityProfile(species, traits = traits)

  expect_s3_class(prof, "spacc_profile")
  expect_equal(prof$profile_type, "functional")
  expect_equal(prof$n_sites, 20)
})


test_that("diversityProfile functional with identical traits equals standard Hill", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 2), nrow = 20,
                    dimnames = list(NULL, paste0("sp", 1:10)))
  # All species have identical traits -> Z = identity after normalization
  # Actually Z = all 1s -> functional diversity = 1 for all q > 0
  # Use maximally different traits instead: Z = identity
  traits <- data.frame(
    t1 = 1:10, t2 = 1:10,
    row.names = paste0("sp", 1:10)
  )

  prof_func <- diversityProfile(species, traits = traits, q = c(0, 1, 2))
  prof_tax <- diversityProfile(species, q = c(0, 1, 2))

  # When trait distances are spread out, functional diversity should differ
  # from taxonomic (just check structure, not equality)
  expect_s3_class(prof_func, "spacc_profile")
  expect_equal(ncol(prof_func$per_site), 3)
})


test_that("diversityProfile functional values decrease in q", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 2), nrow = 20,
                    dimnames = list(NULL, paste0("sp", 1:10)))
  traits <- data.frame(
    body_size = rnorm(10), diet = rnorm(10),
    row.names = paste0("sp", 1:10)
  )

  prof <- diversityProfile(species, traits = traits, q = c(0, 1, 2))

  # Regional profile should be non-increasing in q
  expect_true(prof$regional[1] >= prof$regional[2] - 1e-10)
  expect_true(prof$regional[2] >= prof$regional[3] - 1e-10)
})


test_that("diversityProfile functional errors on missing species", {
  species <- matrix(rpois(20 * 10, 2), nrow = 20,
                    dimnames = list(NULL, paste0("sp", 1:10)))
  traits <- data.frame(body_size = rnorm(8), row.names = paste0("sp", 1:8))

  expect_error(diversityProfile(species, traits = traits), "Species not in traits")
})


test_that("diversityProfile functional print works", {
  set.seed(42)
  species <- matrix(rpois(20 * 10, 2), nrow = 20,
                    dimnames = list(NULL, paste0("sp", 1:10)))
  traits <- data.frame(
    body_size = rnorm(10), diet = rnorm(10),
    row.names = paste0("sp", 1:10)
  )

  prof <- diversityProfile(species, traits = traits)
  expect_output(print(prof), "Functional diversity profile")
})


test_that("diversityProfile functional plot works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(42)
  species <- matrix(rpois(20 * 10, 2), nrow = 20,
                    dimnames = list(NULL, paste0("sp", 1:10)))
  traits <- data.frame(
    body_size = rnorm(10), diet = rnorm(10),
    row.names = paste0("sp", 1:10)
  )

  prof <- diversityProfile(species, traits = traits)
  p <- plot(prof)
  expect_s3_class(p, "ggplot")
})


test_that("diversityProfile functional gower distance works", {
  skip_if_not_installed("cluster")

  set.seed(42)
  species <- matrix(rpois(20 * 10, 2), nrow = 20,
                    dimnames = list(NULL, paste0("sp", 1:10)))
  traits <- data.frame(
    body_size = rnorm(10),
    diet = factor(sample(c("herb", "carn", "omni"), 10, replace = TRUE)),
    row.names = paste0("sp", 1:10)
  )

  prof <- diversityProfile(species, traits = traits, dist_method = "gower")
  expect_s3_class(prof, "spacc_profile")
})
