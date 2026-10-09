# Site-visitation orders.
#
# Every accumulation method is a rule for ordering sites; the species data are
# only read afterwards, along those orders. `.accumulation_orders()` is the one
# place that turns a method name into an n_orders x n_sites matrix of 1-based
# site indices, so every front door (spacc(), spaccHill(), spaccBeta(), ...)
# accumulates along identical traversals.

.sample_convex_hull <- function(points, n) {
  centered <- sweep(points, 2L, colMeans(points), "-")
  decomposition <- svd(centered, nu = 0L, nv = 2L)
  tolerance <- max(dim(centered)) * max(decomposition$d) * .Machine$double.eps
  rank <- sum(decomposition$d > tolerance)

  if (rank < 2L) {
    axis <- decomposition$v[, 1]
    projections <- drop(centered %*% axis)
    sampled <- stats::runif(n, min(projections), max(projections))
    return(outer(sampled, axis) + matrix(colMeans(points), n, 2L, byrow = TRUE))
  }

  hull <- points[grDevices::chull(points), , drop = FALSE]
  anchor <- hull[1, ]
  triangles <- lapply(seq.int(2L, nrow(hull) - 1L), function(i) {
    rbind(anchor, hull[i, ], hull[i + 1L, ])
  })
  areas <- vapply(triangles, function(triangle) {
    abs((triangle[2, 1] - triangle[1, 1]) * (triangle[3, 2] - triangle[1, 2]) -
        (triangle[3, 1] - triangle[1, 1]) * (triangle[2, 2] - triangle[1, 2])) / 2
  }, numeric(1))
  selected <- sample.int(length(triangles), n, replace = TRUE, prob = areas)
  first_weight <- sqrt(stats::runif(n))
  second_weight <- stats::runif(n)
  focus <- matrix(NA_real_, n, 2L)
  for (i in seq_len(n)) {
    triangle <- triangles[[selected[i]]]
    focus[i, ] <- (1 - first_weight[i]) * triangle[1, ] +
      first_weight[i] * (1 - second_weight[i]) * triangle[2, ] +
      first_weight[i] * second_weight[i] * triangle[3, ]
  }
  focus
}

.project_laea <- function(points) {
  radians <- points * pi / 180
  longitude <- radians[, 1]
  latitude <- radians[, 2]
  center_vector <- c(
    mean(cos(latitude) * cos(longitude)),
    mean(cos(latitude) * sin(longitude)),
    mean(sin(latitude))
  )
  center_longitude <- atan2(center_vector[2], center_vector[1])
  center_latitude <- atan2(center_vector[3], sqrt(sum(center_vector[1:2]^2)))
  delta_longitude <- atan2(sin(longitude - center_longitude),
                           cos(longitude - center_longitude))
  denominator <- 1 + sin(center_latitude) * sin(latitude) +
    cos(center_latitude) * cos(latitude) * cos(delta_longitude)
  if (any(denominator <= .Machine$double.eps)) {
    stop("Geographic coordinates span too large a region for focal-point sampling.",
         call. = FALSE)
  }
  scale <- sqrt(2 / denominator)
  projected <- cbind(
    scale * cos(latitude) * sin(delta_longitude),
    scale * (cos(center_latitude) * sin(latitude) -
      sin(center_latitude) * cos(latitude) * cos(delta_longitude))
  )
  list(points = projected, longitude = center_longitude, latitude = center_latitude)
}

.inverse_laea <- function(points, center_longitude, center_latitude) {
  radius <- sqrt(rowSums(points^2))
  angular_distance <- 2 * asin(pmin(radius, 2) / 2)
  latitude <- rep(center_latitude, nrow(points))
  longitude <- rep(center_longitude, nrow(points))
  away_from_center <- radius > sqrt(.Machine$double.eps)
  latitude[away_from_center] <- asin(
    cos(angular_distance[away_from_center]) * sin(center_latitude) +
      points[away_from_center, 2] * sin(angular_distance[away_from_center]) *
      cos(center_latitude) / radius[away_from_center]
  )
  longitude[away_from_center] <- center_longitude + atan2(
    points[away_from_center, 1] * sin(angular_distance[away_from_center]),
    radius[away_from_center] * cos(center_latitude) *
      cos(angular_distance[away_from_center]) -
      points[away_from_center, 2] * sin(center_latitude) *
      sin(angular_distance[away_from_center])
  )
  cbind(
    ((longitude + pi) %% (2 * pi) - pi) * 180 / pi,
    latitude * 180 / pi
  )
}

.sample_focal_points <- function(coords, n, domain = NULL,
                                 distance = "euclidean") {
  if (!is.null(domain)) {
    check_suggests("sf")
    if (!inherits(domain, c("sf", "sfc"))) {
      stop("`focal_domain` must be an sf or sfc polygon.", call. = FALSE)
    }
    geometry <- sf::st_geometry(domain)
    if (!all(as.character(sf::st_geometry_type(geometry)) %in% c("POLYGON", "MULTIPOLYGON"))) {
      stop("`focal_domain` must contain polygon geometries.", call. = FALSE)
    }
    sampled <- sf::st_sample(sf::st_union(geometry), size = n, type = "random", exact = TRUE)
    sampled_coords <- sf::st_coordinates(sampled)
    return(data.frame(x = sampled_coords[, 1], y = sampled_coords[, 2]))
  }

  points <- unique(as.matrix(coords[, c("x", "y"), drop = FALSE]))
  storage.mode(points) <- "double"
  if (any(!is.finite(points))) {
    stop("Coordinates must be finite to sample focal points.", call. = FALSE)
  }
  if (nrow(points) == 1L) {
    return(data.frame(x = rep(points[1, 1], n), y = rep(points[1, 2], n)))
  }

  if (distance == "haversine") {
    projection <- .project_laea(points)
    sampled <- .sample_convex_hull(projection$points, n)
    focus <- .inverse_laea(sampled, projection$longitude, projection$latitude)
  } else {
    focus <- .sample_convex_hull(points, n)
  }
  data.frame(x = focus[, 1], y = focus[, 2])
}

.validate_focal_points <- function(focal_points) {
  if (inherits(focal_points, "sf")) {
    coordinates <- sf::st_coordinates(focal_points)
    focal_points <- data.frame(x = coordinates[, 1], y = coordinates[, 2])
  }
  if (!is.data.frame(focal_points) || !all(c("x", "y") %in% names(focal_points))) {
    stop("`focal_points` must have x and y columns.", call. = FALSE)
  }
  focal_points <- focal_points[, c("x", "y"), drop = FALSE]
  if (!is.numeric(focal_points$x) || !is.numeric(focal_points$y) ||
      any(!is.finite(as.matrix(focal_points)))) {
    stop("`focal_points` coordinates must be finite numeric values.", call. = FALSE)
  }
  focal_points
}

.knn_orders <- function(coords, n_seeds, distance, focal_points = NULL,
                        focal_domain = NULL, domain_indices = NULL, n_cores = 1L) {
  if (!is.null(focal_points) && !is.null(focal_domain)) {
    stop("Supply either `focal_points` or `focal_domain`, not both.", call. = FALSE)
  }
  if (is.null(focal_points)) {
    domain_coords <- if (is.null(domain_indices)) coords else coords[domain_indices, , drop = FALSE]
    focal_points <- .sample_focal_points(domain_coords, n_seeds, focal_domain,
                                         distance)
  } else {
    focal_points <- .validate_focal_points(focal_points)
    n_seeds <- nrow(focal_points)
  }
  orders <- cpp_focal_order(coords$x, coords$y, focal_points$x, focal_points$y,
                            distance, n_cores) + 1L
  list(orders = orders, focal_points = focal_points, n_seeds = n_seeds)
}

# Fixed-focus orders with every observed site as its own focus (one row per site).
.observed_site_orders <- function(coords, distance = "euclidean", n_cores = 1L) {
  cpp_focal_order(coords$x, coords$y, coords$x, coords$y, distance, n_cores) + 1L
}

# Should the nearest-neighbour methods query a spatial tree instead of a
# distance matrix? A supplied matrix (spacc_dist input, spatiotemporal
# composite) is always used as given.
.use_tree <- function(backend, n_sites, dist_mat) {
  if (!is.null(dist_mat)) return(FALSE)
  switch(backend,
    kdtree = TRUE,
    exact = FALSE,
    auto = n_sites > 500L,
    stop("Unknown backend: ", backend, call. = FALSE)
  )
}

#' Site-visitation orders for an accumulation method
#'
#' @param method One of `"knn"`, `"kncn"`, `"nn_walk"`, `"random"`,
#'   `"gaussian"`, `"cone"`, `"collector"`.
#' @param coords Data frame with `x` and `y` columns.
#' @param n_seeds Number of orders to draw (ignored when `all_sites = TRUE`).
#' @param distance `"euclidean"` or `"haversine"`.
#' @param dist_mat Optional precomputed site distance matrix. Methods that need
#'   one compute it when `NULL`.
#' @param focal_points,focal_domain Passed to [.knn_orders()] for `"knn"`.
#' @param all_sites Use every site once as seed (or focus) instead of sampling.
#' @param seed_pool Optional 1-based site indices from which seeds are drawn.
#' @param backend `"auto"`, `"exact"` or `"kdtree"` for `"nn_walk"` and `"kncn"`.
#' @param sigma Gaussian bandwidth; `NULL` uses the median positive distance.
#' @param cone_width Cone half-width in radians.
#' @param n_cores Threads for the parallel order generators.
#' @return List with `orders` (1-based matrix, one row per order), `n_seeds`,
#'   `focal_points`, `backend` label and resolved `sigma`.
#' @noRd
.accumulation_orders <- function(method, coords, n_seeds, distance = "euclidean",
                                 dist_mat = NULL, focal_points = NULL,
                                 focal_domain = NULL, all_sites = FALSE,
                                 seed_pool = NULL, backend = "auto",
                                 sigma = NULL, cone_width = pi / 4,
                                 n_cores = 1L) {
  n_sites <- nrow(coords)
  result <- list(focal_points = NULL, backend = NA_character_, sigma = NULL)

  if (method == "knn") {
    knn <- if (all_sites) {
      list(orders = .observed_site_orders(coords, distance, n_cores),
           focal_points = coords[, c("x", "y"), drop = FALSE],
           n_seeds = n_sites)
    } else {
      .knn_orders(coords, n_seeds, distance, focal_points, focal_domain,
                  domain_indices = seed_pool, n_cores = n_cores)
    }
    result[names(knn)] <- knn
    result$backend <- "fixed_focus"
    return(result)
  }

  if (method == "collector") {
    result$orders <- matrix(seq_len(n_sites), nrow = 1L)
    result$n_seeds <- 1L
    return(result)
  }

  seeds <- if (all_sites) {
    seq_len(n_sites)
  } else {
    pool <- if (is.null(seed_pool)) seq_len(n_sites) else seed_pool
    pool[sample.int(length(pool), n_seeds, replace = TRUE)]
  }
  seeds0 <- seeds - 1L
  n_orders <- length(seeds)

  needs_matrix <- method == "gaussian" ||
    (method == "nn_walk" && !.use_tree(backend, n_sites, dist_mat))
  if (needs_matrix && is.null(dist_mat)) {
    dist_mat <- cpp_distance_matrix(coords$x, coords$y, distance)
  }

  orders <- switch(method,
    nn_walk = if (.use_tree(backend, n_sites, dist_mat)) {
      result$backend <- "kdtree"
      cpp_nn_walk_tree_order(coords$x, coords$y, seeds0, distance, n_cores)
    } else {
      result$backend <- "exact"
      cpp_nn_walk_order(dist_mat, seeds0, n_cores)
    },
    kncn = if (.use_tree(backend, n_sites, dist_mat)) {
      result$backend <- "kdtree"
      cpp_kncn_tree_order(coords$x, coords$y, seeds0, distance, n_cores)
    } else {
      result$backend <- "exact"
      cpp_kncn_order(coords$x, coords$y, seeds0, distance, n_cores)
    },
    random = t(vapply(seq_len(n_orders), function(i) sample.int(n_sites) - 1L,
                      integer(n_sites))),
    gaussian = {
      if (is.null(sigma)) sigma <- stats::median(dist_mat[dist_mat > 0])
      result$sigma <- sigma
      cpp_gaussian_order(dist_mat, seeds0, sigma)
    },
    cone = cpp_cone_order(coords$x, coords$y, seeds0,
                          stats::runif(n_orders, 0, 2 * pi), cone_width),
    stop("Unknown accumulation method: ", method, call. = FALSE)
  )

  result$orders <- orders + 1L
  result$n_seeds <- n_orders
  result
}

# Last value of every row: the end point of each curve.
.final_per_seed <- function(mat) {
  mat[cbind(seq_len(nrow(mat)), ncol(mat))]
}

# Per-site map values: run the method from every site once and keep the final
# value of each curve. `compute(orders)` returns a named list of per-order
# vectors (one element per mapped column).
.site_map_values <- function(method, coords, distance, dist_mat, n_cores, compute) {
  map_orders <- .accumulation_orders(method, coords, nrow(coords), distance,
                                     dist_mat, all_sites = TRUE, n_cores = n_cores)
  values <- compute(map_orders$orders)
  cbind(
    data.frame(site_id = seq_len(nrow(coords)), x = coords$x, y = coords$y),
    as.data.frame(values, optional = TRUE)
  )
}
