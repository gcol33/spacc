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
                        focal_domain = NULL, domain_indices = NULL) {
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
  orders <- cpp_focal_order(coords$x, coords$y, focal_points$x, focal_points$y, distance) + 1L
  list(orders = orders, focal_points = focal_points, n_seeds = n_seeds)
}

.observed_site_orders <- function(coords, distance = "euclidean") {
  cpp_focal_order(coords$x, coords$y, coords$x, coords$y, distance) + 1L
}

.accumulation_orders <- function(method, coords, n_seeds, distance,
                                 dist_mat = NULL, focal_points = NULL,
                                 focal_domain = NULL, all_sites = FALSE) {
  if (method == "knn") {
    if (all_sites) {
      return(list(
        orders = .observed_site_orders(coords, distance),
        focal_points = coords[, c("x", "y"), drop = FALSE],
        n_seeds = nrow(coords)
      ))
    }
    return(.knn_orders(coords, n_seeds, distance, focal_points, focal_domain))
  }

  seeds <- if (all_sites) seq_len(nrow(coords)) - 1L else
    sample.int(nrow(coords), n_seeds, replace = TRUE) - 1L
  orders <- switch(method,
    nn_walk = {
      if (is.null(dist_mat)) {
        dist_mat <- cpp_distance_matrix(coords$x, coords$y, distance)
      }
      cpp_nn_walk_order(dist_mat, seeds) + 1L
    },
    kncn = cpp_kncn_order(coords$x, coords$y, seeds) + 1L,
    random = t(vapply(seq_len(length(seeds)), function(i) sample.int(nrow(coords)),
                      integer(nrow(coords)))),
    stop("Unknown accumulation method: ", method, call. = FALSE)
  )
  list(orders = orders, focal_points = NULL, n_seeds = nrow(orders))
}
