#' Spatial Species Accumulation Curves
#'
#' Compute species accumulation curves using various spatial sampling methods
#' with C++ backend for performance.
#'
#' @param x A site-by-species matrix (rows = sites, cols = species) with
#'   presence/absence (0/1) or abundance data. Can also be a data.frame.
#' @param coords A data.frame with columns `x` and `y` containing site coordinates,
#'   or a `spacc_dist` object from [distances()].
#' @param n_seeds Integer. Number of random starting points for uncertainty
#'   quantification. Default 50.
#' @param method Character. Accumulation method:
#'   - `"knn"`: k-Nearest Neighbor (always visit closest unvisited)
#'   - `"kncn"`: k-Nearest Centroid Neighbor (visit closest to centroid)
#'   - `"random"`: Random order (null model)
#'   - `"radius"`: Expand by distance from seed
#'   - `"gaussian"`: Probabilistic selection weighted by distance
#'   - `"cone"`: Directional expansion within angular constraint
#'   - `"collector"`: Sites in data order (no randomization, single curve)
#' @param distance Character. Distance method: `"euclidean"` or `"haversine"`.
#' @param support Optional. Spatial support for core/halo classification via
#'   [areaOfEffect::aoe()]. Can be:
#'   - `"auto"`: Auto-detect countries and run per-country accumulation,
#'     returning a grouped `spacc` object with one curve per country
#'   - Country name or ISO code: `"France"`, `"FR"`, `"FRA"`
#'   - Vector of countries: `c("France", "Germany")`
#'   - An `sf` polygon object
#'   - An `aoe_result` object (pre-computed)
#'   When provided, seeds are sampled only from "core" sites (inside support),
#'   while accumulation can expand into "halo" sites (buffer zone).
#' @param include_halo Logical. When `support` is provided, should halo sites
#'   be included in accumulation? Default `TRUE` (ecological boundary).
#'   Set to `FALSE` for political/hard boundary.
#' @param backend Character. Nearest-neighbor backend for `knn` and `kncn`:
#'   - `"auto"` (default): Uses exact (brute-force) for ≤500 sites,
#'     spatial tree for >500 sites.
#'   - `"exact"`: Always use brute-force with precomputed distance matrix.
#'   - `"kdtree"`: Always use spatial tree. Uses k-d tree (nanoflann) for
#'     Euclidean distances and ball tree for haversine distances. Faster for
#'     large datasets, no distance matrix needed.
#' @param sigma Numeric. Bandwidth for Gaussian method. Default auto-calculated.
#' @param cone_width Numeric. Half-width in radians for cone method. Default pi/4.
#' @param parallel Logical. Use parallel processing? Default `TRUE`.
#' @param n_cores Integer. Number of cores. Default `NULL` uses `detectCores() - 1`.
#' @param progress Logical. Show progress bar? Default `TRUE`.
#' @param groups Optional. A factor, character, or integer vector of length
#'   `ncol(x)` assigning each species (column) to a group. When provided,
#'   separate accumulation curves are computed for each group using the
#'   **same spatial site ordering**, and a grouped `spacc` object is returned.
#'   Useful for comparing native vs alien species, families, or any
#'   categorical split. Default `NULL` (no grouping).
#' @param time Optional. Numeric vector of length `nrow(x)` giving a temporal
#'   coordinate for each site. When provided, a combined spatiotemporal distance
#'   matrix is computed as `w_space * d_spatial + w_time * d_temporal` and used
#'   for accumulation. Forces exact (brute-force) backend since spatial trees
#'   cannot handle composite distances. Only supported for methods that use a
#'   distance matrix: `"knn"`, `"radius"`, `"gaussian"`.
#' @param w_space Numeric. Weight for spatial distance when `time` is provided.
#'   Default 1.
#' @param w_time Numeric. Weight for temporal distance when `time` is provided.
#'   Default 1.
#' @param seed Integer. Random seed for reproducibility. Default `NULL`.
#'
#' @return When `groups = NULL`, an object of class `spacc` containing:
#'   \item{curves}{Matrix of cumulative species counts (n_seeds x n_sites)}
#'   \item{coords}{Original coordinates}
#'   \item{n_seeds}{Number of seeds used}
#'   \item{method}{Method used}
#'   \item{n_species}{Total species in dataset}
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' sac <- spacc(species, coords)
#' plot(sac)
#'
#' # Different methods
#' sac_knn <- spacc(species, coords, method = "knn")
#' sac_gauss <- spacc(species, coords, method = "gaussian", sigma = 10)
#' sac_cone <- spacc(species, coords, method = "cone", cone_width = pi/6)
#'
#' # Compare to null model
#' sac_rand <- spacc(species, coords, method = "random")
#' comp <- compare(sac_knn, sac_rand)
#'
#' # With spatial support (seeds from France, accumulate into neighbors)
#' sac_france <- spacc(species, coords, support = "France")
#'
#' # Hard boundary (France only, no halo)
#' sac_france_only <- spacc(species, coords, support = "France", include_halo = FALSE)
#'
#' # Grouped accumulation (e.g., native vs alien)
#' status <- ifelse(grepl("alien", colnames(species)), "alien", "native")
#' sac_grouped <- spacc(species, coords, groups = status, seed = 42)
#' plot(sac_grouped)  # Overlaid curves per group
#' }
#'
#' @export
spacc <- function(x,
                  coords,
                  n_seeds = 50L,
                  method = c("knn", "kncn", "random", "radius", "gaussian", "cone", "collector"),
                  distance = c("euclidean", "haversine"),
                  backend = c("auto", "exact", "kdtree"),
                  support = NULL,
                  include_halo = TRUE,
                  sigma = NULL,
                  cone_width = pi / 4,
                  parallel = TRUE,
                  n_cores = NULL,
                  progress = TRUE,
                  groups = NULL,
                  time = NULL,
                  w_space = 1,
                  w_time = 1,
                  seed = NULL) {

  method <- match.arg(method)
  distance <- match.arg(distance)
  backend <- match.arg(backend)

  # Set RNG seed if provided
  if (!is.null(seed)) set.seed(seed)

  # Resolve cores
  n_cores <- resolve_cores(n_cores, parallel)

  # Input validation
  x <- as.matrix(x)

  # Handle coords: either data.frame or spacc_dist
  if (inherits(coords, "spacc_dist")) {
    dist_mat <- as.matrix(coords)
    coord_data <- attr(coords, "coords")
  } else {
    stopifnot(
      "coords must have x and y columns" = all(c("x", "y") %in% names(coords))
    )
    coord_data <- coords
    dist_mat <- NULL
  }

  stopifnot(
    "x and coords must have same number of rows" = nrow(x) == nrow(coord_data),
    "n_seeds must be positive integer" = n_seeds > 0
  )

  # Handle groups: split species by group and recurse
  if (!is.null(groups)) {
    groups <- as.character(groups)
    stopifnot(
      "groups must have length equal to ncol(x)" = length(groups) == ncol(x)
    )
    group_levels <- unique(groups)
    if (progress) cli_info(sprintf("Running grouped accumulation (%d groups: %s)",
                                    length(group_levels), paste(group_levels, collapse = ", ")))
    objects <- lapply(group_levels, function(g) {
      cols <- which(groups == g)
      spacc(x[, cols, drop = FALSE], coords,
            n_seeds = n_seeds, method = method, distance = distance,
            backend = backend, support = support, include_halo = include_halo,
            sigma = sigma, cone_width = cone_width,
            parallel = parallel, n_cores = n_cores,
            progress = FALSE, groups = NULL,
            time = time, w_space = w_space, w_time = w_time,
            seed = seed)
    })
    names(objects) <- group_levels

    # Build grouped spacc directly
    base <- objects[[1]]
    return(structure(
      list(
        curves = stats::setNames(lapply(objects, `[[`, "curves"), group_levels),
        group_names = group_levels,
        coords = base$coords,
        n_seeds = base$n_seeds,
        n_sites = base$n_sites,
        n_species = stats::setNames(lapply(objects, `[[`, "n_species"), group_levels),
        method = base$method,
        distance = base$distance,
        backend = base$backend,
        sigma = base$sigma,
        cone_width = base$cone_width,
        time = base$time,
        w_space = base$w_space,
        w_time = base$w_time,
        support = base$support,
        call = match.call()
      ),
      class = "spacc"
    ))
  }

  # Handle support = "auto": split by country
  if (is.character(support) && length(support) == 1 && support == "auto") {
    check_suggests("sf")
    check_suggests("areaOfEffect")

    if (progress) cli_info("Auto-detecting countries from coordinates")

    coords_sf <- sf::st_as_sf(coord_data, coords = c("x", "y"), crs = 4326)
    aoe_result <- areaOfEffect::aoe(coords_sf)

    # Map support_id back to country names
    aoe_countries <- areaOfEffect::countries
    sid_to_name <- stats::setNames(aoe_countries$name, rownames(aoe_countries))
    unique_sids <- unique(aoe_result$support_id)
    country_names <- sid_to_name[unique_sids]

    if (progress) cli_info(sprintf("Found %d countries: %s",
      length(country_names), paste(country_names, collapse = ", ")))

    # Run spacc per country
    objects <- lapply(seq_along(unique_sids), function(i) {
      spacc(x, coords,
            n_seeds = n_seeds, method = method, distance = distance,
            backend = backend, support = country_names[i],
            include_halo = include_halo,
            sigma = sigma, cone_width = cone_width,
            parallel = parallel, n_cores = n_cores,
            progress = FALSE, groups = groups,
            time = time, w_space = w_space, w_time = w_time,
            seed = seed)
    })
    names(objects) <- country_names

    # Handle compound grouping (countries x species groups)
    has_groups <- !is.null(groups)
    if (has_groups) {
      # Flatten nested grouped objects into compound names
      all_curves <- list()
      all_n_species <- list()
      for (cname in country_names) {
        obj <- objects[[cname]]
        for (gname in obj$group_names) {
          compound <- paste(cname, gname, sep = ".")
          all_curves[[compound]] <- obj$curves[[gname]]
          all_n_species[[compound]] <- obj$n_species[[gname]]
        }
      }
      group_labels <- names(all_curves)
    } else {
      all_curves <- stats::setNames(lapply(objects, `[[`, "curves"), country_names)
      all_n_species <- stats::setNames(lapply(objects, `[[`, "n_species"), country_names)
      group_labels <- country_names
    }

    base <- objects[[1]]
    return(structure(
      list(
        curves = all_curves,
        group_names = group_labels,
        coords = coord_data,
        n_seeds = n_seeds,
        n_sites = nrow(x),
        n_species = all_n_species,
        method = base$method,
        distance = base$distance,
        backend = base$backend,
        sigma = base$sigma,
        cone_width = base$cone_width,
        time = time,
        w_space = if (!is.null(time)) w_space else NULL,
        w_time = if (!is.null(time)) w_time else NULL,
        support = list(auto = TRUE, countries = country_names, aoe_result = aoe_result),
        call = match.call()
      ),
      class = "spacc"
    ))
  }

  # Handle areaOfEffect support
  aoe_result <- NULL
  core_indices <- NULL
  original_indices <- seq_len(nrow(x))

  if (!is.null(support)) {
    check_suggests("sf")
    check_suggests("areaOfEffect")

    if (progress) cli_info("Processing spatial support via areaOfEffect")

    # Convert coords to sf if needed
    coords_sf <- sf::st_as_sf(coord_data, coords = c("x", "y"))

    # Get aoe classification
    if (inherits(support, "aoe_result")) {
      aoe_result <- support
    } else {
      aoe_result <- areaOfEffect::aoe(coords_sf, support = support)
    }

    # Get indices of core and halo sites
    core_mask <- aoe_result$aoe_class == "core"
    halo_mask <- aoe_result$aoe_class == "halo"

    # Map back to original indices (aoe may have pruned some points)
    point_ids <- as.integer(aoe_result$point_id)
    core_indices <- point_ids[core_mask]
    halo_indices <- point_ids[halo_mask]

    if (include_halo) {
      # Keep core + halo sites
      keep_indices <- c(core_indices, halo_indices)
    } else {
      # Keep core sites only (hard boundary)
      keep_indices <- core_indices
    }

    # Filter data to kept sites
    x <- x[keep_indices, , drop = FALSE]
    coord_data <- coord_data[keep_indices, , drop = FALSE]
    original_indices <- keep_indices

    # Update core_indices to be relative to filtered data
    core_indices <- which(keep_indices %in% core_indices)

    if (progress) {
      cli_info(sprintf("Support: %d core, %d halo sites (using %d total)",
                        sum(core_mask), sum(halo_mask), length(keep_indices)))
    }

    # Warn if core/halo are unbalanced
    n_core <- sum(core_mask)
    n_halo <- sum(halo_mask)
    if (n_core > 0 && n_halo > 0) {
      ratio <- max(n_core, n_halo) / min(n_core, n_halo)
      if (ratio > 3) {
        warning(sprintf(
          "Unbalanced core/halo: %d core vs %d halo (ratio %.1f:1). Consider adjusting support scale.",
          n_core, n_halo, ratio
        ), call. = FALSE)
      }
    } else if (n_core == 0) {
      stop("No core sites found. Check that your points intersect the support polygon.", call. = FALSE)
    } else if (n_halo == 0 && include_halo) {
      if (progress) cli_info("No halo sites found. All points are inside the support.")
    }
  }

  n_sites <- nrow(x)
  n_species_total <- ncol(x)

  # Handle spatiotemporal distance
  if (!is.null(time)) {
    stopifnot(
      "time must have length equal to nrow(x)" = length(time) == n_sites,
      "time must be numeric" = is.numeric(time),
      "w_space must be a positive number" = is.numeric(w_space) && length(w_space) == 1 && w_space > 0,
      "w_time must be a positive number" = is.numeric(w_time) && length(w_time) == 1 && w_time > 0
    )
    if (!method %in% c("knn", "radius", "gaussian")) {
      stop(sprintf("Spatiotemporal accumulation (time argument) is only supported for methods 'knn', 'radius', and 'gaussian', not '%s'.", method),
           call. = FALSE)
    }
    if (progress) cli_info("Computing spatiotemporal distance matrix")
    d_space <- cpp_distance_matrix(coord_data$x, coord_data$y, distance)
    d_time <- abs(outer(time, time, `-`))
    dist_mat <- w_space * d_space + w_time * d_time
  }

  # Convert to presence/absence if abundance
  species_pa <- (x > 0) * 1L
  storage.mode(species_pa) <- "integer"

  # Resolve backend for knn/kncn
  use_kdtree <- FALSE
  if (method %in% c("knn", "kncn")) {
    if (!is.null(time)) {
      # Spatiotemporal: force exact backend (kdtree can't handle composite distances)
      use_kdtree <- FALSE
    } else if (backend == "auto") {
      use_kdtree <- n_sites > 500L
    } else {
      use_kdtree <- backend == "kdtree"
    }
  }

  # Collector method: no simulation needed
  if (method == "collector") {
    curve <- cpp_collector_single(species_pa)
    curves <- matrix(curve, nrow = 1)
    n_seeds <- 1L
  } else {
    # Compute distance matrix if needed (exact knn, radius, gaussian)
    needs_dist <- (!use_kdtree && method == "knn") || method %in% c("radius", "gaussian")
    if (is.null(dist_mat) && needs_dist) {
      if (progress) cli_info(sprintf("Computing distances (%d x %d)", n_sites, n_sites))
      dist_mat <- cpp_distance_matrix(coord_data$x, coord_data$y, distance)
    }

    # Auto-calculate sigma for Gaussian method
    if (method == "gaussian" && is.null(sigma)) {
      if (is.null(dist_mat)) {
        dist_mat <- cpp_distance_matrix(coord_data$x, coord_data$y, distance)
      }
      sigma <- stats::median(dist_mat[dist_mat > 0])
    }

    # Run accumulation curves
    backend_label <- if (method %in% c("knn", "kncn")) {
      if (use_kdtree) "kdtree" else "exact"
    } else {
      NA_character_
    }
    if (method %in% c("knn", "kncn")) {
      if (progress) cli_info(sprintf("Running %s accumulation (%d seeds, %d cores, %s backend)",
                                      method, n_seeds, n_cores, backend_label))
    } else {
      if (progress) cli_info(sprintf("Running %s accumulation (%d seeds, %d cores)", method, n_seeds, n_cores))
    }

    # Sample seeds (from core sites only if support provided)
    if (!is.null(core_indices) && length(core_indices) > 0) {
      seed_pool <- core_indices - 1L
      explicit_seeds <- sample(seed_pool, n_seeds, replace = TRUE)
    } else {
      explicit_seeds <- NULL
    }

    curves <- switch(method,
      knn = if (use_kdtree) {
        if (is.null(explicit_seeds)) {
          cpp_knn_kdtree_parallel(species_pa, coord_data$x, coord_data$y, n_seeds, n_cores, progress, distance)
        } else {
          cpp_knn_kdtree_parallel_seeds(species_pa, coord_data$x, coord_data$y, explicit_seeds, n_cores, progress, distance)
        }
      } else {
        if (is.null(explicit_seeds)) {
          cpp_knn_parallel(species_pa, dist_mat, n_seeds, n_cores, progress)
        } else {
          cpp_knn_parallel_seeds(species_pa, dist_mat, explicit_seeds, n_cores, progress)
        }
      },
      kncn = if (use_kdtree) {
        cpp_kncn_kdtree_parallel(species_pa, coord_data$x, coord_data$y, n_seeds, n_cores, progress, distance)
      } else {
        cpp_kncn_parallel(species_pa, coord_data$x, coord_data$y, n_seeds, n_cores, progress)
      },
      random = cpp_random_parallel(species_pa, n_seeds, n_cores, progress),
      radius = cpp_radius_parallel(species_pa, dist_mat, n_seeds, n_cores, progress),
      gaussian = cpp_gaussian_parallel(species_pa, dist_mat, n_seeds, sigma, n_cores, progress),
      cone = cpp_cone_parallel(species_pa, coord_data$x, coord_data$y, n_seeds, cone_width, n_cores, progress)
    )
  }

  if (progress) cli_success("Done")

  structure(
    list(
      curves = curves,
      coords = coord_data,
      n_seeds = n_seeds,
      n_sites = n_sites,
      n_species = n_species_total,
      method = method,
      distance = distance,
      backend = backend_label,
      sigma = sigma,
      cone_width = if (method == "cone") cone_width else NULL,
      time = time,
      w_space = if (!is.null(time)) w_space else NULL,
      w_time = if (!is.null(time)) w_time else NULL,
      support = if (!is.null(support)) list(
        aoe_result = aoe_result,
        include_halo = include_halo,
        n_core = length(core_indices),
        n_halo = n_sites - length(core_indices),
        original_indices = original_indices
      ) else NULL,
      call = match.call()
    ),
    class = "spacc"
  )
}


#' Wavefront Expansion Accumulation
#'
#' Accumulate species within an expanding radius from seed points.
#' Models invasion spread from introduction points.
#'
#' @param x A site-by-species matrix.
#' @param coords A data.frame with x and y columns, or a `spacc_dist` object.
#' @param n_seeds Integer. Number of random starting points.
#' @param r0 Numeric. Initial radius. Default 0.
#' @param dr Numeric. Radius increment per step. Default auto-calculated.
#' @param n_steps Integer. Number of expansion steps. Default 50.
#' @param distance Character. Distance method.
#' @param progress Logical. Show progress?
#' @param seed Integer. Random seed.
#'
#' @return An object of class `spacc_wavefront` containing:
#'   \item{curves}{Matrix of species counts (n_seeds x n_steps)}
#'   \item{radius}{Vector of radius values}
#'   \item{sites_included}{Matrix of sites included at each step}
#'
#' @examples
#' \dontrun{
#' wf <- wavefront(species, coords, n_seeds = 20, n_steps = 100)
#' plot(wf)
#' }
#'
#' @export
wavefront <- function(x,
                      coords,
                      n_seeds = 50L,
                      r0 = 0,
                      dr = NULL,
                      n_steps = 50L,
                      distance = c("euclidean", "haversine"),
                      progress = TRUE,
                      seed = NULL) {

  distance <- match.arg(distance)

  if (!is.null(seed)) set.seed(seed)

  x <- as.matrix(x)

  # Handle coords
  if (inherits(coords, "spacc_dist")) {
    dist_mat <- as.matrix(coords)
    coord_data <- attr(coords, "coords")
  } else {
    stopifnot("coords must have x and y columns" = all(c("x", "y") %in% names(coords)))
    coord_data <- coords
    dist_mat <- cpp_distance_matrix(coord_data$x, coord_data$y, distance)
  }

  n_sites <- nrow(x)
  n_species_total <- ncol(x)

  # Auto-calculate dr if not provided
  if (is.null(dr)) {
    max_dist <- max(dist_mat)
    dr <- max_dist / n_steps
  }

  # Convert to presence/absence
  species_pa <- (x > 0) * 1L
  storage.mode(species_pa) <- "integer"

  if (progress) cli_info(sprintf("Running wavefront expansion (%d seeds, %d steps)", n_seeds, n_steps))

  result <- cpp_wavefront_parallel(species_pa, dist_mat, n_seeds, r0, dr, n_steps, 1L, progress)

  if (progress) cli_success("Done")

  structure(
    list(
      curves = result$curves,
      radius = result$radius,
      sites_included = result$sites_included,
      coords = coord_data,
      n_seeds = n_seeds,
      n_sites = n_sites,
      n_species = n_species_total,
      r0 = r0,
      dr = dr,
      n_steps = n_steps,
      distance = distance,
      call = match.call()
    ),
    class = "spacc_wavefront"
  )
}


#' @export
print.spacc_wavefront <- function(x, ...) {
  cat(sprintf("spacc wavefront: %d sites, %d species, %d seeds\n",
              x$n_sites, x$n_species, x$n_seeds))
  cat(sprintf("Radius: %.2f to %.2f (%d steps)\n",
              x$r0, max(x$radius), x$n_steps))
  invisible(x)
}


#' @export
summary.spacc_wavefront <- function(object, ci_level = 0.95, ...) {
  alpha <- (1 - ci_level) / 2

  mean_curve <- colMeans(object$curves)
  mean_sites <- colMeans(object$sites_included)
  lower <- apply(object$curves, 2, stats::quantile, probs = alpha)
  upper <- apply(object$curves, 2, stats::quantile, probs = 1 - alpha)

  data.frame(
    radius = object$radius,
    mean_species = mean_curve,
    lower = lower,
    upper = upper,
    mean_sites = mean_sites
  )
}


#' @export
plot.spacc_wavefront <- function(x, ci = TRUE, ci_alpha = 0.3,
                                  xaxis = c("radius", "sites"), ...) {
  check_suggests("ggplot2")

  xaxis <- match.arg(xaxis)
  summ <- summary(x)

  if (xaxis == "radius") {
    df <- data.frame(
      x = summ$radius,
      mean = summ$mean_species,
      lower = summ$lower,
      upper = summ$upper
    )
    xlab <- "Radius"
  } else {
    df <- data.frame(
      x = summ$mean_sites,
      mean = summ$mean_species,
      lower = summ$lower,
      upper = summ$upper
    )
    xlab <- "Sites included"
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = mean))

  if (ci) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = ci_alpha, fill = "#4CAF50"
    )
  }

  p +
    ggplot2::geom_line(linewidth = 1, color = "#2E7D32") +
    ggplot2::labs(
      x = xlab,
      y = "Cumulative species",
      title = "Wavefront Species Accumulation",
      subtitle = sprintf("%d seeds, r0=%.1f, dr=%.2f", x$n_seeds, x$r0, x$dr)
    ) +
    ggplot2::theme_minimal(base_size = 12)
}


#' Distance-Decay Analysis
#'
#' Analyze how species richness changes with distance from focal points.
#'
#' @param x A site-by-species matrix.
#' @param coords A data.frame with x and y columns.
#' @param n_seeds Integer. Number of focal points.
#' @param breaks Numeric vector. Distance thresholds. Default auto-calculated.
#' @param distance Character. Distance method.
#' @param progress Logical. Show progress?
#' @param seed Integer. Random seed.
#'
#' @return An object of class `spacc_decay` with distance-species relationships.
#'
#' @export
distanceDecay <- function(x,
                           coords,
                           n_seeds = 50L,
                           breaks = NULL,
                           distance = c("euclidean", "haversine"),
                           progress = TRUE,
                           seed = NULL) {

  distance <- match.arg(distance)

  if (!is.null(seed)) set.seed(seed)

  x <- as.matrix(x)
  stopifnot("coords must have x and y columns" = all(c("x", "y") %in% names(coords)))

  n_sites <- nrow(x)

  # Compute distance matrix
  dist_mat <- cpp_distance_matrix(coords$x, coords$y, distance)

  # Auto-calculate breaks if not provided
  if (is.null(breaks)) {
    max_dist <- max(dist_mat)
    breaks <- seq(0, max_dist, length.out = 51)[-1]  # exclude 0
  }

  species_pa <- (x > 0) * 1L
  storage.mode(species_pa) <- "integer"

  if (progress) cli_info(sprintf("Computing distance-decay (%d seeds, %d distances)", n_seeds, length(breaks)))

  curves <- cpp_distance_decay_parallel(species_pa, dist_mat, n_seeds, breaks, 1L, progress)

  if (progress) cli_success("Done")

  structure(
    list(
      curves = curves,
      breaks = breaks,
      n_seeds = n_seeds,
      n_sites = n_sites,
      n_species = ncol(x),
      distance = distance,
      call = match.call()
    ),
    class = "spacc_decay"
  )
}


#' @export
print.spacc_decay <- function(x, ...) {
  cat(sprintf("spacc distance-decay: %d seeds, %d distance breaks\n",
              x$n_seeds, length(x$breaks)))
  cat(sprintf("Distance range: %.2f to %.2f\n", min(x$breaks), max(x$breaks)))
  invisible(x)
}


#' @export
plot.spacc_decay <- function(x, ci = TRUE, ci_alpha = 0.3, ...) {
  check_suggests("ggplot2")

  mean_curve <- colMeans(x$curves)
  lower <- apply(x$curves, 2, stats::quantile, 0.025)
  upper <- apply(x$curves, 2, stats::quantile, 0.975)

  df <- data.frame(
    distance = x$breaks,
    mean = mean_curve,
    lower = lower,
    upper = upper
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = distance, y = mean))

  if (ci) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = ci_alpha, fill = "#4CAF50"
    )
  }

  p +
    ggplot2::geom_line(linewidth = 1, color = "#2E7D32") +
    ggplot2::labs(
      x = "Distance from focal point",
      y = "Cumulative species",
      title = "Distance-Decay Relationship",
      subtitle = sprintf("%d focal points", x$n_seeds)
    ) +
    ggplot2::theme_minimal(base_size = 12)
}


#' Resolve number of cores
#' @noRd
resolve_cores <- function(n_cores, parallel) {
  if (!parallel) return(1L)
  if (is.null(n_cores)) {
    max(1L, parallel::detectCores() - 1L)
  } else {
    as.integer(n_cores)
  }
}
