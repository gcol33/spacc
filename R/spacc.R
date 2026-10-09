#' Spatial Species Accumulation Curves
#'
#' Compute species accumulation curves using various spatial sampling methods
#' with C++ backend for performance.
#'
#' @param x A site-by-species matrix (rows = sites, cols = species) with
#'   presence/absence (0/1) or abundance data. Can also be a data.frame.
#' @param coords Site coordinates. Can be:
#'   - A data.frame with columns `x` and `y`
#'   - An `sf` object with point geometries (CRS is preserved)
#'   - A `spacc_dist` object from [distances()]
#' @param n_seeds Integer. Number of random focal points or starting sites for
#'   uncertainty quantification. Default 50.
#' @param method Character. Accumulation method:
#'   - `"knn"`: fixed-focus spatially constrained rarefaction
#'   - `"kncn"`: k-Nearest Centroid Neighbor (visit closest to centroid)
#'   - `"nn_walk"`: nearest-neighbour walk (move from the current site to the
#'     closest unvisited site)
#'   - `"random"`: Random order (null model)
#'   - `"gaussian"`: Probabilistic selection weighted by distance
#'   - `"cone"`: Directional expansion within angular constraint
#'   - `"collector"`: Sites in data order (no randomization, single curve)
#' @param order Optional user-defined accumulation order(s). When supplied,
#'   `method` is ignored and sites are accumulated in the given sequence(s)
#'   (no distance computation or seed sampling). Can be:
#'   - A numeric vector of length `nrow(x)`: a single ordering (one curve).
#'   - A list of such vectors, or a matrix with one ordering per row: each
#'     produces a curve, treated like seeds for uncertainty quantification.
#'   Each ordering must be a permutation of `seq_len(nrow(x))`. Cannot be
#'   combined with `support` or `time`. Default `NULL`.
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
#' @param backend Character. Nearest-neighbor backend for `nn_walk` and `kncn`:
#'   - `"auto"` (default): Uses exact (brute-force) for <=500 sites,
#'     spatial tree for >500 sites.
#'   - `"exact"`: Always use brute-force with precomputed distance matrix.
#'   - `"kdtree"`: Always use spatial tree. Uses k-d tree (nanoflann) for
#'     Euclidean distances and ball tree for haversine distances. Faster for
#'     large datasets, no distance matrix needed.
#'   A supplied distance matrix (`spacc_dist` coordinates or the
#'   spatiotemporal composite built from `time`) is always used as given, so
#'   those inputs take the exact backend.
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
#'   distance matrix: `"nn_walk"`, `"gaussian"`.
#' @param w_space Numeric. Weight for spatial distance when `time` is provided.
#'   Default 1.
#' @param w_time Numeric. Weight for temporal distance when `time` is provided.
#'   Default 1.
#' @param seed Integer. Random seed for reproducibility. Default `NULL`.
#' @param focal_points Optional data frame or `sf` point object with `x` and
#'   `y` coordinates. For `method = "knn"`, each row is a fixed focal point and
#'   sites are accumulated by increasing distance from it. When supplied,
#'   `n_seeds` is set to the number of focal points. Default `NULL` samples
#'   continuous focal points uniformly from the spatial domain.
#' @param focal_domain Optional `sf` or `sfc` polygon defining the domain from
#'   which continuous `knn` focal points are sampled. When omitted, the convex
#'   hull of the eligible site coordinates is used. For haversine distances,
#'   hull sampling uses a local equal-area projection. Supplying the study
#'   polygon preserves concavities, holes, and disconnected components in
#'   irregular sampling domains. Its coordinates must use the same coordinate
#'   reference system as `coords`.
#'
#' @return When `groups = NULL`, an object of class `spacc` containing:
#'   \item{curves}{Matrix of cumulative species counts (n_seeds x n_sites)}
#'   \item{orders}{Matrix of 1-based site indices (n_seeds x n_sites): the
#'     order in which each curve visited the sites}
#'   \item{coords}{Original coordinates}
#'   \item{n_seeds}{Number of seeds used}
#'   \item{method}{Method used}
#'   \item{n_species}{Total species in dataset}
#'   \item{focal_points}{Continuous focal points used by `method = "knn"`}
#'
#' @details
#' The canonical `knn` method follows Chiarucci et al. (2009): each curve uses
#' one fixed focal point and accumulates sites by increasing distance from that
#' point. The `nn_walk` method is a greedy traversal whose reference point moves
#' to the most recently selected site. The two methods represent different
#' spatial sampling models.
#'
#' Every method is a rule for ordering the sites; the species data are then
#' accumulated along those orders. The `orders` element of the result lets the
#' same traversals be reused, for example through the `order` argument.
#'
#' @examples
#' \donttest{
#' coords <- data.frame(x = runif(50), y = runif(50))
#' species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)
#'
#' # Basic usage
#' sac <- spacc(species, coords)
#' plot(sac)
#'
#' # Different methods
#' sac_knn <- spacc(species, coords, method = "knn")
#' sac_rand <- spacc(species, coords, method = "random")
#' comp <- compare(sac_knn, sac_rand)
#' }
#'
#' @references
#' Arrhenius, O. (1921). Species and area. Journal of Ecology, 9, 95-99.
#'
#' Scheiner, S.M. (2003). Six types of species-area curves. Global Ecology
#' and Biogeography, 12, 441-447.
#'
#' Chiarucci, A., Bacaro, G., Rocchini, D., Ricotta, C., Palmer, M.W. &
#' Scheiner, S.M. (2009). Spatially constrained rarefaction: incorporating the
#' autocorrelated structure of biological communities into sample-based
#' rarefaction. Community Ecology, 10, 209-214.
#'
#' @export
spacc <- function(x,
                  coords,
                  n_seeds = 50L,
                  method = c("knn", "kncn", "nn_walk", "random", "gaussian", "cone", "collector"),
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
                  seed = NULL,
                  focal_points = NULL,
                  focal_domain = NULL,
                  order = NULL) {

  method <- match.arg(method)
  distance <- match.arg(distance)
  backend <- match.arg(backend)
  call <- match.call()

  if (!is.null(seed)) set.seed(seed)
  n_cores <- resolve_cores(n_cores, parallel)

  x <- as.matrix(x)
  if (!is.numeric(x)) {
    stop("Species matrix must be numeric (presence/absence or abundance data)")
  }

  site <- .resolve_coords(coords)
  coord_data <- site$coords
  dist_mat <- site$dist_mat

  stopifnot(
    "x and coords must have same number of rows" = nrow(x) == nrow(coord_data),
    "n_seeds must be positive integer" = n_seeds > 0
  )

  if (!is.null(order)) {
    if (!is.null(support)) {
      stop("`order` cannot be combined with `support`; user-defined orderings index the full dataset.", call. = FALSE)
    }
    if (!is.null(time)) {
      stop("`order` cannot be combined with `time`; user-defined orderings already fix the accumulation sequence.", call. = FALSE)
    }
  }
  if ((!is.null(focal_points) || !is.null(focal_domain)) && method != "knn") {
    stop("`focal_points` and `focal_domain` are only available for `method = \"knn\"`.", call. = FALSE)
  }

  # Species groups: one curve set per group along the same spatial orders
  if (!is.null(groups)) {
    groups <- as.character(groups)
    stopifnot("groups must have length equal to ncol(x)" = length(groups) == ncol(x))
    group_levels <- unique(groups)
    if (progress) cli_info(sprintf("Running grouped accumulation (%d groups: %s)",
                                    length(group_levels), paste(group_levels, collapse = ", ")))
    shared_seed <- if (is.null(seed)) sample.int(.Machine$integer.max, 1L) else seed
    objects <- lapply(group_levels, function(g) {
      spacc(x[, groups == g, drop = FALSE], coords,
            n_seeds = n_seeds, method = method, distance = distance,
            backend = backend, support = support, include_halo = include_halo,
            sigma = sigma, cone_width = cone_width,
            parallel = parallel, n_cores = n_cores,
            progress = FALSE, groups = NULL,
            time = time, w_space = w_space, w_time = w_time,
            seed = shared_seed, focal_points = focal_points,
            focal_domain = focal_domain, order = order)
    })
    return(.grouped_spacc(objects, group_levels, inner_first = TRUE, call = call))
  }

  # support = "auto": one curve set per detected country
  if (is.character(support) && length(support) == 1 && support == "auto") {
    check_suggests("sf")
    check_suggests("areaOfEffect")

    if (progress) cli_info("Auto-detecting countries from coordinates")

    coords_sf <- sf::st_as_sf(coord_data, coords = c("x", "y"), crs = 4326)
    aoe_result <- areaOfEffect::aoe(coords_sf)

    aoe_countries <- areaOfEffect::countries
    sid_to_name <- stats::setNames(aoe_countries$name, rownames(aoe_countries))
    country_names <- unname(sid_to_name[unique(aoe_result$support_id)])

    if (progress) cli_info(sprintf("Found %d countries: %s",
      length(country_names), paste(country_names, collapse = ", ")))

    objects <- lapply(country_names, function(country) {
      spacc(x, coords,
            n_seeds = n_seeds, method = method, distance = distance,
            backend = backend, support = country,
            include_halo = include_halo,
            sigma = sigma, cone_width = cone_width,
            parallel = parallel, n_cores = n_cores,
            progress = FALSE, groups = NULL,
            time = time, w_space = w_space, w_time = w_time,
            seed = seed, focal_points = focal_points,
            focal_domain = focal_domain)
    })
    out <- .grouped_spacc(objects, country_names, inner_first = FALSE, call = call,
                          coords = coord_data, n_sites = nrow(x), n_seeds = n_seeds)
    out$time <- time
    out$w_space <- if (!is.null(time)) w_space else NULL
    out$w_time <- if (!is.null(time)) w_time else NULL
    out$support <- list(auto = TRUE, countries = country_names, aoe_result = aoe_result)
    return(out)
  }

  # Explicit support: keep core (+ halo) sites, draw seeds from the core
  seed_pool <- NULL
  support_info <- NULL
  if (!is.null(support)) {
    kept <- .apply_support(support, coord_data, site$crs, distance, include_halo, progress)
    x <- x[kept$keep, , drop = FALSE]
    coord_data <- coord_data[kept$keep, , drop = FALSE]
    seed_pool <- kept$core
    support_info <- list(
      aoe_result = kept$aoe_result,
      include_halo = include_halo,
      n_core = length(kept$core),
      n_halo = length(kept$keep) - length(kept$core),
      original_indices = kept$keep
    )
  }

  n_sites <- nrow(x)

  # Spatiotemporal composite distance
  if (!is.null(time)) {
    stopifnot(
      "time must have length equal to nrow(x)" = length(time) == n_sites,
      "time must be numeric" = is.numeric(time),
      "w_space must be a positive number" = is.numeric(w_space) && length(w_space) == 1 && w_space > 0,
      "w_time must be a positive number" = is.numeric(w_time) && length(w_time) == 1 && w_time > 0
    )
    if (!method %in% c("nn_walk", "gaussian")) {
      stop(sprintf("Spatiotemporal accumulation (time argument) is only supported for methods 'nn_walk' and 'gaussian', not '%s'.", method),
           call. = FALSE)
    }
    if (progress) cli_info("Computing spatiotemporal distance matrix")
    d_space <- cpp_distance_matrix(coord_data$x, coord_data$y, distance)
    dist_mat <- w_space * d_space + w_time * abs(outer(time, time, `-`))
  }

  species_pa <- (x > 0) * 1L
  storage.mode(species_pa) <- "integer"

  if (!is.null(order)) {
    orders <- .build_orders(order, n_sites)
    ordering <- list(orders = orders, n_seeds = nrow(orders), focal_points = NULL,
                     backend = NA_character_, sigma = NULL)
    method <- "user"
    if (progress) cli_info(sprintf("Running user-defined accumulation (%d ordering%s)",
                                    nrow(orders), if (nrow(orders) == 1L) "" else "s"))
  } else {
    ordering <- .accumulation_orders(method, coord_data, n_seeds, distance,
                                     dist_mat, focal_points, focal_domain,
                                     seed_pool = seed_pool, backend = backend,
                                     sigma = sigma, cone_width = cone_width,
                                     n_cores = n_cores)
    if (progress) {
      backend_note <- if (is.na(ordering$backend)) "" else sprintf(", %s backend", ordering$backend)
      cli_info(sprintf("Running %s accumulation (%d seeds, %d cores%s)",
                       method, ordering$n_seeds, n_cores, backend_note))
    }
  }

  curves <- cpp_order_parallel(species_pa, ordering$orders - 1L, n_cores)

  if (progress) cli_success("Done")

  .new_spacc(curves, coord_data, species_pa, ordering, method, distance,
             cone_width = if (method == "cone") cone_width else NULL,
             time = time, w_space = w_space, w_time = w_time,
             support = support_info, call = call)
}


# Core/halo site selection through areaOfEffect. Returns the kept site indices,
# the core indices relative to the kept set, and the aoe classification.
.apply_support <- function(support, coord_data, crs, distance, include_halo, progress) {
  check_suggests("sf")
  check_suggests("areaOfEffect")

  if (progress) cli_info("Processing spatial support via areaOfEffect")

  # Haversine coordinates without a CRS are taken as WGS84
  coords_sf <- if (!is.null(crs) && !is.na(crs)) {
    sf::st_as_sf(coord_data, coords = c("x", "y"), crs = crs)
  } else if (distance == "haversine") {
    sf::st_as_sf(coord_data, coords = c("x", "y"), crs = 4326)
  } else {
    sf::st_as_sf(coord_data, coords = c("x", "y"))
  }

  aoe_result <- if (inherits(support, "aoe_result")) {
    support
  } else {
    areaOfEffect::aoe(coords_sf, support = support)
  }

  core_mask <- aoe_result$aoe_class == "core"
  halo_mask <- aoe_result$aoe_class == "halo"
  point_ids <- as.integer(aoe_result$point_id)
  core_ids <- point_ids[core_mask]
  halo_ids <- point_ids[halo_mask]
  keep <- if (include_halo) c(core_ids, halo_ids) else core_ids

  n_core <- sum(core_mask)
  n_halo <- sum(halo_mask)
  if (progress) {
    cli_info(sprintf("Support: %d core, %d halo sites (using %d total)",
                     n_core, n_halo, length(keep)))
  }
  if (n_core == 0) {
    stop("No core sites found. Check that your points intersect the support polygon.", call. = FALSE)
  }
  if (n_halo > 0) {
    ratio <- max(n_core, n_halo) / min(n_core, n_halo)
    if (ratio > 3) {
      warning(sprintf(
        "Unbalanced core/halo: %d core vs %d halo (ratio %.1f:1). Consider adjusting support scale.",
        n_core, n_halo, ratio
      ), call. = FALSE)
    }
  } else if (include_halo && progress) {
    cli_info("No halo sites found. All points are inside the support.")
  }

  list(keep = keep, core = which(keep %in% core_ids), aoe_result = aoe_result)
}


# Assemble an ungrouped spacc object.
.new_spacc <- function(curves, coords, species_pa, ordering, method, distance,
                       cone_width = NULL, time = NULL, w_space = NULL, w_time = NULL,
                       support = NULL, call = NULL) {
  structure(
    list(
      curves = curves,
      orders = ordering$orders,
      coords = coords,
      incidence_freq = colSums(species_pa),
      n_seeds = ordering$n_seeds,
      n_sites = nrow(species_pa),
      n_species = ncol(species_pa),
      method = method,
      distance = distance,
      backend = ordering$backend,
      sigma = ordering$sigma,
      cone_width = cone_width,
      time = time,
      w_space = if (!is.null(time)) w_space else NULL,
      w_time = if (!is.null(time)) w_time else NULL,
      support = support,
      focal_points = ordering$focal_points,
      call = call
    ),
    class = "spacc"
  )
}


# Combine named spacc objects into a grouped spacc. Objects that are already
# grouped are flattened to compound names: `inner.outer` when `inner_first`,
# else `outer.inner` (so country x group always reads "country.group").
.grouped_spacc <- function(objects, names, inner_first = TRUE, call = NULL,
                           coords = NULL, n_sites = NULL, n_seeds = NULL) {
  names(objects) <- names
  base <- objects[[1]]

  if (is_grouped(base)) {
    curves <- list()
    n_species <- list()
    for (outer in names) {
      obj <- objects[[outer]]
      for (inner in obj$group_names) {
        compound <- if (inner_first) paste(inner, outer, sep = ".") else paste(outer, inner, sep = ".")
        curves[[compound]] <- obj$curves[[inner]]
        n_species[[compound]] <- obj$n_species[[inner]]
      }
    }
  } else {
    curves <- stats::setNames(lapply(objects, `[[`, "curves"), names)
    n_species <- stats::setNames(lapply(objects, `[[`, "n_species"), names)
  }

  structure(
    list(
      curves = curves,
      group_names = names(curves),
      coords = if (is.null(coords)) base$coords else coords,
      n_seeds = if (is.null(n_seeds)) base$n_seeds else n_seeds,
      n_sites = if (is.null(n_sites)) base$n_sites else n_sites,
      n_species = n_species,
      method = base$method,
      distance = base$distance,
      backend = base$backend,
      sigma = base$sigma,
      cone_width = base$cone_width,
      time = base$time,
      w_space = base$w_space,
      w_time = base$w_time,
      support = base$support,
      focal_points = base$focal_points,
      call = call
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
#' @references
#' Shigesada, N. & Kawasaki, K. (1997). Biological Invasions: Theory and
#' Practice. Oxford University Press.
#'
#' @examples
#' \donttest{
#' coords <- data.frame(x = runif(50), y = runif(50))
#' species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)
#' wf <- spaccWavefront(species, coords, n_seeds = 20, n_steps = 50)
#' plot(wf)
#' }
#'
#' @export
spaccWavefront <- function(x,
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
  site <- .resolve_coords(coords)
  coord_data <- site$coords
  dist_mat <- site$dist_mat
  if (is.null(dist_mat)) {
    dist_mat <- cpp_distance_matrix(coord_data$x, coord_data$y, distance)
  }

  n_sites <- nrow(x)
  if (is.null(dr)) dr <- max(dist_mat) / n_steps
  radius <- r0 + (seq_len(n_steps) - 1) * dr

  species_pa <- (x > 0) * 1L
  storage.mode(species_pa) <- "integer"

  if (progress) cli_info(sprintf("Running wavefront expansion (%d seeds, %d steps)", n_seeds, n_steps))

  seeds <- sample.int(n_sites, n_seeds, replace = TRUE)
  result <- cpp_radius_parallel(species_pa, dist_mat, seeds - 1L, radius)

  if (progress) cli_success("Done")

  structure(
    list(
      curves = result$curves,
      radius = radius,
      sites_included = result$sites_included,
      coords = coord_data,
      n_seeds = n_seeds,
      n_sites = n_sites,
      n_species = ncol(x),
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
  .print_curves(x, "spacc wavefront",
                sprintf("Radius: %.2f to %.2f (%d steps)", x$r0, max(x$radius), x$n_steps))
}


#' @export
summary.spacc_wavefront <- function(object, ci_level = 0.95, ...) {
  stats <- .curve_stats(object$curves, ci_level)
  data.frame(
    radius = object$radius,
    mean_species = stats$mean,
    lower = stats$lower,
    upper = stats$upper,
    mean_sites = colMeans(object$sites_included)
  )
}


#' @export
plot.spacc_wavefront <- function(x, ci = TRUE, ci_alpha = 0.3,
                                  xaxis = c("radius", "sites"), ...) {
  xaxis <- match.arg(xaxis)
  summ <- summary(x)
  .plot_curve(summ,
              x = if (xaxis == "radius") "radius" else "mean_sites",
              y = "mean_species",
              xlab = if (xaxis == "radius") "Radius" else "Sites included",
              ylab = "Cumulative species",
              title = "Wavefront Species Accumulation",
              subtitle = sprintf("%d seeds, r0=%.1f, dr=%.2f", x$n_seeds, x$r0, x$dr),
              ci = ci, ci_alpha = ci_alpha)
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
#' @references
#' Nekola, J.C. & White, P.S. (1999). The distance decay of similarity in
#' biogeography and ecology. Journal of Biogeography, 26, 867-878.
#'
#' Soininen, J., McDonald, R. & Hillebrand, H. (2007). The distance decay of
#' similarity in ecological communities. Ecography, 30, 3-12.
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
  site <- .resolve_coords(coords)
  coord_data <- site$coords
  dist_mat <- site$dist_mat
  if (is.null(dist_mat)) {
    dist_mat <- cpp_distance_matrix(coord_data$x, coord_data$y, distance)
  }

  n_sites <- nrow(x)
  if (is.null(breaks)) {
    breaks <- seq(0, max(dist_mat), length.out = 51)[-1]  # exclude 0
  }

  species_pa <- (x > 0) * 1L
  storage.mode(species_pa) <- "integer"

  if (progress) cli_info(sprintf("Computing distance-decay (%d seeds, %d distances)", n_seeds, length(breaks)))

  seeds <- sample.int(n_sites, n_seeds, replace = TRUE)
  curves <- cpp_radius_parallel(species_pa, dist_mat, seeds - 1L, breaks)$curves

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
  summ <- cbind(distance = x$breaks, .curve_stats(x$curves))
  .plot_curve(summ, x = "distance",
              xlab = "Distance from focal point", ylab = "Cumulative species",
              title = "Distance-Decay Relationship",
              subtitle = sprintf("%d focal points", x$n_seeds),
              ci = ci, ci_alpha = ci_alpha)
}


#' Validate and normalise user-supplied accumulation orders
#'
#' Accepts a single ordering vector, a list of vectors, or a matrix (one
#' ordering per row) and returns a 1-based integer matrix with one ordering
#' per row. Each ordering must be a permutation of `seq_len(n_sites)`.
#' @noRd
.build_orders <- function(order, n_sites) {
  if (is.matrix(order)) {
    orders_mat <- order
  } else if (is.list(order)) {
    lens <- vapply(order, length, integer(1))
    if (any(lens != n_sites)) {
      stop(sprintf("Each ordering in `order` must have length %d (number of sites)", n_sites),
           call. = FALSE)
    }
    orders_mat <- do.call(rbind, order)
  } else {
    orders_mat <- matrix(order, nrow = 1)
  }

  if (ncol(orders_mat) != n_sites) {
    stop(sprintf("`order` must contain %d site indices per ordering, got %d",
                 n_sites, ncol(orders_mat)), call. = FALSE)
  }

  expected <- seq_len(n_sites)
  for (r in seq_len(nrow(orders_mat))) {
    o <- orders_mat[r, ]
    if (anyNA(o) || !setequal(o, expected)) {
      stop(sprintf("Ordering %d in `order` must be a permutation of 1:%d (each site exactly once)",
                   r, n_sites), call. = FALSE)
    }
  }

  storage.mode(orders_mat) <- "integer"
  orders_mat
}


#' Resolve number of cores
#' @noRd
resolve_cores <- function(n_cores, parallel) {
  if (!parallel) return(1L)
  if (is.null(n_cores)) {
    n_cores <- max(1L, parallel::detectCores() - 1L)
  } else {
    n_cores <- as.integer(n_cores)
    if (n_cores < 1L) {
      stop("n_cores must be a positive integer")
    }
  }
  # R CMD check sets _R_CHECK_LIMIT_CORES_ and permits at most two cores.
  # Cap the count and pin the RcppParallel thread pool so the C++ workers,
  # which draw from the global pool, honour the limit.
  if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) {
    n_cores <- min(n_cores, 2L)
    RcppParallel::setThreadOptions(numThreads = n_cores)
  }
  n_cores
}
