#' Spatial Phylogenetic Diversity Accumulation
#'
#' Compute spatial accumulation of phylogenetic diversity metrics (MPD, MNTD, PD).
#'
#' @param x A site-by-species matrix.
#' @param coords A data.frame with columns `x` and `y`, or a `spacc_dist` object.
#' @param tree A phylogenetic tree of class `phylo` (from ape package), or
#'   a pairwise phylogenetic distance matrix.
#' @param metric Character vector. Metrics to compute:
#'   - `"mpd"`: Mean Pairwise Distance
#'   - `"mntd"`: Mean Nearest Taxon Distance
#'   - `"pd"`: Faith's Phylogenetic Diversity (requires tree, not distance matrix)
#'   - `"evoheritage"`: EvoHeritage (Rosindell et al. 2024), PD with
#'     attrition of evolutionary heritage at rate `rho` (requires tree)
#'   - `"rao"`: Rao's quadratic entropy (abundance-weighted mean pairwise
#'     phylogenetic distance). Pass abundance data for weighting; with
#'     presence/absence it reduces to the equal-weight form.
#' @param rho Non-negative numeric vector. Attrition rate(s) for
#'   `"evoheritage"`, per unit of branch length. `rho = 0` gives Faith's PD.
#'   Several values give one curve each, named `evoheritage_rho<value>`.
#' @param n_seeds Integer. Number of random starting points. Default 50.
#' @param method Character. Accumulation method: `"knn"` or `"nn_walk"`.
#' @param distance Character. Site distance method: `"euclidean"` or `"haversine"`.
#' @param parallel Logical. Use parallel processing? Default `TRUE`.
#' @param n_cores Integer. Number of cores.
#' @param progress Logical. Show progress? Default `TRUE`.
#' @param seed Integer. Random seed.
#' @param map Logical. If `TRUE`, run accumulation from every site as seed
#'   and store per-site final values for spatial mapping. Enables
#'   [as_sf()] and `plot(type = "map")`. Default `FALSE`.
#' @param focal_points Optional focal points for `method = "knn"`. See [spacc()].
#' @param focal_domain Optional polygonal focal domain. See [spacc()].
#'
#' @return An object of class `spacc_phylo` containing:
#'   \item{curves}{Named list of matrices, one per metric (n_seeds x n_sites)}
#'   \item{metric}{Metrics computed}
#'   \item{coords, n_seeds, n_sites, method}{Parameters used}
#'
#' @details
#' Phylogenetic diversity metrics incorporate evolutionary relationships:
#'
#' - **MPD (Mean Pairwise Distance)**: Average phylogenetic distance between
#'   all pairs of species. Sensitive to tree-wide patterns.
#'
#' - **MNTD (Mean Nearest Taxon Distance)**: Average distance to closest
#'   relative. Sensitive to terminal clustering.
#'
#' - **PD (Faith's Phylogenetic Diversity)**: Total branch length connecting
#'   species to the root. Requires full tree object.
#'
#' - **EvoHeritage**: Features accumulate along each branch and are lost at
#'   rate `rho`, so a branch of length `L` passes on a fraction `exp(-rho * L)`
#'   of what entered it. A feature counts once if it survives to at least one
#'   present species. The value is in standardised units (the heritage at the
#'   end of a single branch of unit length), so `rho = 0` gives PD and large
#'   `rho` tends to species richness. Scanning `rho` places a community
#'   between its richness and PD curves.
#'
#' Tree-based metrics match the columns of `x` to `tree$tip.label` by name;
#' the tree may contain species absent from `x`.
#'
#' @references
#' Faith, D.P. (1992). Conservation evaluation and phylogenetic diversity.
#' Biological Conservation, 61, 1-10.
#'
#' Rosindell, J., Manson, K., Gumbs, R., Pearse, W.D. & Steel, M. (2024).
#' Phylogenetic biodiversity metrics should account for both accumulation and
#' attrition of evolutionary heritage. Systematic Biology, 73, 158-182.
#' \doi{10.1093/sysbio/syad072}
#'
#' Webb, C.O. (2000). Exploring the phylogenetic structure of ecological
#' communities: an example for rain forest trees. American Naturalist, 156, 145-155.
#'
#' @seealso `picante::mpd()`, `picante::mntd()`, `picante::pd()`
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   # Create random tree
#'   tree <- ape::rtree(30)
#'
#'   coords <- data.frame(x = runif(50), y = runif(50))
#'   species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)
#'   colnames(species) <- tree$tip.label
#'
#'   phylo <- spaccPhylo(species, coords, tree, metric = c("mpd", "mntd"))
#'   plot(phylo)
#' }
#' }
#'
#' @export
spaccPhylo <- function(x,
                       coords,
                       tree,
                       metric = c("mpd", "mntd"),
                       rho = 0,
                       n_seeds = 50L,
                       method = c("knn", "nn_walk"),
                       distance = c("euclidean", "haversine"),
                       parallel = TRUE,
                       n_cores = NULL,
                       progress = TRUE,
                       seed = NULL,
                       map = FALSE,
                       focal_points = NULL,
                       focal_domain = NULL) {

  method <- match.arg(method)
  distance <- match.arg(distance)
  metric <- match.arg(metric, c("mpd", "mntd", "pd", "evoheritage", "rao"),
                      several.ok = TRUE)
  stopifnot("rho must be non-negative and finite" =
              is.numeric(rho) && length(rho) >= 1 && all(is.finite(rho) & rho >= 0))

  if (!is.null(seed)) set.seed(seed)

  n_cores <- resolve_cores(n_cores, parallel)

  x <- as.matrix(x)
  site <- .resolve_coords(coords)
  coord_data <- site$coords
  stopifnot("x and coords must have same rows" = nrow(x) == nrow(coord_data))

  n_sites <- nrow(x)
  n_species <- ncol(x)
  ordering <- .accumulation_orders(method, coord_data, n_seeds, distance,
                                   site$dist_mat, focal_points, focal_domain,
                                   n_cores = n_cores)

  phylo_dist_mat <- .phylo_distances(tree, colnames(x))
  tree_metrics <- intersect(metric, c("pd", "evoheritage"))
  if (is.matrix(tree) && length(tree_metrics) > 0) {
    warning(sprintf("%s requires a full tree object; removing it from metrics",
                    paste(sQuote(tree_metrics, FALSE), collapse = " and ")))
    metric <- setdiff(metric, tree_metrics)
    tree_metrics <- character(0)
  }
  stopifnot("Phylo distance matrix must match species" = ncol(phylo_dist_mat) == n_species)

  # One curve per rho value; a single rate keeps the plain label.
  if ("evoheritage" %in% metric && length(rho) > 1) {
    at <- match("evoheritage", metric)
    metric <- append(metric[-at], paste0("evoheritage_rho", rho), after = at - 1L)
  }

  # Keep abundances as doubles. Presence metrics (mpd/mntd/pd) read
  # cumulative > 0, so binarising is unnecessary; Rao uses the abundance
  # weights and supports non-integer values (e.g. cover, biomass).
  species_pa <- x
  storage.mode(species_pa) <- "double"

  if (progress) cli_info(sprintf("Computing phylogenetic diversity (%s, %d seeds)",
                                  paste(metric, collapse = ", "), n_seeds))

  # Tree data for the metrics that walk the edges (PD, EvoHeritage)
  tree_edge <- NULL
  tree_edge_length <- NULL
  tree_n_tips <- 0L
  species_tip <- integer(0)
  if (length(tree_metrics) > 0) {
    if (is.null(tree$edge.length)) stop("tree must have edge lengths for PD and EvoHeritage")
    tree_edge <- tree$edge
    tree_edge_length <- tree$edge.length
    tree_n_tips <- length(tree$tip.label)
    species_tip <- .species_tips(tree, colnames(x), n_species)
  }
  evo_rho <- if ("evoheritage" %in% tree_metrics) as.numeric(rho) else numeric(0)

  run <- function(orders) {
    if (length(metric) == 0) return(list())
    cpp_phylo_order_parallel(species_pa, orders - 1L, phylo_dist_mat, metric, n_cores,
                             tree_edge, tree_edge_length, tree_n_tips,
                             species_tip, evo_rho)
  }

  result <- run(ordering$orders)

  if (progress) cli_success("Done")

  site_values <- NULL
  if (map && length(metric) > 0) {
    if (progress) cli_info("Computing per-site phylo map values (all sites as seeds)")
    site_values <- .site_map_values(method, coord_data, distance, site$dist_mat, n_cores,
                                    function(orders) lapply(run(orders), .final_per_seed))
    if (progress) cli_success("Map values computed")
  }

  structure(
    list(
      curves = result,
      metric = metric,
      rho = if (length(evo_rho) > 0) evo_rho else NULL,
      coords = coord_data,
      site_values = site_values,
      n_seeds = ordering$n_seeds,
      n_sites = n_sites,
      n_species = n_species,
      method = method,
      distance = distance,
      focal_points = ordering$focal_points,
      call = match.call()
    ),
    class = "spacc_phylo"
  )
}


# 0-based tip index of each species column. Species are matched to tips by
# name; an unnamed community matrix must list the tips in tree order.
.species_tips <- function(tree, species_names, n_species) {
  if (is.null(species_names)) {
    if (n_species != length(tree$tip.label)) {
      stop("x has no column names, so it must have one column per tree tip, in tip order")
    }
    return(seq_len(n_species) - 1L)
  }
  tip <- match(species_names, tree$tip.label)
  if (anyNA(tip)) stop("Some species in x are not in the tree")
  tip - 1L
}


# Pairwise phylogenetic distances between the species in `species_names`,
# from a phylo tree (cophenetic distances) or a supplied matrix.
.phylo_distances <- function(tree, species_names) {
  if (inherits(tree, "phylo")) {
    check_suggests("ape")
    phylo_dist <- as.matrix(ape::cophenetic.phylo(tree))
    if (!is.null(species_names)) {
      if (!all(species_names %in% tree$tip.label)) {
        stop("Some species in x are not in the tree")
      }
      phylo_dist <- phylo_dist[species_names, species_names]
    }
    return(phylo_dist)
  }
  if (is.matrix(tree)) return(tree)
  stop("tree must be a phylo object or distance matrix")
}


# Reorder a species-by-traits matrix to the species in `species_names`.
.match_traits <- function(traits, species_names, n_species) {
  traits <- as.matrix(traits)
  if (!is.null(species_names) && !is.null(rownames(traits))) {
    if (!all(species_names %in% rownames(traits))) {
      stop("Some species in x are not in traits")
    }
    traits <- traits[species_names, , drop = FALSE]
  }
  stopifnot("Traits must have one row per species" = nrow(traits) == n_species)
  traits
}


#' Spatial Functional Diversity Accumulation
#'
#' Compute spatial accumulation of functional diversity metrics based on traits.
#'
#' @param x A site-by-species matrix (abundance data recommended).
#' @param coords A data.frame with columns `x` and `y`, or a `spacc_dist` object.
#' @param traits A species-by-traits matrix. Row names should match species (columns of x).
#' @param metric Character vector. Metrics to compute:
#'   - `"fdis"`: Functional Dispersion (mean distance to centroid)
#'   - `"fric"`: Functional Richness (convex hull volume approximation)
#'   - `"rao"`: Rao's quadratic entropy (abundance-weighted mean pairwise
#'     Euclidean trait distance)
#' @param n_seeds Integer. Number of random starting points. Default 50.
#' @param method Character. Accumulation method: `"knn"` or `"nn_walk"`.
#' @param distance Character. Site distance method: `"euclidean"` or `"haversine"`.
#' @param parallel Logical. Use parallel processing? Default `TRUE`.
#' @param n_cores Integer. Number of cores.
#' @param progress Logical. Show progress? Default `TRUE`.
#' @param seed Integer. Random seed.
#' @param map Logical. If `TRUE`, run accumulation from every site as seed
#'   and store per-site final values for spatial mapping. Enables
#'   [as_sf()] and `plot(type = "map")`. Default `FALSE`.
#' @param focal_points Optional focal points for `method = "knn"`. See [spacc()].
#' @param focal_domain Optional polygonal focal domain. See [spacc()].
#'
#' @return An object of class `spacc_func` containing:
#'   \item{curves}{Named list of matrices, one per metric (n_seeds x n_sites)}
#'   \item{metric}{Metrics computed}
#'   \item{coords, n_seeds, n_sites, method}{Parameters used}
#'
#' @details
#' Functional diversity metrics quantify trait space occupation:
#'
#' - **FDis (Functional Dispersion)**: Abundance-weighted mean distance from
#'   the community centroid in trait space. Captures functional divergence.
#'
#' - **FRic (Functional Richness)**: Volume of trait space occupied (convex hull).
#'   Requires more species than traits to compute.
#'
#' @references
#' Laliberté, E. & Legendre, P. (2010). A distance-based framework for measuring
#' functional diversity from multiple traits. Ecology, 91, 299-305.
#'
#' @seealso `FD::dbFD()` for comprehensive functional diversity analysis
#'
#' @examples
#' \donttest{
#' coords <- data.frame(x = runif(50), y = runif(50))
#' species <- matrix(rpois(50 * 20, 2), nrow = 50)
#'
#' # Trait matrix (species x traits)
#' traits <- matrix(rnorm(20 * 3), nrow = 20)
#' rownames(traits) <- paste0("sp", 1:20)
#' colnames(species) <- rownames(traits)
#'
#' func <- spaccFunc(species, coords, traits, metric = c("fdis", "fric"))
#' plot(func)
#' }
#'
#' @export
spaccFunc <- function(x,
                      coords,
                      traits,
                      metric = c("fdis", "fric"),
                      n_seeds = 50L,
                      method = c("knn", "nn_walk"),
                      distance = c("euclidean", "haversine"),
                      parallel = TRUE,
                      n_cores = NULL,
                      progress = TRUE,
                      seed = NULL,
                      map = FALSE,
                      focal_points = NULL,
                      focal_domain = NULL) {

  method <- match.arg(method)
  distance <- match.arg(distance)
  metric <- match.arg(metric, c("fdis", "fric", "rao"), several.ok = TRUE)

  if (!is.null(seed)) set.seed(seed)

  n_cores <- resolve_cores(n_cores, parallel)

  x <- as.matrix(x)
  site <- .resolve_coords(coords)
  coord_data <- site$coords
  stopifnot("x and coords must have same rows" = nrow(x) == nrow(coord_data))

  n_sites <- nrow(x)
  n_species <- ncol(x)
  ordering <- .accumulation_orders(method, coord_data, n_seeds, distance,
                                   site$dist_mat, focal_points, focal_domain,
                                   n_cores = n_cores)

  traits <- .match_traits(traits, colnames(x), n_species)

  # Keep abundance data as doubles (supports non-integer abundances,
  # e.g. cover or biomass, for FDis and Rao weighting)
  storage.mode(x) <- "double"

  if (progress) cli_info(sprintf("Computing functional diversity (%s, %d seeds)",
                                  paste(metric, collapse = ", "), n_seeds))

  run <- function(orders) {
    cpp_func_order_parallel(x, orders - 1L, traits, metric, n_cores)
  }

  result <- run(ordering$orders)

  if (progress) cli_success("Done")

  site_values <- NULL
  if (map) {
    if (progress) cli_info("Computing per-site functional map values (all sites as seeds)")
    site_values <- .site_map_values(method, coord_data, distance, site$dist_mat, n_cores,
                                    function(orders) lapply(run(orders), .final_per_seed))
    if (progress) cli_success("Map values computed")
  }

  structure(
    list(
      curves = result,
      metric = metric,
      traits = traits,
      coords = coord_data,
      site_values = site_values,
      n_seeds = ordering$n_seeds,
      n_sites = n_sites,
      n_species = n_species,
      n_traits = ncol(traits),
      method = method,
      distance = distance,
      focal_points = ordering$focal_points,
      call = match.call()
    ),
    class = "spacc_func"
  )
}


# Shared curve/map plot for the per-metric accumulation classes.
.plot_metric_curves <- function(x, type, ci, ci_alpha, metric, point_size, palette,
                                map_title, ylab, title) {
  if (type == "map") {
    if (is.null(metric)) metric <- x$metric[1]
    return(.plot_site_map(x, metric,
                          title = sprintf(map_title, toupper(metric)),
                          subtitle = sprintf("%d sites, %s method", x$n_sites, x$method),
                          point_size = point_size, palette = palette))
  }
  summ <- summary(x)
  summ$metric <- factor(toupper(summ$metric))
  .plot_curves_by(summ, "metric", "Metric",
                  xlab = "Sites accumulated", ylab = ylab, title = title,
                  ci = ci, ci_alpha = ci_alpha)
}


# S3 methods for spacc_phylo -----------------------------------------------

#' @export
print.spacc_phylo <- function(x, ...) {
  .print_curves(x, "spacc phylogenetic diversity",
                sprintf("Metrics: %s", paste(x$metric, collapse = ", ")))
}


#' @export
summary.spacc_phylo <- function(object, ci_level = 0.95, ...) {
  .curve_frame(object$curves, "metric", object$metric, ci_level)
}


#' @export
plot.spacc_phylo <- function(x, type = c("curve", "map"), ci = TRUE, ci_alpha = 0.2,
                              metric = NULL, point_size = 3, palette = "viridis", ...) {
  .plot_metric_curves(x, match.arg(type), ci, ci_alpha, metric, point_size, palette,
                      map_title = "Phylogenetic diversity map: %s",
                      ylab = "Phylogenetic diversity",
                      title = "Spatial Phylogenetic Diversity Accumulation")
}


#' @export
as_sf.spacc_phylo <- function(x, crs = NULL) .as_sf_site_values(x, crs)


# S3 methods for spacc_func -----------------------------------------------

#' @export
print.spacc_func <- function(x, ...) {
  cat(sprintf("spacc functional diversity: %d sites, %d species, %d traits, %d seeds\n",
              x$n_sites, x$n_species, x$n_traits, x$n_seeds))
  cat(sprintf("Metrics: %s\n", paste(x$metric, collapse = ", ")))
  invisible(x)
}


#' @export
summary.spacc_func <- function(object, ci_level = 0.95, ...) {
  .curve_frame(object$curves, "metric", object$metric, ci_level)
}


#' @export
plot.spacc_func <- function(x, type = c("curve", "map"), ci = TRUE, ci_alpha = 0.2,
                             metric = NULL, point_size = 3, palette = "viridis", ...) {
  .plot_metric_curves(x, match.arg(type), ci, ci_alpha, metric, point_size, palette,
                      map_title = "Functional diversity map: %s",
                      ylab = "Functional diversity",
                      title = "Spatial Functional Diversity Accumulation")
}


#' @export
as_sf.spacc_func <- function(x, crs = NULL) .as_sf_site_values(x, crs)
