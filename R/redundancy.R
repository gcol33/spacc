#' Spatial Phylogenetic or Functional Redundancy Accumulation
#'
#' Track how much of the taxonomic diversity added along the spatial
#' accumulation is new phylogenetic or functional diversity, as the ratio of
#' the taxonomic Hill number to the phylogenetic or functional Hill number of
#' the accumulated community.
#'
#' @param x A site-by-species matrix (abundance or presence/absence).
#' @param coords A data.frame with columns `x` and `y`, or a `spacc_dist` object.
#' @param tree A phylogenetic tree of class `phylo` with branch lengths.
#'   Supply either `tree` or `traits`.
#' @param traits A species-by-traits matrix or data.frame, row names matching
#'   the columns of `x`.
#' @param q Numeric vector. Orders of diversity. Default `c(0, 1, 2)`.
#' @param n_seeds Integer. Number of random starting points. Default 50.
#' @param method Character. Accumulation method: `"knn"` or `"nn_walk"`.
#' @param distance Character. Site distance method: `"euclidean"` or `"haversine"`.
#' @param dist_method Character. Trait distance for `traits`: `"euclidean"`
#'   on standardised traits (default) or `"gower"`.
#' @param normalize Logical. Rescale trait distances to a maximum of 1.
#'   Default `TRUE`.
#' @param parallel Logical. Use parallel processing? Default `TRUE`.
#' @param n_cores Integer. Number of cores.
#' @param progress Logical. Show progress? Default `TRUE`.
#' @param seed Integer. Random seed.
#' @param focal_points Optional focal points for `method = "knn"`. See [spacc()].
#' @param focal_domain Optional polygonal focal domain. See [spacc()].
#'
#' @return An object of class `spacc_redundancy` containing:
#'   \item{curves}{Named list (one per q) of n_seeds x n_sites matrices of the
#'     ratio taxonomic / phylogenetic or functional Hill number}
#'   \item{taxonomic}{The taxonomic Hill numbers, same layout}
#'   \item{dimension}{The phylogenetic or functional Hill numbers, same layout}
#'   \item{type}{`"phylogenetic"` or `"functional"`}
#'   \item{q, coords, n_seeds, n_sites, n_species, method, orders}{Parameters
#'     and the site orders used}
#'
#' @details
#' At each accumulation step the function computes, for every order q, the
#' taxonomic Hill number \eqn{^qD} of the accumulated community and either
#' its phylogenetic Hill number \eqn{^qD(T)} (Chao et al. 2010; effective
#' number of equally abundant, equally distinct lineages, with \eqn{T} the
#' abundance-weighted mean root-to-tip depth) or its similarity-sensitive
#' functional Hill number \eqn{^qD^Z} (Leinster & Cobbold 2012; with
#' \eqn{Z = 1 - d} from trait distances). The ratio
#' \deqn{R_q = {}^qD / {}^qD(T) \quad \textrm{or} \quad {}^qD / {}^qD^Z}
#' is at least 1, equal to 1 when species are maximally distinct (a star
#' tree, or every pair of species at the maximal trait distance), and grows
#' as the species added are close relatives or functionally similar to those
#' already present. A ratio rising along the curve says the area is adding
#' species faster than lineages or trait combinations.
#'
#' @references
#' Chao, A., Chiu, C.-H. & Jost, L. (2010). Phylogenetic diversity measures
#' based on Hill numbers. Philosophical Transactions of the Royal Society B,
#' 365, 3599-3609.
#'
#' Leinster, T. & Cobbold, C.A. (2012). Measuring diversity: the importance
#' of species similarity. Ecology, 93, 477-489.
#'
#' Alberdi, A. & Gilbert, M.T.P. (2019). A guide to the application of Hill
#' numbers to DNA-based diversity analyses. Molecular Ecology Resources, 19,
#' 804-817. \doi{10.1111/1755-0998.13014}
#'
#' @seealso [spaccHill()], [spaccPhylo()], [spaccFunc()], [diversityProfile()]
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   tree <- ape::rcoal(25)
#'   coords <- data.frame(x = runif(40), y = runif(40))
#'   species <- matrix(rpois(40 * 25, 1), nrow = 40,
#'                     dimnames = list(NULL, tree$tip.label))
#'   red <- spaccRedundancy(species, coords, tree = tree, n_seeds = 10)
#'   plot(red)
#' }
#' }
#'
#' @export
spaccRedundancy <- function(x,
                            coords,
                            tree = NULL,
                            traits = NULL,
                            q = c(0, 1, 2),
                            n_seeds = 50L,
                            method = c("knn", "nn_walk"),
                            distance = c("euclidean", "haversine"),
                            dist_method = c("euclidean", "gower"),
                            normalize = TRUE,
                            parallel = TRUE,
                            n_cores = NULL,
                            progress = TRUE,
                            seed = NULL,
                            focal_points = NULL,
                            focal_domain = NULL) {

  method <- match.arg(method)
  distance <- match.arg(distance)
  dist_method <- match.arg(dist_method)
  if (is.null(tree) == is.null(traits)) {
    stop("Supply either `tree` (phylogenetic) or `traits` (functional).", call. = FALSE)
  }
  stopifnot("q must be numeric and >= 0" = is.numeric(q) && all(q >= 0))
  if (!is.null(seed)) set.seed(seed)
  n_cores <- resolve_cores(n_cores, parallel)

  x <- as.matrix(x)
  storage.mode(x) <- "double"
  site <- .resolve_coords(coords)
  coord_data <- site$coords
  stopifnot("x and coords must have same rows" = nrow(x) == nrow(coord_data))
  n_species <- ncol(x)

  if (!is.null(tree)) {
    if (!inherits(tree, "phylo")) stop("tree must be a phylo object", call. = FALSE)
    if (is.null(tree$edge.length)) stop("tree must have branch lengths", call. = FALSE)
    type <- "phylogenetic"
    dimension_args <- list(tree_edge = tree$edge, tree_edge_length = tree$edge.length,
                           tree_n_tips = length(tree$tip.label),
                           species_tip = .species_tips(tree, colnames(x), n_species))
  } else {
    type <- "functional"
    dimension_args <- list(similarity_mat = .trait_similarity(traits, colnames(x), dist_method,
                                                              normalize))
  }

  if (progress) cli_info(sprintf("Computing %s redundancy (q = %s, %d seeds)",
                                  type, paste(q, collapse = ", "), n_seeds))
  ordering <- .accumulation_orders(method, coord_data, n_seeds, distance,
                                   site$dist_mat, focal_points, focal_domain,
                                   n_cores = n_cores)
  result <- do.call(cpp_redundancy_order_parallel,
                    c(list(x, ordering$orders - 1L, q, n_cores), dimension_args))
  if (progress) cli_success("Done")

  structure(
    list(
      curves = result$redundancy,
      taxonomic = result$taxonomic,
      dimension = result$dimension,
      type = type,
      q = q,
      coords = coord_data,
      n_seeds = ordering$n_seeds,
      n_sites = nrow(x),
      n_species = n_species,
      method = method,
      distance = distance,
      orders = ordering$orders,
      focal_points = ordering$focal_points,
      call = match.call()
    ),
    class = "spacc_redundancy"
  )
}


#' @export
print.spacc_redundancy <- function(x, ...) {
  final <- vapply(x$curves, function(m) mean(m[, x$n_sites]), numeric(1))
  .print_curves(x, sprintf("spacc %s redundancy", x$type),
                sprintf("Orders (q): %s", paste(x$q, collapse = ", ")),
                sprintf("  q = %s: mean final redundancy = %.2f", x$q, final))
}


#' @export
summary.spacc_redundancy <- function(object, ci_level = 0.95, ...) {
  comps <- list(redundancy = object$curves, taxonomic = object$taxonomic,
                dimension = object$dimension)
  do.call(rbind, lapply(names(comps), function(comp) {
    cbind(component = comp, .curve_frame(comps[[comp]], "q", object$q, ci_level))
  }))
}


#' @export
as.data.frame.spacc_redundancy <- function(x, row.names = NULL, optional = FALSE, ...) {
  summary(x)
}


#' @export
plot.spacc_redundancy <- function(x, ci = TRUE, ci_alpha = 0.2, ...) {
  summ <- summary(x)
  summ <- summ[summ$component == "redundancy", ]
  summ$q_label <- .q_label(summ$q)
  .plot_curves_by(summ, "q_label", "Order",
                  xlab = "Sites accumulated",
                  ylab = "Redundancy ratio",
                  title = sprintf("Spatial %s Redundancy", tools::toTitleCase(x$type)),
                  ci = ci, ci_alpha = ci_alpha) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50")
}
