#' Spatial Accumulation with Hill Numbers
#'
#' Compute spatial species accumulation curves using Hill numbers (effective
#' number of species) instead of raw richness. Hill numbers unify diversity
#' measures: q=0 is richness, q=1 is exponential Shannon, q=2 is inverse Simpson.
#'
#' @param x A site-by-species matrix (rows = sites, cols = species) with
#'   presence/absence (0/1) or abundance data.
#' @param coords A data.frame with columns `x` and `y` containing site coordinates,
#'   or a `spacc_dist` object from [distances()].
#' @param q Numeric vector. Orders of diversity to compute. Default `c(0, 1, 2)`.
#'   - q = 0: Species richness
#'   - q = 1: Exponential of Shannon entropy (effective common species)
#'   - q = 2: Inverse Simpson (effective dominant species)
#' @param n_seeds Integer. Number of random starting points. Default 50.
#' @param method Character. Accumulation method: `"knn"` or `"nn_walk"`.
#' @param distance Character. Distance method: `"euclidean"` or `"haversine"`.
#' @param parallel Logical. Use parallel processing? Default `TRUE`.
#' @param n_cores Integer. Number of cores. Default `NULL` uses `detectCores() - 1`.
#' @param progress Logical. Show progress bar? Default `TRUE`.
#' @param seed Integer. Random seed for reproducibility.
#' @param map Logical. If `TRUE`, run accumulation from every site as seed
#'   and store per-site final Hill numbers for spatial mapping. Enables
#'   [as_sf()] and `plot(type = "map")`. Default `FALSE`.
#' @param focal_points Optional focal points for `method = "knn"`. See [spacc()].
#' @param focal_domain Optional polygonal focal domain. See [spacc()].
#'
#' @return An object of class `spacc_hill` containing:
#'   \item{curves}{Named list of matrices, one per q value (n_seeds x n_sites)}
#'   \item{q}{Vector of q values used}
#'   \item{coords}{Original coordinates}
#'   \item{n_seeds}{Number of seeds}
#'   \item{n_sites}{Number of sites}
#'   \item{n_species}{Total species}
#'   \item{method}{Method used}
#'
#' @details
#' Hill numbers (Chao et al. 2014) provide a unified framework for diversity
#' measurement. Unlike raw richness (q=0), higher-order Hill numbers (q=1, q=2)
#' down-weight rare species, providing different perspectives on diversity.
#'
#' The spatial accumulation of Hill numbers can reveal scale-dependent diversity
#' patterns missed by richness alone.
#'
#' @references
#' Chao, A., Gotelli, N.J., Hsieh, T.C., Sander, E.L., Ma, K.H., Colwell, R.K.
#' & Ellison, A.M. (2014). Rarefaction and extrapolation with Hill numbers:
#' a framework for sampling and estimation in species diversity studies.
#' Ecological Monographs, 84, 45-67.
#'
#' @seealso [spacc()] for richness-only accumulation, [iNEXT::iNEXT()] for
#'   non-spatial Hill number rarefaction
#'
#' @examples
#' \donttest{
#' # Compare diversity at different orders
#' coords <- data.frame(x = runif(50), y = runif(50))
#' species <- matrix(rpois(50 * 30, 2), nrow = 50)
#'
#' hill <- spaccHill(species, coords, q = c(0, 1, 2))
#' plot(hill)
#'
#' # Extract summary at final site
#' summary(hill)
#' }
#'
#' @export
spaccHill <- function(x,
                      coords,
                      q = c(0, 1, 2),
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

  if (!is.null(seed)) set.seed(seed)

  n_cores <- resolve_cores(n_cores, parallel)

  x <- as.matrix(x)
  stopifnot(
    "q must be numeric and >= 0" = is.numeric(q) && all(q >= 0),
    "n_seeds must be positive" = n_seeds > 0
  )

  site <- .resolve_coords(coords)
  coord_data <- site$coords
  dist_mat <- site$dist_mat
  stopifnot("x and coords must have same rows" = nrow(x) == nrow(coord_data))

  n_sites <- nrow(x)
  n_species <- ncol(x)
  storage.mode(x) <- "integer"

  if (progress) cli_info(sprintf("Computing Hill numbers (q = %s, %d seeds)",
                                  paste(q, collapse = ", "), n_seeds))

  ordering <- .accumulation_orders(method, coord_data, n_seeds, distance,
                                   dist_mat, focal_points, focal_domain,
                                   n_cores = n_cores)
  curves <- cpp_order_hill_parallel(x, ordering$orders - 1L, q, n_cores)

  if (progress) cli_success("Done")

  site_values <- NULL
  if (map) {
    if (progress) cli_info("Computing per-site Hill map values (all sites as seeds)")
    site_values <- .site_map_values(method, coord_data, distance, dist_mat, n_cores,
      function(orders) {
        map_curves <- cpp_order_hill_parallel(x, orders - 1L, q, n_cores)
        stats::setNames(lapply(map_curves, .final_per_seed), paste0("q", q))
      })
    if (progress) cli_success("Map values computed")
  }

  structure(
    list(
      curves = curves,
      q = q,
      coords = coord_data,
      site_values = site_values,
      n_seeds = ordering$n_seeds,
      n_sites = n_sites,
      n_species = n_species,
      method = method,
      distance = distance,
      focal_points = ordering$focal_points,
      orders = ordering$orders,
      call = match.call()
    ),
    class = "spacc_hill"
  )
}


#' @export
print.spacc_hill <- function(x, ...) {
  .print_curves(x, "spacc Hill numbers",
                sprintf("Orders (q): %s", paste(x$q, collapse = ", ")),
                sprintf("Method: %s", x$method))
}


#' @export
summary.spacc_hill <- function(object, ci_level = 0.95, ...) {
  .curve_frame(object$curves, "q", object$q, ci_level)
}


#' @export
plot.spacc_hill <- function(x, type = c("curve", "map"), ci = TRUE, ci_alpha = 0.2,
                             q = NULL, point_size = 3, palette = "viridis", ...) {
  type <- match.arg(type)

  if (type == "map") {
    if (is.null(q)) q <- x$q[1]
    return(.plot_site_map(x, paste0("q", q),
                          title = sprintf("Hill number map (q = %s)", q),
                          subtitle = sprintf("%d sites, %s method", x$n_sites, x$method),
                          point_size = point_size, palette = palette))
  }

  summ <- summary(x)
  summ$q_label <- .q_label(summ$q)
  .plot_curves_by(summ, "q_label", "Order",
                  xlab = "Sites accumulated", ylab = "Hill number (effective species)",
                  title = "Spatial Hill Number Accumulation",
                  ci = ci, ci_alpha = ci_alpha)
}


#' @export
as_sf.spacc_hill <- function(x, crs = NULL) .as_sf_site_values(x, crs)


# ============================================================================
# HILL BETA DIVERSITY (Jost 2007 along accumulation)
# ============================================================================

#' Spatial Hill Number Beta Diversity
#'
#' Compute multisite beta diversity as gamma/alpha decomposition of Hill
#' numbers along the spatial accumulation curve (Jost 2007 framework).
#'
#' @param x A site-by-species matrix (rows = sites, cols = species).
#' @param coords A data.frame with columns `x` and `y`, or a `spacc_dist` object.
#' @param q Numeric vector. Orders of diversity. Default `c(0, 1, 2)`.
#' @param n_seeds Integer. Number of random starting points. Default 50.
#' @param distance Character. `"euclidean"` or `"haversine"`.
#' @param parallel Logical. Use parallel processing? Default `TRUE`.
#' @param n_cores Integer. Number of cores. Default `NULL`.
#' @param progress Logical. Show progress? Default `TRUE`.
#' @param seed Integer. Random seed.
#' @param method Character. Accumulation method: `"knn"` or `"nn_walk"`.
#' @param focal_points Optional focal points for `method = "knn"`. See [spacc()].
#' @param focal_domain Optional polygonal focal domain. See [spacc()].
#' @param standardize Character. `"none"` (default) decomposes the observed
#'   samples; `"coverage"` standardises gamma and alpha to a common sample
#'   coverage at every step (iNEXT.beta3D, Chao et al. 2023). Needs integer
#'   abundances.
#' @param coverage Numeric in (0, 1]. Fixed coverage for
#'   `standardize = "coverage"`, the same at every step. `NULL` (default)
#'   uses, at each step, the lower of the observed coverages of the pooled and
#'   joint samples.
#'
#' @return An object of class `spacc_hill_beta` containing:
#'   \item{gamma}{Named list of n_seeds x n_sites matrices (one per q)}
#'   \item{alpha}{Named list of n_seeds x n_sites matrices (one per q)}
#'   \item{beta}{Named list of n_seeds x n_sites matrices (one per q)}
#'   \item{dissimilarity}{With `standardize = "coverage"`: list of the four
#'     dissimilarities `1-CqN`, `1-UqN`, `1-VqN`, `1-SqN`, each a named list
#'     of n_seeds x n_sites matrices (one per q); `NULL` otherwise}
#'   \item{coverage}{With `standardize = "coverage"`: n_seeds x n_sites
#'     matrix of the coverage used at each step}
#'   \item{q}{Vector of q values}
#'   \item{coords}{Original coordinates}
#'   \item{n_seeds, n_sites, n_species}{Dimensions}
#'
#' @details
#' At each accumulation step k, the function computes:
#' - **Gamma**: Hill number of the pooled community (all k sites combined)
#' - **Alpha**: Generalized mean of per-site Hill numbers (Jost's power mean)
#' - **Beta**: gamma / alpha (effective number of distinct communities)
#'
#' Beta = 1 means all sites are identical; beta = k means all sites are
#' completely different. This provides the Hill-number analogue of the
#' Baselga-based `spaccBeta()`.
#'
#' @section Coverage standardisation:
#' Observed beta depends on sampling completeness: gamma and alpha approach
#' their asymptotes at different rates, so the ratio moves along the curve
#' even where composition does not change. With `standardize = "coverage"`,
#' gamma is the Hill number of the pooled sample and alpha the Hill number of
#' the joint sample (every species-by-site cell treated as a unit) divided by
#' k (Chiu et al. 2014), both rarefied or extrapolated to the same sample
#' coverage C with the estimators of Chao et al. (2014); beta is their ratio
#' (Chao et al. 2023). By default C is, at each step, the lower of the two
#' observed coverages, so both values come from the observed samples or
#' their rarefaction. A fixed `coverage` makes the curve comparable across
#' steps and datasets, at the cost of extrapolation where a step's sample is
#' less complete; for q = 0, Chao et al. (2023) advise against extrapolating
#' beyond the coverage of twice the sample size.
#'
#' The standardised beta is transformed to four dissimilarities between 0 and 1
#' (Chao et al. 2023): `1-CqN` (Sorensen-type, the classic Sorensen index at
#' q = 0), `1-UqN` (Jaccard-type), `1-VqN` and `1-SqN`. They are `NA` at the
#' first step, where there is one site.
#'
#' @references
#' Jost, L. (2007). Partitioning diversity into independent alpha and beta
#' components. Ecology, 88, 2427-2439.
#'
#' Chiu, C.-H., Jost, L. & Chao, A. (2014). Phylogenetic beta diversity,
#' similarity, and differentiation measures based on Hill numbers.
#' Ecological Monographs, 84, 21-44.
#'
#' Chao, A., Gotelli, N.J., Hsieh, T.C., et al. (2014). Rarefaction and
#' extrapolation with Hill numbers: a framework for sampling and estimation in
#' species diversity studies. Ecological Monographs, 84, 45-67.
#'
#' Chao, A., Thorn, S., Chiu, C.-H., et al. (2023). Rarefaction and
#' extrapolation with beta diversity under a framework of Hill numbers: the
#' iNEXT.beta3D standardization. Ecological Monographs, 93, e1588.
#' \doi{10.1002/ecm.1588}
#'
#' @seealso [spaccBeta()] for P/A-based Baselga partitioning,
#'   [spaccHill()] for Hill accumulation without beta decomposition
#'
#' @examples
#' \donttest{
#' coords <- data.frame(x = runif(40), y = runif(40))
#' species <- matrix(rpois(40 * 20, 2), nrow = 40)
#'
#' hb <- spaccHillBeta(species, coords, n_seeds = 10, progress = FALSE)
#' plot(hb)
#' }
#'
#' @export
spaccHillBeta <- function(x,
                           coords,
                           q = c(0, 1, 2),
                           n_seeds = 50L,
                           distance = c("euclidean", "haversine"),
                           parallel = TRUE,
                           n_cores = NULL,
                           progress = TRUE,
                           seed = NULL,
                           method = c("knn", "nn_walk"),
                           focal_points = NULL,
                           focal_domain = NULL,
                           standardize = c("none", "coverage"),
                           coverage = NULL) {

  method <- match.arg(method)
  distance <- match.arg(distance)
  standardize <- match.arg(standardize)
  if (!is.null(seed)) set.seed(seed)
  n_cores <- resolve_cores(n_cores, parallel)

  x <- as.matrix(x)
  stopifnot(
    "q must be numeric and >= 0" = is.numeric(q) && all(q >= 0),
    "n_seeds must be positive" = n_seeds > 0
  )
  if (!is.null(coverage)) {
    if (standardize != "coverage") stop("`coverage` requires standardize = \"coverage\"")
    stopifnot("coverage must be a single value in (0, 1]" =
                is.numeric(coverage) && length(coverage) == 1 && coverage > 0 && coverage <= 1)
  }
  if (standardize == "coverage") {
    if (any(abs(x - round(x)) > 0)) stop("standardize = \"coverage\" needs integer abundances")
    if (!is.null(coverage) && any(q != round(q))) {
      stop("A fixed `coverage` may require extrapolation, which is defined for integer q only")
    }
  }

  site <- .resolve_coords(coords)
  coord_data <- site$coords
  stopifnot("x and coords must have same rows" = nrow(x) == nrow(coord_data))

  n_sites <- nrow(x)
  n_species <- ncol(x)
  storage.mode(x) <- "integer"

  if (progress) cli_info(sprintf("Computing Hill beta (q = %s, %d seeds)",
                                  paste(q, collapse = ", "), n_seeds))

  ordering <- .accumulation_orders(method, coord_data, n_seeds, distance,
                                   site$dist_mat, focal_points, focal_domain,
                                   n_cores = n_cores)
  result <- cpp_order_hill_beta_parallel(x, ordering$orders - 1L, q, n_cores,
                                         standardize == "coverage",
                                         if (is.null(coverage)) NA_real_ else coverage)

  if (progress) cli_success("Done")

  dissimilarity <- NULL
  if (standardize == "coverage") {
    dissimilarity <- result[c("1-CqN", "1-UqN", "1-VqN", "1-SqN")]
  }

  structure(
    list(
      gamma = result$gamma,
      alpha = result$alpha,
      beta = result$beta,
      dissimilarity = dissimilarity,
      standardize = standardize,
      coverage = result$coverage,
      target_coverage = coverage,
      q = q,
      coords = coord_data,
      n_seeds = ordering$n_seeds,
      n_sites = n_sites,
      n_species = n_species,
      method = method,
      focal_points = ordering$focal_points,
      distance = distance,
      call = match.call()
    ),
    class = "spacc_hill_beta"
  )
}


# S3 methods for spacc_hill_beta ---------------------------------------------

#' @export
print.spacc_hill_beta <- function(x, ...) {
  final_beta <- vapply(seq_along(x$q), function(i) mean(x$beta[[i]][, x$n_sites]), numeric(1))
  standardized <- identical(x$standardize, "coverage")
  header <- sprintf("Orders (q): %s", paste(x$q, collapse = ", "))
  if (standardized) {
    header <- c(header, sprintf("Coverage-standardised (%s)",
                                if (is.null(x$target_coverage)) "per-step reference coverage"
                                else sprintf("C = %g", x$target_coverage)))
  }
  .print_curves(x, "spacc Hill beta diversity", header,
                sprintf("  q = %s: mean final beta = %.2f", x$q, final_beta))
}


# Components carried by a spacc_hill_beta object, as named lists of per-q
# curve matrices.
.hill_beta_components <- function(x) {
  c(list(gamma = x$gamma, alpha = x$alpha, beta = x$beta), x$dissimilarity)
}


#' @export
summary.spacc_hill_beta <- function(object, ci_level = 0.95, ...) {
  comps <- .hill_beta_components(object)
  do.call(rbind, lapply(names(comps), function(comp) {
    cbind(component = comp,
          .curve_frame(comps[[comp]], "q", object$q, ci_level))
  }))
}


#' @export
as.data.frame.spacc_hill_beta <- function(x, row.names = NULL, optional = FALSE, ...) {
  summary(x)
}


#' @export
plot.spacc_hill_beta <- function(x, component = c("beta", "gamma", "alpha", "all",
                                                  "1-CqN", "1-UqN", "1-VqN", "1-SqN"),
                                  ci = TRUE, ci_alpha = 0.2, ...) {
  component <- match.arg(component)
  summ <- summary(x)
  summ$q_label <- .q_label(summ$q)

  if (startsWith(component, "1-")) {
    if (!component %in% names(x$dissimilarity)) {
      stop(sprintf("Component '%s' needs standardize = \"coverage\"", component), call. = FALSE)
    }
    return(.plot_curves_by(summ[summ$component == component & !is.na(summ$mean), ],
                           "q_label", "Order",
                           xlab = "Sites accumulated", ylab = component,
                           title = sprintf("Coverage-standardised dissimilarity %s", component),
                           ci = ci, ci_alpha = ci_alpha))
  }
  if (component != "all") {
    return(.plot_curves_by(summ[summ$component == component, ], "q_label", "Order",
                           xlab = "Sites accumulated",
                           ylab = paste0(tools::toTitleCase(component), " diversity"),
                           title = sprintf("Hill %s Diversity along Accumulation",
                                           tools::toTitleCase(component)),
                           ci = ci, ci_alpha = ci_alpha))
  }
  summ <- summ[summ$component %in% c("gamma", "alpha", "beta"), ]

  summ$component <- factor(summ$component, levels = c("gamma", "alpha", "beta"))
  .plot_curves_by(summ, "component", "Component",
                  xlab = "Sites accumulated", ylab = "Diversity",
                  title = "Hill Number Diversity Partitioning",
                  ci = ci, ci_alpha = ci_alpha) +
    ggplot2::facet_wrap(~ q_label, scales = "free_y")
}
