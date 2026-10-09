#' Spatial Accumulation of Interaction Diversity
#'
#' Accumulate an interaction network over space: interaction richness and
#' interaction Hill numbers, sample coverage, the richness of the two
#' interacting levels (for example plants and pollinators) and connectance,
#' all along the same site orders.
#'
#' @param interactions A data frame of observed interactions with columns
#'   `site`, `lower`, `upper` and optionally `count` (interaction frequency,
#'   default 1 per row). `lower` and `upper` name the species of the two
#'   levels; an interaction is a `lower`-`upper` pair.
#' @param coords A data frame with columns `x` and `y`, one row per site. Sites
#'   in `interactions` are matched to `coords$site` if present, otherwise to
#'   the row names of `coords`, otherwise taken as row numbers.
#' @param q Numeric vector. Orders of diversity. Default `c(0, 1, 2)`.
#' @param n_seeds Integer. Number of random starting points. Default 50.
#' @param method Character. Accumulation method: `"knn"` or `"nn_walk"`.
#' @param distance Character. `"euclidean"` or `"haversine"`.
#' @param parallel Logical. Use parallel processing? Default `TRUE`.
#' @param n_cores Integer. Number of cores.
#' @param progress Logical. Show progress? Default `TRUE`.
#' @param seed Integer. Random seed.
#' @param focal_points Optional focal points for `method = "knn"`. See [spacc()].
#' @param focal_domain Optional polygonal focal domain. See [spacc()].
#'
#' @return An object of class `spacc_network` containing:
#'   \item{interaction, lower, upper}{Named lists (one per q) of n_seeds x
#'     n_sites matrices: Hill numbers of the interaction frequencies and of
#'     the summed frequencies of each level's species}
#'   \item{coverage}{List with the n_seeds x n_sites sample coverage matrices
#'     of `interaction`, `lower` and `upper`}
#'   \item{connectance}{n_seeds x n_sites matrix of interaction richness
#'     divided by the product of the two species richnesses}
#'   \item{q, coords, n_seeds, n_sites, n_interactions, n_lower, n_upper,
#'     method, orders}{Parameters, dimensions and site orders}
#'
#' @details
#' Each unique `lower`-`upper` pair is a unit of diversity and its frequency
#' its abundance, so the interaction Hill number of order q is the ordinary
#' Hill number of the interaction frequencies of the accumulated sites, the
#' effective number of equally frequent interactions, and the sample coverage
#' is the Chao & Jost (2012) coverage of those frequencies (Chiu et al. 2023).
#' The species curves of the two levels use each species' summed interaction
#' frequency. Connectance is the realised share of the possible interactions
#' between the species detected so far.
#'
#' @references
#' Chiu, C.-H., Chao, A., Vogel, S., Kriegel, P. & Thorn, S. (2023).
#' Quantifying and estimating ecological network diversity based on
#' incomplete sampling data. Philosophical Transactions of the Royal Society
#' B, 378, 20220183. \doi{10.1098/rstb.2022.0183}
#'
#' @seealso [spaccHillCoverage()] for the same quantities on species
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' coords <- data.frame(x = runif(30), y = runif(30))
#' interactions <- data.frame(
#'   site = sample(30, 400, replace = TRUE),
#'   lower = paste0("plant", sample(12, 400, replace = TRUE)),
#'   upper = paste0("bee", sample(20, 400, replace = TRUE))
#' )
#' net <- spaccNetwork(interactions, coords, n_seeds = 10)
#' net
#' plot(net)
#' }
#'
#' @export
spaccNetwork <- function(interactions,
                         coords,
                         q = c(0, 1, 2),
                         n_seeds = 50L,
                         method = c("knn", "nn_walk"),
                         distance = c("euclidean", "haversine"),
                         parallel = TRUE,
                         n_cores = NULL,
                         progress = TRUE,
                         seed = NULL,
                         focal_points = NULL,
                         focal_domain = NULL) {

  method <- match.arg(method)
  distance <- match.arg(distance)
  stopifnot("q must be numeric and >= 0" = is.numeric(q) && all(q >= 0))
  if (!is.null(seed)) set.seed(seed)
  n_cores <- resolve_cores(n_cores, parallel)

  site <- .resolve_coords(coords)
  coord_data <- site$coords
  mats <- .network_matrices(interactions, coords, nrow(coord_data))

  if (progress) cli_info(sprintf("Accumulating %d interactions between %d and %d species (%d seeds)",
                                  ncol(mats$interaction), ncol(mats$lower), ncol(mats$upper),
                                  n_seeds))
  ordering <- .accumulation_orders(method, coord_data, n_seeds, distance,
                                   site$dist_mat, focal_points, focal_domain,
                                   n_cores = n_cores)
  orders <- ordering$orders - 1L
  q_run <- sort(unique(c(0, q)))
  run <- function(m) cpp_order_hill_coverage_parallel(m, orders, q_run, n_cores)
  inter <- run(mats$interaction)
  lower <- run(mats$lower)
  upper <- run(mats$upper)
  if (progress) cli_success("Done")

  keep <- match(q, q_run)
  connectance <- inter$hills[[1]] / (lower$hills[[1]] * upper$hills[[1]])

  structure(
    list(
      interaction = inter$hills[keep],
      lower = lower$hills[keep],
      upper = upper$hills[keep],
      coverage = list(interaction = inter$coverage, lower = lower$coverage,
                      upper = upper$coverage),
      connectance = connectance,
      q = q,
      coords = coord_data,
      n_seeds = ordering$n_seeds,
      n_sites = nrow(coord_data),
      n_interactions = ncol(mats$interaction),
      n_lower = ncol(mats$lower),
      n_upper = ncol(mats$upper),
      method = method,
      distance = distance,
      orders = ordering$orders,
      focal_points = ordering$focal_points,
      call = match.call()
    ),
    class = "spacc_network"
  )
}


# Site x interaction, site x lower-species and site x upper-species frequency
# matrices from a long interaction table.
.network_matrices <- function(interactions, coords, n_sites) {
  interactions <- as.data.frame(interactions)
  needed <- c("site", "lower", "upper")
  if (!all(needed %in% names(interactions))) {
    stop("interactions needs columns site, lower and upper", call. = FALSE)
  }
  count <- if ("count" %in% names(interactions)) interactions$count else rep(1, nrow(interactions))
  if (any(!is.finite(count) | count < 0 | count != round(count))) {
    stop("count must be non-negative integers", call. = FALSE)
  }

  site_ids <- if (is.data.frame(coords) && "site" %in% names(coords)) {
    coords$site
  } else if (is.data.frame(coords) && !is.null(rownames(coords)) &&
             !identical(rownames(coords), as.character(seq_len(nrow(coords))))) {
    rownames(coords)
  } else {
    seq_len(n_sites)
  }
  row <- match(as.character(interactions$site), as.character(site_ids))
  if (anyNA(row)) {
    stop(sprintf("Sites in interactions not found in coords: %s",
                 paste(utils::head(unique(interactions$site[is.na(row)]), 5), collapse = ", ")),
         call. = FALSE)
  }

  tabulate_by <- function(key) {
    levels <- sort(unique(key))
    m <- matrix(0L, n_sites, length(levels), dimnames = list(NULL, levels))
    agg <- stats::aggregate(count, by = list(row = row, col = match(key, levels)), FUN = sum)
    m[cbind(agg$row, agg$col)] <- as.integer(agg$x)
    m
  }
  pair <- paste(interactions$lower, interactions$upper, sep = " -- ")
  list(interaction = tabulate_by(pair),
       lower = tabulate_by(as.character(interactions$lower)),
       upper = tabulate_by(as.character(interactions$upper)))
}


#' @export
print.spacc_network <- function(x, ...) {
  last <- function(m) mean(m[, x$n_sites])
  .print_curves(x, "spacc interaction network",
                sprintf("%d interactions between %d lower and %d upper species",
                        x$n_interactions, x$n_lower, x$n_upper),
                sprintf("  q = %s: mean final interaction diversity = %.1f", x$q,
                        vapply(x$interaction, last, numeric(1))),
                sprintf("Final coverage: interactions %.3f, lower %.3f, upper %.3f",
                        last(x$coverage$interaction), last(x$coverage$lower),
                        last(x$coverage$upper)),
                sprintf("Final connectance: %.3f", last(x$connectance)))
}


#' @export
summary.spacc_network <- function(object, ci_level = 0.95, ...) {
  hill <- do.call(rbind, lapply(c("interaction", "lower", "upper"), function(comp) {
    cbind(component = comp, .curve_frame(object[[comp]], "q", object$q, ci_level))
  }))
  extra <- list(coverage_interaction = object$coverage$interaction,
                coverage_lower = object$coverage$lower,
                coverage_upper = object$coverage$upper,
                connectance = object$connectance)
  other <- do.call(rbind, lapply(names(extra), function(comp) {
    cbind(component = comp, .curve_frame(list(extra[[comp]]), "q", NA_real_, ci_level))
  }))
  rbind(hill, other)
}


#' @export
as.data.frame.spacc_network <- function(x, row.names = NULL, optional = FALSE, ...) {
  summary(x)
}


#' @export
plot.spacc_network <- function(x, component = c("interaction", "species", "connectance",
                                                "coverage"),
                               q = 0, ci = TRUE, ci_alpha = 0.2, ...) {
  component <- match.arg(component)
  summ <- summary(x)
  if (component == "interaction") {
    df <- summ[summ$component == "interaction", ]
    df$label <- .q_label(df$q)
    return(.plot_curves_by(df, "label", "Order", xlab = "Sites accumulated",
                           ylab = "Interaction diversity",
                           title = "Spatial Accumulation of Interactions",
                           ci = ci, ci_alpha = ci_alpha))
  }
  if (component == "species") {
    if (!q %in% x$q) stop(sprintf("q = %g was not computed", q), call. = FALSE)
    df <- summ[summ$component %in% c("interaction", "lower", "upper") & summ$q == q, ]
    df$label <- factor(df$component, levels = c("interaction", "lower", "upper"),
                       labels = c("Interactions", "Lower level", "Upper level"))
    return(.plot_curves_by(df, "label", "Curve", xlab = "Sites accumulated",
                           ylab = sprintf("Hill number (q = %g)", q),
                           title = "Interactions and Species",
                           ci = ci, ci_alpha = ci_alpha))
  }
  if (component == "connectance") {
    df <- summ[summ$component == "connectance", ]
    return(.plot_curve(df, x = "sites", y = "mean", lower = "lower", upper = "upper",
                       xlab = "Sites accumulated", ylab = "Connectance",
                       title = "Connectance along Accumulation", ci = ci, ci_alpha = ci_alpha))
  }
  df <- summ[startsWith(summ$component, "coverage_"), ]
  df$label <- factor(sub("coverage_", "", df$component), levels = c("interaction", "lower", "upper"),
                     labels = c("Interactions", "Lower level", "Upper level"))
  .plot_curves_by(df, "label", "Curve", xlab = "Sites accumulated", ylab = "Sample coverage",
                  title = "Sample Coverage along Accumulation", ci = ci, ci_alpha = ci_alpha)
}
