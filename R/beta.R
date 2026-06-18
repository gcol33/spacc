#' Spatial Beta Diversity Accumulation
#'
#' Analyze how beta diversity changes as sites are accumulated spatially.
#' Partitions beta diversity into turnover (species replacement) and
#' nestedness (species loss) components following Baselga (2010).
#'
#' @param x A site-by-species matrix (presence/absence or abundance).
#' @param coords A data.frame with columns `x` and `y`, or a `spacc_dist` object.
#' @param traits Optional species-by-traits matrix or data.frame (row names
#'   matching species). When supplied, functional beta diversity is computed.
#' @param tree Optional phylogenetic tree of class `phylo`, or a pairwise
#'   phylogenetic distance matrix. When supplied, phylogenetic beta diversity is
#'   computed. Supply at most one of `traits` or `tree`.
#' @param n_seeds Integer. Number of random starting points. Default 50.
#' @param method Character. Accumulation method. Default `"knn"`.
#' @param index Character. Dissimilarity index: `"sorensen"` (default) or `"jaccard"`.
#' @param distance Character. Distance method: `"euclidean"` or `"haversine"`.
#' @param parallel Logical. Use parallel processing? Default `TRUE`.
#' @param n_cores Integer. Number of cores.
#' @param progress Logical. Show progress? Default `TRUE`.
#' @param seed Integer. Random seed.
#' @param map Logical. If `TRUE`, run accumulation from every site as seed
#'   and store per-site final beta values for spatial mapping. Enables
#'   [as_sf()] and `plot(type = "map")`. Default `FALSE`.
#'
#' @return An object of class `spacc_beta` containing:
#'   \item{beta_total}{Matrix of total beta diversity (n_seeds x n_sites-1)}
#'   \item{beta_turnover}{Matrix of turnover component}
#'   \item{beta_nestedness}{Matrix of nestedness component}
#'   \item{distance}{Matrix of cumulative distances}
#'   \item{n_seeds, n_sites, method, index}{Parameters used}
#'
#' @details
#' At each step of spatial accumulation, beta diversity is calculated between
#' the accumulated species pool and the newly added site. This reveals how
#' species composition changes as you expand spatially.
#'
#' **Interpretation:**
#' - High turnover: New sites contribute different species (replacement)
#' - High nestedness: New sites contribute subsets of existing species (loss)
#'
#' The sum of turnover and nestedness equals total beta diversity.
#'
#' Supplying `traits` computes functional beta diversity, weighting the partition
#' by the trait distinctiveness of the species exchanged between the accumulated
#' pool and each new site. Supplying `tree` computes phylogenetic beta diversity
#' (PhyloSor), weighting by shared branch length. The taxonomic, functional, and
#' phylogenetic variants all return a `spacc_beta` object whose `beta_type` field
#' records which was computed. `map = TRUE` is supported for the taxonomic
#' variant only.
#'
#' @references
#' Baselga, A. (2010). Partitioning the turnover and nestedness components
#' of beta diversity. Global Ecology and Biogeography, 19, 134-143.
#'
#' @seealso [betapart::beta.pair()] for pairwise beta diversity
#'
#' @examples
#' \donttest{
#' coords <- data.frame(x = runif(50), y = runif(50))
#' species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)
#'
#' beta <- spaccBeta(species, coords, n_seeds = 30)
#' plot(beta)
#'
#' # Compare Sorensen vs Jaccard
#' beta_jac <- spaccBeta(species, coords, index = "jaccard")
#'
#' # Functional beta diversity (pass a species-by-traits matrix)
#' traits <- matrix(rnorm(30 * 3), nrow = 30)
#' rownames(traits) <- colnames(species) <- paste0("sp", 1:30)
#' beta_func <- spaccBeta(species, coords, traits = traits)
#' }
#'
#' @export
spaccBeta <- function(x,
                      coords,
                      traits = NULL,
                      tree = NULL,
                      n_seeds = 50L,
                      method = "knn",
                      index = c("sorensen", "jaccard"),
                      distance = c("euclidean", "haversine"),
                      parallel = TRUE,
                      n_cores = NULL,
                      progress = TRUE,
                      seed = NULL,
                      map = FALSE) {

  index <- match.arg(index)
  distance <- match.arg(distance)

  if (!is.null(traits) && !is.null(tree)) {
    stop("Supply either `traits` (functional beta) or `tree` (phylogenetic beta), not both.")
  }
  if (map && (!is.null(traits) || !is.null(tree))) {
    stop("`map = TRUE` is supported for taxonomic beta only (no `traits`/`tree`).")
  }

  if (!is.null(seed)) set.seed(seed)

  n_cores <- resolve_cores(n_cores, parallel)

  x <- as.matrix(x)
  species_pa <- (x > 0) * 1L
  storage.mode(species_pa) <- "integer"

  # Handle coords -> distance matrix (shared across variants)
  if (inherits(coords, "spacc_dist")) {
    dist_mat <- as.matrix(coords)
    coord_data <- attr(coords, "coords")
  } else {
    stopifnot("coords must have x and y columns" = all(c("x", "y") %in% names(coords)))
    coord_data <- coords
    if (progress) cli_info("Computing site distances")
    dist_mat <- cpp_distance_matrix(coord_data$x, coord_data$y, distance)
  }

  stopifnot("x and coords must have same rows" = nrow(x) == nrow(coord_data))

  n_sites <- nrow(x)
  n_species <- ncol(x)

  # Dispatch on the requested dimension
  if (!is.null(tree)) {
    comp <- beta_phylogenetic(x, species_pa, tree, dist_mat, n_seeds, index, progress)
    beta_type <- "phylogenetic"
  } else if (!is.null(traits)) {
    comp <- beta_functional(x, species_pa, traits, dist_mat, n_seeds, index, progress)
    beta_type <- "functional"
  } else {
    comp <- beta_taxonomic(species_pa, dist_mat, n_seeds, index, n_cores,
                           progress, map, coord_data, n_sites)
    beta_type <- "taxonomic"
  }

  if (progress) cli_success("Done")

  structure(
    list(
      beta_total = comp$beta_total,
      beta_turnover = comp$beta_turnover,
      beta_nestedness = comp$beta_nestedness,
      distance = comp$distance,
      coords = coord_data,
      site_values = comp$site_values,
      n_seeds = n_seeds,
      n_sites = n_sites,
      n_species = n_species,
      method = method,
      index = index,
      beta_type = beta_type,
      call = match.call()
    ),
    class = "spacc_beta"
  )
}


#' @export
print.spacc_beta <- function(x, ...) {
  cat(sprintf("spacc beta diversity: %d sites, %d seeds\n",
              x$n_sites, x$n_seeds))
  cat(sprintf("Index: %s, Method: %s\n", x$index, x$method))

  # Final values
  final_total <- mean(x$beta_total[, ncol(x$beta_total)])
  final_turn <- mean(x$beta_turnover[, ncol(x$beta_turnover)])
  final_nest <- mean(x$beta_nestedness[, ncol(x$beta_nestedness)])

 cat(sprintf("Mean final beta: %.3f (turnover: %.3f, nestedness: %.3f)\n",
              final_total, final_turn, final_nest))
  invisible(x)
}


#' @export
summary.spacc_beta <- function(object, ci_level = 0.95, ...) {
  alpha <- (1 - ci_level) / 2
  n_steps <- ncol(object$beta_total)

  data.frame(
    step = 1:n_steps,
    mean_distance = colMeans(object$distance),
    beta_total = colMeans(object$beta_total),
    beta_total_lower = apply(object$beta_total, 2, stats::quantile, alpha),
    beta_total_upper = apply(object$beta_total, 2, stats::quantile, 1 - alpha),
    beta_turnover = colMeans(object$beta_turnover),
    beta_nestedness = colMeans(object$beta_nestedness)
  )
}


#' @export
plot.spacc_beta <- function(x, type = c("curve", "map"), partition = TRUE,
                            xaxis = c("sites", "distance"),
                            ci = TRUE, ci_alpha = 0.2,
                            component = c("beta_total", "beta_turnover", "beta_nestedness"),
                            point_size = 3, palette = "viridis", ...) {
  type <- match.arg(type)

  if (type == "map") {
    if (is.null(x$site_values)) stop("No map data. Rerun spaccBeta() with map = TRUE.")
    component <- match.arg(component)
    return(plot_spatial_map(x$site_values, component,
                            title = sprintf("Beta diversity map: %s", component),
                            subtitle = sprintf("%d sites, %s (%s)", x$n_sites, x$index, x$method),
                            point_size = point_size, palette = palette))
  }

  check_suggests("ggplot2")

  xaxis <- match.arg(xaxis)
  summ <- summary(x)

  if (xaxis == "sites") {
    summ$x <- summ$step + 1  # step 1 = site 2
    xlab <- "Sites accumulated"
  } else {
    summ$x <- summ$mean_distance
    xlab <- "Cumulative distance"
  }

  if (partition) {
    # Reshape for stacked/faceted plot
    df <- data.frame(
      x = rep(summ$x, 3),
      value = c(summ$beta_total, summ$beta_turnover, summ$beta_nestedness),
      component = rep(c("Total", "Turnover", "Nestedness"), each = nrow(summ))
    )
    df$component <- factor(df$component, levels = c("Total", "Turnover", "Nestedness"))

    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = value, color = component)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::scale_color_manual(values = c("Total" = "#2E7D32",
                                              "Turnover" = "#1565C0",
                                              "Nestedness" = "#C62828"))
  } else {
    p <- ggplot2::ggplot(summ, ggplot2::aes(x = x, y = beta_total))

    if (ci) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = beta_total_lower, ymax = beta_total_upper),
        alpha = ci_alpha, fill = "#4CAF50"
      )
    }

    p <- p + ggplot2::geom_line(linewidth = 1, color = "#2E7D32")
  }

  p +
    ggplot2::labs(
      x = xlab,
      y = paste0("Beta diversity (", x$index, ")"),
      color = "Component",
      title = "Spatial Beta Diversity Accumulation"
    ) +
    spacc_theme() +
    ggplot2::theme(legend.position = "bottom")
}


#' @export
as_sf.spacc_beta <- function(x, crs = NULL) {
  if (is.null(x$site_values)) stop("No map data. Rerun spaccBeta() with map = TRUE.")
  as_sf_from_df(x$site_values, crs = crs)
}


# ============================================================================
# INTERNAL COMPUTE HELPERS (taxonomic / functional / phylogenetic)
# ============================================================================

# Taxonomic beta via the C++ backend; optionally per-site map values.
beta_taxonomic <- function(species_pa, dist_mat, n_seeds, index, n_cores,
                           progress, map, coord_data, n_sites) {
  use_jaccard <- index == "jaccard"

  if (progress) cli_info(sprintf("Computing beta diversity (%s, %d seeds)",
                                  index, n_seeds))

  result <- cpp_beta_knn_parallel(species_pa, dist_mat, n_seeds,
                                   use_jaccard, n_cores, progress)

  site_values <- NULL
  if (map) {
    if (progress) cli_info("Computing per-site beta map values (all sites as seeds)")
    map_result <- cpp_beta_knn_parallel(species_pa, dist_mat, n_sites,
                                         use_jaccard, n_cores, progress)
    n_steps <- ncol(map_result$beta_total)
    site_values <- data.frame(
      site_id = seq_len(n_sites),
      x = coord_data$x,
      y = coord_data$y,
      beta_total = map_result$beta_total[, n_steps],
      beta_turnover = map_result$beta_turnover[, n_steps],
      beta_nestedness = map_result$beta_nestedness[, n_steps]
    )
    if (progress) cli_success("Map values computed")
  }

  list(
    beta_total = result$beta_total,
    beta_turnover = result$beta_turnover,
    beta_nestedness = result$beta_nestedness,
    distance = result$distance,
    site_values = site_values
  )
}


# Functional beta: Baselga partition weighted by trait distinctiveness of the
# species exchanged between the accumulated pool and each new site.
beta_functional <- function(x, species_pa, traits, dist_mat, n_seeds, index,
                            progress) {
  traits <- as.matrix(traits)
  n_sites <- nrow(species_pa)
  n_species <- ncol(species_pa)

  # Match traits to species
  if (!is.null(colnames(x)) && !is.null(rownames(traits))) {
    if (!all(colnames(x) %in% rownames(traits))) {
      stop("Some species in x are not in traits")
    }
    traits <- traits[colnames(x), , drop = FALSE]
  }
  stopifnot("Traits must have one row per species" = nrow(traits) == n_species)

  # Compute pairwise trait distances between species
  trait_dist <- as.matrix(stats::dist(traits))

  if (progress) cli_info(sprintf("Computing functional beta diversity (%s, %d seeds)",
                                  index, n_seeds))

  use_jaccard <- index == "jaccard"
  n_steps <- n_sites - 1

  beta_total <- matrix(0, n_seeds, n_steps)
  beta_turn <- matrix(0, n_seeds, n_steps)
  beta_nest <- matrix(0, n_seeds, n_steps)
  distances <- matrix(0, n_seeds, n_steps)

  for (s in seq_len(n_seeds)) {
    seed_site <- sample(n_sites, 1) - 1L
    visited <- logical(n_sites)
    current <- seed_site + 1L
    visited[current] <- TRUE

    accumulated_sp <- which(species_pa[current, ] > 0)

    for (step in seq_len(n_steps)) {
      dists <- dist_mat[current, ]
      dists[visited] <- Inf
      next_site <- which.min(dists)
      distances[s, step] <- dists[next_site]
      visited[next_site] <- TRUE

      new_sp <- which(species_pa[next_site, ] > 0)

      # Functional beta: based on mean trait distance between shared/unique species
      shared <- intersect(accumulated_sp, new_sp)
      only_acc <- setdiff(accumulated_sp, new_sp)
      only_new <- setdiff(new_sp, accumulated_sp)

      a <- length(shared)
      b <- length(only_acc)
      c_count <- length(only_new)

      if (a + b + c_count > 0) {
        if (use_jaccard) {
          beta_total[s, step] <- (b + c_count) / (a + b + c_count)
          min_bc <- min(b, c_count)
          if (a + min_bc > 0) {
            beta_turn[s, step] <- 2 * min_bc / (a + 2 * min_bc)
          }
        } else {
          beta_total[s, step] <- (b + c_count) / (2 * a + b + c_count)
          min_bc <- min(b, c_count)
          if (a + min_bc > 0) {
            beta_turn[s, step] <- min_bc / (a + min_bc)
          }
        }
        beta_nest[s, step] <- beta_total[s, step] - beta_turn[s, step]

        # Weight by mean trait distance of unique species
        if (length(c(only_acc, only_new)) > 0 && length(shared) > 0) {
          unique_sp <- c(only_acc, only_new)
          mean_trait_dist <- mean(trait_dist[unique_sp, shared])
          # Scale beta by trait distinctiveness (normalized to 0-1)
          max_dist <- max(trait_dist)
          if (max_dist > 0) {
            trait_weight <- mean_trait_dist / max_dist
            beta_total[s, step] <- beta_total[s, step] * trait_weight
            beta_turn[s, step] <- beta_turn[s, step] * trait_weight
            beta_nest[s, step] <- beta_total[s, step] - beta_turn[s, step]
          }
        }
      }

      accumulated_sp <- union(accumulated_sp, new_sp)
      current <- next_site
    }
  }

  distances <- t(apply(distances, 1, cumsum))

  list(
    beta_total = beta_total,
    beta_turnover = beta_turn,
    beta_nestedness = beta_nest,
    distance = distances,
    site_values = NULL
  )
}


# Phylogenetic beta: PhyloSor partition weighted by shared branch length.
beta_phylogenetic <- function(x, species_pa, tree, dist_mat, n_seeds, index,
                              progress) {
  n_sites <- nrow(species_pa)

  # Get phylogenetic distances
  if (inherits(tree, "phylo")) {
    check_suggests("ape")
    phylo_dist <- as.matrix(ape::cophenetic.phylo(tree))
    if (!is.null(colnames(x))) {
      if (!all(colnames(x) %in% tree$tip.label)) {
        stop("Some species in x are not in the tree")
      }
      phylo_dist <- phylo_dist[colnames(x), colnames(x)]
    }
  } else if (is.matrix(tree)) {
    phylo_dist <- tree
  } else {
    stop("tree must be a phylo object or distance matrix")
  }

  if (progress) cli_info(sprintf("Computing phylogenetic beta diversity (%s, %d seeds)",
                                  index, n_seeds))

  use_jaccard <- index == "jaccard"
  n_steps <- n_sites - 1

  beta_total <- matrix(0, n_seeds, n_steps)
  beta_turn <- matrix(0, n_seeds, n_steps)
  beta_nest <- matrix(0, n_seeds, n_steps)
  distances <- matrix(0, n_seeds, n_steps)

  for (s in seq_len(n_seeds)) {
    seed_site <- sample(n_sites, 1) - 1L
    visited <- logical(n_sites)
    current <- seed_site + 1L
    visited[current] <- TRUE

    accumulated_sp <- which(species_pa[current, ] > 0)

    for (step in seq_len(n_steps)) {
      dists <- dist_mat[current, ]
      dists[visited] <- Inf
      next_site <- which.min(dists)
      distances[s, step] <- dists[next_site]
      visited[next_site] <- TRUE

      new_sp <- which(species_pa[next_site, ] > 0)

      shared <- intersect(accumulated_sp, new_sp)
      only_acc <- setdiff(accumulated_sp, new_sp)
      only_new <- setdiff(new_sp, accumulated_sp)

      # Phylogenetic beta: weight a, b, c by mean phylogenetic distances
      if (length(c(accumulated_sp, new_sp)) >= 2) {
        # Compute phylogenetically weighted components
        a_phylo <- if (length(shared) >= 2) sum(phylo_dist[shared, shared]) / 2 else 0
        b_phylo <- if (length(only_acc) > 0 && length(shared) > 0) {
          sum(phylo_dist[only_acc, shared])
        } else 0
        c_phylo <- if (length(only_new) > 0 && length(shared) > 0) {
          sum(phylo_dist[only_new, shared])
        } else 0

        total_phylo <- a_phylo + b_phylo + c_phylo
        if (total_phylo > 0) {
          if (use_jaccard) {
            beta_total[s, step] <- (b_phylo + c_phylo) / total_phylo
            min_bc <- min(b_phylo, c_phylo)
            if (a_phylo + min_bc > 0) {
              beta_turn[s, step] <- 2 * min_bc / (a_phylo + 2 * min_bc)
            }
          } else {
            beta_total[s, step] <- (b_phylo + c_phylo) / (2 * a_phylo + b_phylo + c_phylo)
            min_bc <- min(b_phylo, c_phylo)
            if (a_phylo + min_bc > 0) {
              beta_turn[s, step] <- min_bc / (a_phylo + min_bc)
            }
          }
          beta_nest[s, step] <- beta_total[s, step] - beta_turn[s, step]
        }
      }

      accumulated_sp <- union(accumulated_sp, new_sp)
      current <- next_site
    }
  }

  distances <- t(apply(distances, 1, cumsum))

  list(
    beta_total = beta_total,
    beta_turnover = beta_turn,
    beta_nestedness = beta_nest,
    distance = distances,
    site_values = NULL
  )
}


# ============================================================================
# DEPRECATED: superseded by spaccBeta(traits = ) / spaccBeta(tree = )
# ============================================================================

#' @rdname spaccBeta
#' @export
spaccBetaFunc <- function(x, coords, traits, n_seeds = 50L, method = "knn",
                          index = c("sorensen", "jaccard"),
                          distance = c("euclidean", "haversine"),
                          parallel = TRUE, n_cores = NULL, progress = TRUE,
                          seed = NULL) {
  .Deprecated("spaccBeta(traits = ...)")
  spaccBeta(x, coords, traits = traits, n_seeds = n_seeds, method = method,
            index = index, distance = distance, parallel = parallel,
            n_cores = n_cores, progress = progress, seed = seed)
}


#' @rdname spaccBeta
#' @export
spaccBetaPhylo <- function(x, coords, tree, n_seeds = 50L, method = "knn",
                           index = c("sorensen", "jaccard"),
                           distance = c("euclidean", "haversine"),
                           parallel = TRUE, n_cores = NULL, progress = TRUE,
                           seed = NULL) {
  .Deprecated("spaccBeta(tree = ...)")
  spaccBeta(x, coords, tree = tree, n_seeds = n_seeds, method = method,
            index = index, distance = distance, parallel = parallel,
            n_cores = n_cores, progress = progress, seed = seed)
}
