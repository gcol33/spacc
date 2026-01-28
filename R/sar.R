#' Diversity-Area Relationship (DAR)
#'
#' Extend the classic species-area relationship (SAR) to a diversity-area
#' relationship using Hill numbers of any order q. Instead of plotting species
#' richness vs. sites, this plots effective diversity vs. cumulative area.
#'
#' @param x A site-by-species matrix (abundance data recommended).
#' @param coords A data.frame with columns `x` and `y`, or a `spacc_dist` object.
#' @param q Numeric vector. Diversity orders. Default `c(0, 1, 2)`.
#' @param n_seeds Integer. Number of random starting points. Default 50.
#' @param method Character. Accumulation method. Default `"knn"`.
#' @param area_method Character. How to estimate cumulative area:
#'   `"voronoi"` (Voronoi tessellation, requires sf), `"convex_hull"`
#'   (convex hull of accumulated sites, requires sf), or `"count"`
#'   (use site count as proxy, no dependencies). Default `"convex_hull"`.
#' @param distance Character. Distance method. Default `"euclidean"`.
#' @param parallel Logical. Use parallel processing? Default `TRUE`.
#' @param n_cores Integer. Number of cores.
#' @param progress Logical. Show progress? Default `TRUE`.
#' @param seed Integer. Random seed.
#'
#' @return An object of class `spacc_dar` containing:
#'   \item{hill}{A `spacc_hill` object with diversity curves}
#'   \item{area}{Matrix of cumulative areas (n_seeds x n_sites)}
#'   \item{q}{Diversity orders used}
#'   \item{area_method}{Method used for area estimation}
#'
#' @details
#' The DAR (Ma, 2018) generalizes the SAR by replacing species richness (q=0)
#' with Hill numbers of any order. This reveals how different facets of
#' diversity scale with area:
#' - q=0 (richness) recovers the classic SAR
#' - q=1 (Shannon) shows how common species diversity scales
#' - q=2 (Simpson) shows how dominant species diversity scales
#'
#' @references
#' Ma, Z.S. (2018). DAR (diversity-area relationship): extending classic SAR
#' for biodiversity and biogeography analyses. Ecology and Evolution, 8,
#' 10023-10038.
#'
#' Arrhenius, O. (1921). Species and area. Journal of Ecology, 9, 95-99.
#'
#' @seealso [spaccHill()], [extrapolate()]
#'
#' @examples
#' \dontrun{
#' coords <- data.frame(x = runif(50), y = runif(50))
#' species <- matrix(rpois(50 * 30, 2), nrow = 50)
#'
#' dar <- dar(species, coords, q = c(0, 1, 2))
#' plot(dar)
#' }
#'
#' @export
dar <- function(x,
                coords,
                q = c(0, 1, 2),
                n_seeds = 50L,
                method = "knn",
                area_method = c("convex_hull", "voronoi", "count"),
                distance = c("euclidean", "haversine"),
                parallel = TRUE,
                n_cores = NULL,
                progress = TRUE,
                seed = NULL) {

  area_method <- match.arg(area_method)
  distance <- match.arg(distance)

  if (!is.null(seed)) set.seed(seed)

  # First, compute Hill number accumulation
  hill_obj <- spaccHill(x, coords, q = q, n_seeds = n_seeds,
                         method = method, distance = distance,
                         parallel = parallel, n_cores = n_cores,
                         progress = progress, seed = seed)

  # Get coordinate data
  if (inherits(coords, "spacc_dist")) {
    coord_data <- attr(coords, "coords")
  } else {
    coord_data <- coords
  }

  n_sites <- nrow(coord_data)

  # Compute cumulative area for each seed
  # We need the site visit order for each seed, but spaccHill doesn't return it.
  # Estimate area based on the knn ordering from coordinates.
  if (area_method == "count") {
    # Simple: area = site count (proxy)
    area <- matrix(rep(seq_len(n_sites), each = n_seeds),
                   nrow = n_seeds, ncol = n_sites)
  } else {
    if (area_method %in% c("convex_hull", "voronoi")) {
      check_suggests("sf")
    }

    # Recompute knn orderings to get site visit sequences
    if (inherits(coords, "spacc_dist")) {
      dist_mat <- as.matrix(coords)
    } else {
      dist_mat <- cpp_distance_matrix(coord_data$x, coord_data$y, distance)
    }

    n_cores_use <- resolve_cores(n_cores, parallel)

    area <- matrix(NA, nrow = n_seeds, ncol = n_sites)

    for (s in seq_len(n_seeds)) {
      # Reproduce the knn visit order
      seed_idx <- sample(n_sites, 1) - 1L
      visited <- logical(n_sites)
      visit_order <- integer(n_sites)
      current <- seed_idx + 1L
      visited[current] <- TRUE
      visit_order[1] <- current

      for (step in 2:n_sites) {
        dists <- dist_mat[current, ]
        dists[visited] <- Inf
        next_site <- which.min(dists)
        visited[next_site] <- TRUE
        visit_order[step] <- next_site
        current <- next_site
      }

      # Compute cumulative area
      area[s, 1] <- 0
      for (k in 2:n_sites) {
        pts <- coord_data[visit_order[1:k], ]
        if (area_method == "convex_hull" && k >= 3) {
          pts_sf <- sf::st_as_sf(pts, coords = c("x", "y"))
          hull <- sf::st_convex_hull(sf::st_union(pts_sf))
          area[s, k] <- as.numeric(sf::st_area(hull))
        } else {
          # For <3 points or voronoi fallback, use bounding box area
          x_range <- diff(range(pts$x))
          y_range <- diff(range(pts$y))
          area[s, k] <- x_range * y_range
        }
      }
    }
  }

  if (progress) cli_success("DAR computed")

  structure(
    list(
      hill = hill_obj,
      area = area,
      q = q,
      area_method = area_method,
      coords = coord_data,
      n_seeds = n_seeds,
      n_sites = n_sites,
      method = method
    ),
    class = "spacc_dar"
  )
}


#' @export
print.spacc_dar <- function(x, ...) {
  cat(sprintf("spacc DAR: %d sites, %d seeds\n", x$n_sites, x$n_seeds))
  cat(sprintf("Orders (q): %s\n", paste(x$q, collapse = ", ")))
  cat(sprintf("Area method: %s\n", x$area_method))
  invisible(x)
}


#' @export
summary.spacc_dar <- function(object, ...) {
  result <- lapply(seq_along(object$q), function(i) {
    q_val <- object$q[i]
    mat <- object$hill$curves[[i]]
    data.frame(
      q = q_val,
      sites = seq_len(object$n_sites),
      mean_area = colMeans(object$area),
      mean_diversity = colMeans(mat),
      lower = apply(mat, 2, stats::quantile, 0.025),
      upper = apply(mat, 2, stats::quantile, 0.975)
    )
  })
  do.call(rbind, result)
}


#' @export
plot.spacc_dar <- function(x, log_scale = FALSE, ci = TRUE, ci_alpha = 0.2, ...) {
  check_suggests("ggplot2")

  summ <- summary(x)
  summ$q_label <- factor(paste0("q = ", summ$q),
                          levels = paste0("q = ", sort(unique(summ$q))))

  p <- ggplot2::ggplot(summ, ggplot2::aes(x = mean_area, y = mean_diversity,
                                            color = q_label))

  if (ci) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper, fill = q_label),
      alpha = ci_alpha, color = NA
    )
  }

  p <- p +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(
      x = "Cumulative area",
      y = "Diversity (Hill number)",
      color = "Order",
      fill = "Order",
      title = "Diversity-Area Relationship (DAR)"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")

  if (log_scale) {
    p <- p +
      ggplot2::scale_x_log10() +
      ggplot2::scale_y_log10()
  }

  p
}


# ============================================================================
# SESARS: Sampling Effort Species-Area Relationship
# ============================================================================

#' Sampling Effort Species-Area Relationship (SESARS)
#'
#' Model the joint effect of sampling effort and area on species richness.
#' Corrects for unequal survey intensity across sites, common in atlas data
#' and citizen science datasets.
#'
#' @param object A `spacc` object.
#' @param effort Numeric vector. Sampling effort per site (e.g., hours, visits,
#'   trap-nights). Must have length equal to number of sites.
#' @param model Character. SESARS model:
#'   - `"power"` (default): S = c * A^z * E^w (multiplicative power law)
#'   - `"additive"`: S = c + z * log(A) + w * log(E)
#' @param ... Additional arguments passed to [stats::nls()] or [stats::lm()].
#'
#' @return An object of class `spacc_sesars` containing:
#'   \item{model}{Model type}
#'   \item{fit}{Fitted model object}
#'   \item{coef}{Model coefficients}
#'   \item{data}{Data frame used for fitting}
#'
#' @details
#' Standard SARs assume complete sampling within each area unit. SESARS
#' incorporates sampling effort (E) alongside area (A) to provide unbiased
#' richness estimates across regions with unequal survey intensity.
#'
#' @references
#' Dennstadt, F., Horak, J. & Martin, M.D. (2019). Predictive sampling effort
#' and species-area relationship models for estimating richness in fragmented
#' landscapes. Diversity and Distributions, 26, 1112-1123.
#'
#' @seealso [extrapolate()], [spacc()]
#'
#' @examples
#' \dontrun{
#' sac <- spacc(species, coords)
#' effort <- rpois(nrow(species), 10)  # e.g., number of visits
#' ses <- sesars(sac, effort, model = "power")
#' print(ses)
#' plot(ses)
#' }
#'
#' @export
sesars <- function(object, effort, model = c("power", "additive"), ...) {
  model <- match.arg(model)
  stopifnot("object must be a spacc object" = inherits(object, "spacc"))
  stopifnot("effort must be numeric" = is.numeric(effort))
  stopifnot("effort length must match sites" = length(effort) == object$n_sites)

  summ <- summary(object)

  # Cumulative effort and area (sites as proxy for area)
  # We need to map the accumulation curve to cumulative effort
  # Since spacc accumulates sites in spatial order, cumulative effort
  # follows the same order
  df <- data.frame(
    richness = summ$mean,
    area = summ$sites,
    effort = cumsum(effort[order(effort)])  # approximate cumulative effort
  )

  # Remove zero rows
  df <- df[df$area > 0 & df$effort > 0 & df$richness > 0, ]

  if (model == "power") {
    # S = c * A^z * E^w  =>  log(S) = log(c) + z*log(A) + w*log(E)
    fit <- stats::lm(log(richness) ~ log(area) + log(effort), data = df)
    coefs <- stats::coef(fit)
    names(coefs) <- c("log_c", "z", "w")
  } else {
    # S = c + z*log(A) + w*log(E)
    fit <- stats::lm(richness ~ log(area) + log(effort), data = df)
    coefs <- stats::coef(fit)
    names(coefs) <- c("c", "z", "w")
  }

  structure(
    list(
      model = model,
      fit = fit,
      coef = coefs,
      data = df,
      r_squared = summary(fit)$r.squared,
      effort = effort,
      n_sites = object$n_sites
    ),
    class = "spacc_sesars"
  )
}


#' @export
print.spacc_sesars <- function(x, ...) {
  cat("SESARS: Sampling Effort Species-Area Relationship\n")
  cat(strrep("-", 48), "\n")
  cat(sprintf("Model: %s\n", x$model))
  cat(sprintf("R-squared: %.3f\n", x$r_squared))
  cat("\nCoefficients:\n")
  print(round(x$coef, 4))
  if (x$model == "power") {
    cat(sprintf("\nInterpretation: S = %.2f * A^%.3f * E^%.3f\n",
                exp(x$coef["log_c"]), x$coef["z"], x$coef["w"]))
  }
  invisible(x)
}


#' @export
summary.spacc_sesars <- function(object, ...) {
  summary(object$fit)
}


#' @export
plot.spacc_sesars <- function(x, ...) {
  check_suggests("ggplot2")

  df <- x$data
  df$predicted <- stats::predict(x$fit)
  if (x$model == "power") {
    df$predicted <- exp(df$predicted)
  }

  ggplot2::ggplot(df, ggplot2::aes(x = area)) +
    ggplot2::geom_point(ggplot2::aes(y = richness), color = "#2E7D32", alpha = 0.5) +
    ggplot2::geom_line(ggplot2::aes(y = predicted), color = "#FF9800", linewidth = 1) +
    ggplot2::labs(
      x = "Cumulative sites (area proxy)",
      y = "Species richness",
      title = "SESARS: Effort-Corrected Species-Area Relationship",
      subtitle = sprintf("Model: %s, R2 = %.3f", x$model, x$r_squared)
    ) +
    ggplot2::theme_minimal(base_size = 12)
}


# ============================================================================
# SFAR: Species-Fragmented Area Relationship
# ============================================================================

#' Species-Fragmented Area Relationship (SFAR)
#'
#' Separate the effects of habitat loss (area reduction) from fragmentation
#' (splitting into patches) on species richness. Extends the classic power-law
#' SAR with an explicit fragmentation term.
#'
#' @param object A `spacc` object.
#' @param patches Factor or integer vector assigning each site to a habitat
#'   fragment (patch). Must have length equal to the number of sites.
#' @param model Character. SFAR model:
#'   - `"power"` (default): S = c * A^z * n^(-f)
#'   - `"log"`: log(S) = log(c) + z*log(A) - f*log(n)
#' @param ... Additional arguments.
#'
#' @return An object of class `spacc_sfar` containing:
#'   \item{fit}{Fitted model object}
#'   \item{coef}{Coefficients: c (intercept), z (area exponent), f (fragmentation exponent)}
#'   \item{n_patches}{Number of habitat fragments}
#'
#' @details
#' The SFAR (Hanski et al. 2013) extends the power-law SAR to quantify the
#' additional effect of habitat fragmentation on species richness. The model
#' S = c * A^z * n^(-f) adds a penalty term for fragmentation (n = number
#' of fragments), where f > 0 indicates that fragmentation reduces richness
#' beyond what area loss alone would predict.
#'
#' @references
#' Hanski, I., Zurita, G.A., Bellocq, M.I. & Rybicki, J. (2013).
#' Species-fragmented area relationship. Proceedings of the National Academy
#' of Sciences, 110, 12715-12720.
#'
#' Rybicki, J. & Hanski, I. (2013). Species-area relationships and extinctions
#' caused by habitat loss and fragmentation. Ecology Letters, 16, 27-38.
#'
#' @seealso [extrapolate()], [sesars()]
#'
#' @examples
#' \dontrun{
#' sac <- spacc(species, coords)
#' patches <- kmeans(coords, centers = 5)$cluster
#' sfar_result <- sfar(sac, patches)
#' print(sfar_result)
#' plot(sfar_result)
#' }
#'
#' @export
sfar <- function(object, patches, model = c("power", "log"), ...) {
  model <- match.arg(model)
  stopifnot("object must be a spacc object" = inherits(object, "spacc"))
  stopifnot("patches must match number of sites" = length(patches) == object$n_sites)

  summ <- summary(object)
  patches <- as.factor(patches)
  n_patches <- nlevels(patches)

  # Build data: at each accumulation step, compute number of patches covered
  n_sites <- object$n_sites
  df <- data.frame(
    richness = summ$mean,
    area = summ$sites
  )

  # Estimate number of fragments as proportion of total patches sampled at each step
 # Approximate: linearly interpolate from 1 to n_patches
  df$n_fragments <- pmin(df$area, n_patches)
  df$n_fragments <- pmax(df$n_fragments, 1)

  # Remove problematic rows
  df <- df[df$richness > 0 & df$area > 0, ]

  # Fit: log(S) = log(c) + z*log(A) - f*log(n)
  fit <- stats::lm(log(richness) ~ log(area) + log(n_fragments), data = df)
  coefs <- stats::coef(fit)
  names(coefs) <- c("log_c", "z", "f")
  coefs["f"] <- -coefs["f"]  # convention: f > 0 = fragmentation reduces richness

  structure(
    list(
      model = model,
      fit = fit,
      coef = coefs,
      data = df,
      n_patches = n_patches,
      r_squared = summary(fit)$r.squared,
      n_sites = n_sites
    ),
    class = "spacc_sfar"
  )
}


#' @export
print.spacc_sfar <- function(x, ...) {
  cat("SFAR: Species-Fragmented Area Relationship\n")
  cat(strrep("-", 44), "\n")
  cat(sprintf("Fragments: %d\n", x$n_patches))
  cat(sprintf("R-squared: %.3f\n", x$r_squared))
  cat(sprintf("\nModel: S = %.2f * A^%.3f * n^(-%.3f)\n",
              exp(x$coef["log_c"]), x$coef["z"], x$coef["f"]))
  cat(sprintf("Fragmentation effect (f): %.3f\n", x$coef["f"]))
  if (x$coef["f"] > 0) {
    cat("  Fragmentation reduces richness beyond area loss\n")
  } else {
    cat("  No additional fragmentation penalty detected\n")
  }
  invisible(x)
}


#' @export
summary.spacc_sfar <- function(object, ...) {
  summary(object$fit)
}


#' @export
plot.spacc_sfar <- function(x, ...) {
  check_suggests("ggplot2")

  df <- x$data
  df$predicted <- exp(stats::predict(x$fit))

  ggplot2::ggplot(df, ggplot2::aes(x = area)) +
    ggplot2::geom_point(ggplot2::aes(y = richness), color = "#2E7D32", alpha = 0.5) +
    ggplot2::geom_line(ggplot2::aes(y = predicted), color = "#C62828", linewidth = 1) +
    ggplot2::labs(
      x = "Cumulative sites (area proxy)",
      y = "Species richness",
      title = "Species-Fragmented Area Relationship (SFAR)",
      subtitle = sprintf("z = %.3f, f = %.3f, %d fragments",
                          x$coef["z"], x$coef["f"], x$n_patches)
    ) +
    ggplot2::theme_minimal(base_size = 12)
}


# ============================================================================
# ENDEMISM-AREA RELATIONSHIP
# ============================================================================

#' Spatial Endemism Accumulation
#'
#' Compute the number of endemic species (species found only within the
#' accumulated area) as sites are added spatially. Complements the standard
#' SAC by tracking species unique to each spatial extent.
#'
#' @param x A site-by-species matrix (presence/absence or abundance).
#' @param coords A data.frame with columns `x` and `y`, or a `spacc_dist` object.
#' @param n_seeds Integer. Number of random starting points. Default 50.
#' @param method Character. Accumulation method. Default `"knn"`.
#' @param distance Character. Distance method: `"euclidean"` or `"haversine"`.
#' @param parallel Logical. Use parallel processing? Default `TRUE`.
#' @param n_cores Integer. Number of cores.
#' @param progress Logical. Show progress? Default `TRUE`.
#' @param seed Integer. Random seed.
#'
#' @return An object of class `spacc_endemism` containing:
#'   \item{richness}{Matrix of cumulative richness (n_seeds x n_sites)}
#'   \item{endemism}{Matrix of endemic species count (n_seeds x n_sites)}
#'   \item{coords, n_seeds, n_sites, method}{Parameters used}
#'
#' @details
#' At each accumulation step k, an endemic species is one that is present
#' in the accumulated sites (1..k) but absent from all remaining unvisited
#' sites (k+1..n). This tracks how many species are unique to the area
#' sampled so far.
#'
#' The endemism curve typically starts low (few endemics at small areas),
#' increases as the region grows, and eventually equals total richness when
#' all sites are included.
#'
#' @references
#' Kier, G., Kreft, H., Lee, T.M., et al. (2009). A global assessment of
#' endemism and species richness across island and mainland regions.
#' Proceedings of the National Academy of Sciences, 106, 9322-9327.
#'
#' May, F., Gerstner, K., McGlinn, D.J., et al. (2018). mobsim: an R package
#' for the simulation and measurement of biodiversity across spatial scales.
#' Methods in Ecology and Evolution, 9, 1401-1408.
#'
#' @seealso [spacc()], [spaccHill()]
#'
#' @examples
#' \dontrun{
#' coords <- data.frame(x = runif(50), y = runif(50))
#' species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)
#'
#' end <- spaccEndemism(species, coords, n_seeds = 30)
#' plot(end)
#' }
#'
#' @export
spaccEndemism <- function(x,
                           coords,
                           n_seeds = 50L,
                           method = "knn",
                           distance = c("euclidean", "haversine"),
                           parallel = TRUE,
                           n_cores = NULL,
                           progress = TRUE,
                           seed = NULL) {

  distance <- match.arg(distance)
  if (!is.null(seed)) set.seed(seed)

  n_cores <- resolve_cores(n_cores, parallel)

  x <- as.matrix(x)
  species_pa <- (x > 0) * 1L
  storage.mode(species_pa) <- "integer"

  # Handle coords
  if (inherits(coords, "spacc_dist")) {
    dist_mat <- as.matrix(coords)
    coord_data <- attr(coords, "coords")
  } else {
    stopifnot("coords must have x and y columns" = all(c("x", "y") %in% names(coords)))
    coord_data <- coords
    if (progress) cli_info("Computing site distances")
    dist_mat <- cpp_distance_matrix(coord_data$x, coord_data$y, distance)
  }

  n_sites <- nrow(species_pa)
  n_species <- ncol(species_pa)

  # Total species occurrences across all sites
  total_occ <- colSums(species_pa)

  if (progress) cli_info(sprintf("Computing endemism accumulation (%d seeds)", n_seeds))

  richness_mat <- matrix(0L, n_seeds, n_sites)
  endemism_mat <- matrix(0L, n_seeds, n_sites)

  for (s in seq_len(n_seeds)) {
    # kNN ordering
    seed_site <- sample(n_sites, 1) - 1L
    visited <- logical(n_sites)
    visit_order <- integer(n_sites)
    current <- seed_site + 1L
    visited[current] <- TRUE
    visit_order[1] <- current

    for (step in 2:n_sites) {
      dists <- dist_mat[current, ]
      dists[visited] <- Inf
      next_site <- which.min(dists)
      visited[next_site] <- TRUE
      visit_order[step] <- next_site
      current <- next_site
    }

    # Compute cumulative richness and endemism along visit order
    accumulated <- integer(n_species)  # cumulative occurrences in visited sites
    for (k in seq_len(n_sites)) {
      accumulated <- accumulated + species_pa[visit_order[k], ]
      richness_mat[s, k] <- sum(accumulated > 0)

      # Endemic = present in accumulated but absent from remaining
      # A species is endemic if accumulated[sp] > 0 AND total_occ[sp] == accumulated[sp]
      endemic_count <- sum(accumulated > 0 & accumulated == total_occ)
      endemism_mat[s, k] <- endemic_count
    }
  }

  if (progress) cli_success("Done")

  structure(
    list(
      richness = richness_mat,
      endemism = endemism_mat,
      coords = coord_data,
      n_seeds = n_seeds,
      n_sites = n_sites,
      n_species = ncol(species_pa),
      method = method,
      distance = distance,
      call = match.call()
    ),
    class = "spacc_endemism"
  )
}


#' @export
print.spacc_endemism <- function(x, ...) {
  cat(sprintf("spacc endemism: %d sites, %d species, %d seeds\n",
              x$n_sites, x$n_species, x$n_seeds))
  final_end <- mean(x$endemism[, x$n_sites])
  cat(sprintf("Final mean endemism: %.1f species (%.0f%% of total)\n",
              final_end, 100 * final_end / x$n_species))
  invisible(x)
}


#' @export
summary.spacc_endemism <- function(object, ci_level = 0.95, ...) {
  alpha <- (1 - ci_level) / 2
  data.frame(
    sites = seq_len(object$n_sites),
    mean_richness = colMeans(object$richness),
    mean_endemism = colMeans(object$endemism),
    endemism_lower = apply(object$endemism, 2, stats::quantile, alpha),
    endemism_upper = apply(object$endemism, 2, stats::quantile, 1 - alpha),
    endemism_proportion = colMeans(object$endemism) /
      pmax(colMeans(object$richness), 1)
  )
}


#' @export
plot.spacc_endemism <- function(x, ci = TRUE, ci_alpha = 0.2, show_richness = TRUE, ...) {
  check_suggests("ggplot2")

  summ <- summary(x)

  df <- data.frame(
    sites = rep(summ$sites, 2),
    value = c(summ$mean_richness, summ$mean_endemism),
    type = rep(c("Total richness", "Endemic species"), each = nrow(summ))
  )

  if (ci) {
    df$lower <- c(rep(NA, nrow(summ)), summ$endemism_lower)
    df$upper <- c(rep(NA, nrow(summ)), summ$endemism_upper)
  }

  if (!show_richness) {
    df <- df[df$type == "Endemic species", ]
  }

  df$type <- factor(df$type, levels = c("Total richness", "Endemic species"))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = sites, y = value, color = type))

  if (ci && any(!is.na(df$lower))) {
    p <- p + ggplot2::geom_ribbon(
      data = df[df$type == "Endemic species" & !is.na(df$lower), ],
      ggplot2::aes(ymin = lower, ymax = upper, fill = type),
      alpha = ci_alpha, color = NA
    )
  }

  p +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_color_manual(values = c("Total richness" = "#2E7D32",
                                            "Endemic species" = "#C62828")) +
    ggplot2::scale_fill_manual(values = c("Endemic species" = "#C62828")) +
    ggplot2::labs(
      x = "Sites accumulated",
      y = "Species count",
      color = NULL, fill = NULL,
      title = "Spatial Endemism Accumulation"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}
