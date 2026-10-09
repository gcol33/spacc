#' Diversity-Area Relationship (DAR)
#'
#' Diversity against cumulative area along a spatial accumulation, for Hill
#' numbers (Ma 2018) or for any curve spacc computes: species richness,
#' phylogenetic and functional metrics, redundancy.
#'
#' @param x A site-by-species matrix (abundance data recommended), or a
#'   curve object: `spacc`, `spacc_hill`, `spacc_phylo`, `spacc_func` or
#'   `spacc_redundancy`. A matrix is first passed to [spaccHill()] with `q`.
#' @param coords A data.frame with columns `x` and `y`, or a `spacc_dist`
#'   object (matrix input only).
#' @param q Numeric vector. Diversity orders (matrix input only). Default
#'   `c(0, 1, 2)`.
#' @param n_seeds Integer. Number of random starting points. Default 50.
#' @param method Character. Accumulation method: `"knn"` or `"nn_walk"`.
#' @param area_method Character. How to measure the area covered by the first
#'   k sites of an order: `"convex_hull"` (area of their convex hull, 0 until
#'   three sites span a polygon), `"voronoi"` (summed area of their Voronoi
#'   cells within the bounding box of all sites) or `"count"` (number of
#'   sites). Areas are in squared coordinate units. Default `"convex_hull"`.
#' @param distance Character. Distance method. Default `"euclidean"`.
#' @param parallel Logical. Use parallel processing? Default `TRUE`.
#' @param n_cores Integer. Number of cores.
#' @param progress Logical. Show progress? Default `TRUE`.
#' @param seed Integer. Random seed.
#' @param focal_points Optional focal points for `method = "knn"`. See [spacc()].
#' @param focal_domain Optional polygonal focal domain. See [spacc()].
#' @param ... Arguments passed to methods.
#'
#' @return An object of class `spacc_dar` containing:
#'   \item{curves}{Named list of diversity curves (n_seeds x n_sites), one per
#'     metric or order}
#'   \item{metric}{Names of the curves}
#'   \item{area}{Matrix of cumulative areas (n_seeds x n_sites)}
#'   \item{power_law}{Data frame with the power-law fit \eqn{D = c A^z} of
#'     each mean curve (`metric`, `c`, `z`)}
#'   \item{q}{Diversity orders for Hill-number input, else `NULL`}
#'   \item{area_method}{Method used for area estimation}
#'   \item{source}{Class of the curve object the DAR was built from}
#'   \item{object}{That curve object}
#'
#' @details
#' The DAR (Ma, 2018) generalizes the SAR by replacing species richness (q=0)
#' with Hill numbers of any order. This reveals how different facets of
#' diversity scale with area:
#' - q=0 (richness) recovers the classic SAR
#' - q=1 (Shannon) shows how common species diversity scales
#' - q=2 (Simpson) shows how dominant species diversity scales
#'
#' The area is accumulated along the same site orders as the diversity curve,
#' so any curve object gives its own diversity-area relationship: functional
#' richness, dispersion and other facets can scale with area differently from
#' richness (Dias et al. 2025). The power-law exponent `z` is the slope of
#' log mean diversity against log mean area over the steps where both are
#' positive.
#'
#' @references
#' Ma, Z.S. (2018). DAR (diversity-area relationship): extending classic SAR
#' for biodiversity and biogeography analyses. Ecology and Evolution, 8,
#' 10023-10038.
#'
#' Arrhenius, O. (1921). Species and area. Journal of Ecology, 9, 95-99.
#'
#' Dias, R.A., Bastazini, V.A.G., Gianuca, A.T. & Matthews, T.J. (2025). The
#' spatial scaling of multiple dimensions of functional diversity in habitat
#' islands. Functional Ecology. \doi{10.1111/1365-2435.70158}
#'
#' @seealso [spaccHill()], [spaccFunc()], [spaccPhylo()], [extrapolate()]
#'
#' @examples
#' \donttest{
#' coords <- data.frame(x = runif(50), y = runif(50))
#' species <- matrix(rpois(50 * 30, 2), nrow = 50)
#'
#' d <- dar(species, coords, q = c(0, 1, 2))
#' plot(d)
#'
#' traits <- matrix(rnorm(30 * 2), nrow = 30)
#' func <- spaccFunc(species, coords, traits, metric = c("fdis", "fric"))
#' dar(func)
#' }
#'
#' @export
dar <- function(x, ...) UseMethod("dar")


#' @rdname dar
#' @export
dar.default <- function(x,
                        coords,
                        q = c(0, 1, 2),
                        n_seeds = 50L,
                        method = "knn",
                        area_method = c("convex_hull", "voronoi", "count"),
                        distance = c("euclidean", "haversine"),
                        parallel = TRUE,
                        n_cores = NULL,
                        progress = TRUE,
                        seed = NULL,
                        focal_points = NULL,
                        focal_domain = NULL,
                        ...) {
  area_method <- match.arg(area_method)
  distance <- match.arg(distance)
  hill_obj <- spaccHill(x, coords, q = q, n_seeds = n_seeds,
                        method = method, distance = distance,
                        parallel = parallel, n_cores = n_cores,
                        progress = progress, seed = seed,
                        focal_points = focal_points,
                        focal_domain = focal_domain)
  dar(hill_obj, area_method = area_method)
}


#' @rdname dar
#' @export
dar.spacc <- function(x, area_method = c("convex_hull", "voronoi", "count"), ...) {
  if (is_grouped(x)) stop("dar() needs an ungrouped spacc object", call. = FALSE)
  .dar_build(x, list(richness = x$curves), match.arg(area_method))
}


#' @rdname dar
#' @export
dar.spacc_hill <- function(x, area_method = c("convex_hull", "voronoi", "count"), ...) {
  .dar_build(x, x$curves, match.arg(area_method), q = x$q)
}


#' @rdname dar
#' @export
dar.spacc_phylo <- function(x, area_method = c("convex_hull", "voronoi", "count"), ...) {
  .dar_build(x, x$curves, match.arg(area_method))
}


#' @rdname dar
#' @export
dar.spacc_func <- function(x, area_method = c("convex_hull", "voronoi", "count"), ...) {
  .dar_build(x, x$curves, match.arg(area_method))
}


#' @rdname dar
#' @export
dar.spacc_redundancy <- function(x, area_method = c("convex_hull", "voronoi", "count"), ...) {
  .dar_build(x, x$curves, match.arg(area_method), q = x$q)
}


# Cumulative area of the first k sites along each order (n_orders x n_sites).
.cumulative_area <- function(coords, orders, area_method) {
  if (area_method == "count") {
    return(matrix(rep(seq_len(ncol(orders)), each = nrow(orders)), nrow = nrow(orders)))
  }
  cpp_area_along_orders(as.numeric(coords$x), as.numeric(coords$y), orders - 1L, area_method)
}


# Power-law fit D = c A^z to a mean diversity curve against mean area.
.power_law <- function(area, diversity) {
  keep <- is.finite(area) & is.finite(diversity) & area > 0 & diversity > 0
  if (sum(keep) < 3) return(c(c = NA_real_, z = NA_real_))
  fit <- stats::lm(log(diversity[keep]) ~ log(area[keep]))
  c(c = exp(unname(stats::coef(fit)[1])), z = unname(stats::coef(fit)[2]))
}


.dar_build <- function(obj, curves, area_method, q = NULL) {
  if (is.null(obj$orders) || is.null(obj$coords)) {
    stop("This object carries no site orders; recompute it with the current spacc",
         call. = FALSE)
  }
  area <- .cumulative_area(obj$coords, obj$orders, area_method)
  mean_area <- colMeans(area)
  fits <- lapply(curves, function(m) .power_law(mean_area, colMeans(m)))
  power_law <- data.frame(metric = names(curves),
                          c = vapply(fits, function(f) f[["c"]], numeric(1)),
                          z = vapply(fits, function(f) f[["z"]], numeric(1)),
                          row.names = NULL)
  structure(
    list(
      curves = curves,
      metric = names(curves),
      area = area,
      power_law = power_law,
      q = q,
      area_method = area_method,
      source = class(obj)[1],
      object = obj,
      coords = obj$coords,
      n_seeds = nrow(area),
      n_sites = ncol(area),
      method = obj$method
    ),
    class = "spacc_dar"
  )
}


#' @export
print.spacc_dar <- function(x, ...) {
  cat(sprintf("spacc DAR: %d sites, %d seeds (from %s)\n", x$n_sites, x$n_seeds, x$source))
  if (!is.null(x$q)) cat(sprintf("Orders (q): %s\n", paste(x$q, collapse = ", ")))
  cat(sprintf("Area method: %s\n", x$area_method))
  cat("Power law D = c * A^z:\n")
  for (i in seq_len(nrow(x$power_law))) {
    cat(sprintf("  %s: z = %.3f\n", x$power_law$metric[i], x$power_law$z[i]))
  }
  invisible(x)
}


#' @export
coef.spacc_dar <- function(object, ...) object$power_law


#' @export
summary.spacc_dar <- function(object, ci_level = 0.95, sd = FALSE, ...) {
  key <- if (!is.null(object$q)) "q" else "metric"
  values <- if (!is.null(object$q)) object$q else object$metric
  summ <- .curve_frame(object$curves, key, values, ci_level, sd = sd)
  names(summ)[names(summ) == "mean"] <- "mean_diversity"
  summ$mean_area <- rep(colMeans(object$area), length(object$curves))
  cols <- c(key, "sites", "mean_area", "mean_diversity", "lower", "upper", if (sd) "sd")
  summ[, cols]
}


#' @export
plot.spacc_dar <- function(x, log_scale = FALSE, ci = TRUE, ci_alpha = 0.2, ...) {
  summ <- summary(x)
  if (!is.null(x$q)) {
    summ$label <- .q_label(summ$q)
    legend <- "Order"
  } else {
    summ$label <- factor(summ$metric, levels = x$metric)
    legend <- "Metric"
  }
  if (log_scale) summ <- summ[summ$mean_area > 0 & summ$mean_diversity > 0, ]

  p <- .plot_curves_by(summ, "label", legend,
                       x = "mean_area", y = "mean_diversity",
                       xlab = "Cumulative area", ylab = "Diversity",
                       title = "Diversity-Area Relationship (DAR)",
                       ci = ci, ci_alpha = ci_alpha)

  if (log_scale) {
    p <- p + ggplot2::scale_x_log10() + ggplot2::scale_y_log10()
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
#' \donttest{
#' coords <- data.frame(x = runif(50), y = runif(50))
#' species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)
#' sac <- spacc(species, coords)
#' effort <- rpois(50, 10)
#' ses <- sesars(sac, effort, model = "power")
#' print(ses)
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
    spacc_theme()
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
#' \donttest{
#' coords <- data.frame(x = runif(50), y = runif(50))
#' species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)
#' sac <- spacc(species, coords)
#' patches <- kmeans(coords, centers = 5)$cluster
#' sfar_result <- sfar(sac, patches)
#' print(sfar_result)
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
    spacc_theme()
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
#' @param map Logical. If `TRUE`, compute per-site endemism by running
#'   accumulation from each site as seed. Stores a `site_values` data.frame
#'   enabling `plot(type = "map")` and [as_sf()]. Default `FALSE`.
#' @param progress Logical. Show progress? Default `TRUE`.
#' @param seed Integer. Random seed.
#' @param focal_points Optional focal points for `method = "knn"`. See [spacc()].
#' @param focal_domain Optional polygonal focal domain. See [spacc()].
#'
#' @return An object of class `spacc_endemism` containing:
#'   \item{richness}{Matrix of cumulative richness (n_seeds x n_sites)}
#'   \item{endemism}{Matrix of endemic species count (n_seeds x n_sites)}
#'   \item{site_values}{Per-site endemism data.frame (if `map = TRUE`)}
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
#' \donttest{
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
                           method = c("knn", "nn_walk"),
                           distance = c("euclidean", "haversine"),
                           map = FALSE,
                           parallel = TRUE,
                           n_cores = NULL,
                           progress = TRUE,
                           seed = NULL,
                           focal_points = NULL,
                           focal_domain = NULL) {

  method <- match.arg(method)
  distance <- match.arg(distance)
  if (!is.null(seed)) set.seed(seed)

  n_cores <- resolve_cores(n_cores, parallel)

  x <- as.matrix(x)
  species_pa <- (x > 0) * 1L
  storage.mode(species_pa) <- "integer"

  site <- .resolve_coords(coords)
  coord_data <- site$coords
  stopifnot("x and coords must have same rows" = nrow(x) == nrow(coord_data))

  n_sites <- nrow(species_pa)
  n_species <- ncol(species_pa)
  ordering <- .accumulation_orders(method, coord_data, n_seeds, distance,
                                   site$dist_mat, focal_points, focal_domain,
                                   n_cores = n_cores)
  n_seeds <- ordering$n_seeds

  # A species is endemic to the accumulated area once all of its occurrences
  # have been visited.
  total_occ <- colSums(species_pa)
  endemism_along <- function(visit_order) {
    accumulated <- integer(n_species)
    out <- matrix(0L, 2L, n_sites)
    for (k in seq_len(n_sites)) {
      accumulated <- accumulated + species_pa[visit_order[k], ]
      out[1L, k] <- sum(accumulated > 0)
      out[2L, k] <- sum(accumulated > 0 & accumulated == total_occ)
    }
    out
  }

  if (progress) cli_info(sprintf("Computing endemism accumulation (%d seeds)", n_seeds))

  richness_mat <- matrix(0L, n_seeds, n_sites)
  endemism_mat <- matrix(0L, n_seeds, n_sites)
  for (s in seq_len(n_seeds)) {
    along <- endemism_along(ordering$orders[s, ])
    richness_mat[s, ] <- along[1L, ]
    endemism_mat[s, ] <- along[2L, ]
  }

  site_values <- NULL
  if (map) {
    if (progress) cli_info("Computing per-site endemism map")
    site_values <- .site_map_values(method, coord_data, distance, site$dist_mat, n_cores,
      function(orders) {
        site_endemism <- vapply(seq_len(nrow(orders)), function(s) {
          endemism_along(orders[s, ])[2L, n_sites]
        }, numeric(1))
        list(endemism = site_endemism,
             endemism_prop = site_endemism / sum(total_occ > 0))
      })
  }

  if (progress) cli_success("Done")

  structure(
    list(
      richness = richness_mat,
      endemism = endemism_mat,
      site_values = site_values,
      coords = coord_data,
      n_seeds = n_seeds,
      n_sites = n_sites,
      n_species = ncol(species_pa),
      method = method,
      distance = distance,
      focal_points = ordering$focal_points,
      call = match.call()
    ),
    class = "spacc_endemism"
  )
}


#' @export
print.spacc_endemism <- function(x, ...) {
  final_end <- mean(x$endemism[, x$n_sites])
  .print_curves(x, "spacc endemism",
                sprintf("Final mean endemism: %.1f species (%.0f%% of total)",
                        final_end, 100 * final_end / x$n_species))
}


#' @export
summary.spacc_endemism <- function(object, ci_level = 0.95, ...) {
  endemism <- .curve_stats(object$endemism, ci_level)
  mean_richness <- colMeans(object$richness)
  data.frame(
    sites = seq_len(object$n_sites),
    mean_richness = mean_richness,
    mean_endemism = endemism$mean,
    endemism_lower = endemism$lower,
    endemism_upper = endemism$upper,
    endemism_proportion = endemism$mean / pmax(mean_richness, 1)
  )
}


#' @export
plot.spacc_endemism <- function(x, type = c("curve", "map"),
                                 ci = TRUE, ci_alpha = 0.2, show_richness = TRUE,
                                 metric = c("endemism", "endemism_prop"),
                                 point_size = 3, palette = "viridis", ...) {
  type <- match.arg(type)

  if (type == "map") {
    metric <- match.arg(metric)
    return(.plot_site_map(x, metric,
                          title = sprintf("Per-site %s", gsub("_", " ", metric)),
                          subtitle = sprintf("%d sites", x$n_sites),
                          point_size = point_size, palette = palette))
  }

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
    ggplot2::scale_fill_manual(values = c("Endemic species" = "#C62828"),
                               guide = "none") +
    ggplot2::labs(
      x = "Sites accumulated",
      y = "Species count",
      color = NULL, fill = NULL,
      title = "Spatial Endemism Accumulation"
    ) +
    spacc_theme() +
    ggplot2::theme(legend.position = "bottom")
}


#' @export
as_sf.spacc_endemism <- function(x, crs = NULL) .as_sf_site_values(x, crs)
