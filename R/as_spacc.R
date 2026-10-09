#' Convert an Accumulation Curve from Another Package to a spacc Object
#'
#' Bring an existing species accumulation curve into the spacc object system so
#' it can be summarised, plotted, extrapolated and compared with [compare()]
#' alongside spatial curves.
#'
#' @param x An object to convert: a `specaccum` object from
#'   `vegan::specaccum()`, a numeric vector holding one accumulation curve, or
#'   a data frame with such a curve in one column (for example the output of
#'   `Rarefy::directionalSAC()` or `Rarefy::rare_alpha()`).
#' @param column For a data frame, the name of the column holding the curve.
#'   `NULL` (default) uses `N_SCR` (the spatially constrained curve of
#'   `Rarefy::directionalSAC()`) or `Rarefaction` (`Rarefy::rare_alpha()`),
#'   whichever is present.
#' @param ... Ignored, for method extensibility.
#'
#' @details
#' A `specaccum` object built with `method = "random"` stores one richness curve
#' per permutation in its `perm` matrix. Those curves become the rows of the
#' returned `curves` matrix and behave exactly like spacc seeds, so the
#' across-curve confidence band and the [compare()] permutation, bootstrap, and
#' AUC tests carry over unchanged. A `specaccum` object without a permutation
#' matrix (for example `method = "exact"` or `"collector"`), a numeric vector
#' and a data frame column each contribute a single mean curve, entry k being
#' the value after k sites.
#'
#' `Rarefy::directionalSAC(community, gradient)` with a site distance matrix
#' as `gradient` orders the sites by distance from each plot in turn and
#' averages the collector curves; this is `spacc(x, coords, method = "knn",
#' focal_points = coords)`, so the two curves agree up to the order of sites
#' at tied distances.
#'
#' Coordinates cannot be recovered from these objects, so the result carries
#' no `coords`: it is a curve container for plotting and comparison, not a
#' re-runnable spatial object.
#'
#' @return An object of class `spacc` with `print()`, `summary()`, `plot()`,
#'   `as.data.frame()`, and [compare()] support.
#'
#' @references
#' Thouverai, E., Pavoine, S., Tordoni, E., Chiarucci, A., Ricotta, C.,
#' Rocchini, D. & Bacaro, G. (2026). Rarefy: an R package for the calculation
#' of taxonomic, functional and phylogenetic sample-based spatially explicit
#' rarefaction curves. Ecological Informatics, 103803.
#' \doi{10.1016/j.ecoinf.2026.103803}
#'
#' @examples
#' \donttest{
#' if (requireNamespace("vegan", quietly = TRUE)) {
#'   data(BCI, package = "vegan")
#'   classic <- vegan::specaccum(BCI, method = "random")
#'   sp <- as_spacc(classic)
#'   summary(sp)
#'   plot(sp)
#' }
#' }
#' as_spacc(c(5, 8, 10, 11, 12))
#'
#' @seealso [compare()] to test a converted curve against a spatial one.
#'
#' @export
as_spacc <- function(x, ...) {
  UseMethod("as_spacc")
}


#' @rdname as_spacc
#' @export
as_spacc.default <- function(x, ...) {
  stop(
    "as_spacc() cannot convert an object of class ",
    paste(class(x), collapse = "/"),
    ". Supply a vegan::specaccum() object, a numeric curve or a data frame.",
    call. = FALSE
  )
}


#' @rdname as_spacc
#' @export
as_spacc.specaccum <- function(x, ...) {
  # method = "random" stores one curve per permutation (sites x perms);
  # other methods expose only the mean richness curve.
  if (!is.null(x$perm) && is.matrix(x$perm) && ncol(x$perm) > 0) {
    curves <- t(x$perm)
  } else if (!is.null(x$richness)) {
    curves <- matrix(x$richness, nrow = 1)
  } else {
    stop("specaccum object has neither a `perm` matrix nor a `richness` vector.",
         call. = FALSE)
  }
  .spacc_from_curves(curves, if (!is.null(x$method)) as.character(x$method) else "vegan",
                     match.call())
}


#' @rdname as_spacc
#' @export
as_spacc.numeric <- function(x, ...) {
  .spacc_from_curves(matrix(x, nrow = 1), "external", match.call())
}


#' @rdname as_spacc
#' @export
as_spacc.data.frame <- function(x, column = NULL, ...) {
  if (is.null(column)) {
    column <- intersect(c("N_SCR", "Rarefaction"), names(x))[1]
    if (is.na(column)) {
      stop("Name the column holding the curve with `column`.", call. = FALSE)
    }
  }
  if (!column %in% names(x)) stop(sprintf("No column '%s' in x", column), call. = FALSE)
  .spacc_from_curves(matrix(x[[column]], nrow = 1), "external", match.call())
}


# A coordinate-free spacc object around a matrix of curves (one per row).
.spacc_from_curves <- function(curves, method, call) {
  storage.mode(curves) <- "double"
  n_sites <- ncol(curves)
  if (n_sites < 1L || anyNA(curves)) {
    stop("The curve is empty or contains NA.", call. = FALSE)
  }
  structure(
    list(
      curves = curves,
      coords = NULL,
      incidence_freq = NULL,
      n_seeds = nrow(curves),
      n_sites = n_sites,
      n_species = as.integer(round(max(curves[, n_sites]))),
      method = method,
      distance = NA_character_,
      backend = NA_character_,
      call = call
    ),
    class = "spacc"
  )
}
