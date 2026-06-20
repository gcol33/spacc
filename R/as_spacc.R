#' Convert an Accumulation Curve from Another Package to a spacc Object
#'
#' Bring an existing species accumulation curve into the spacc object system so
#' it can be summarised, plotted, and compared with [compare()] alongside
#' spatial curves. A \pkg{vegan} `specaccum` object is currently supported.
#'
#' @param x An object to convert. A `specaccum` object from
#'   `vegan::specaccum()` is supported.
#' @param ... Ignored, for method extensibility.
#'
#' @details
#' A `specaccum` object built with `method = "random"` stores one richness curve
#' per permutation in its `perm` matrix. Those curves become the rows of the
#' returned `curves` matrix and behave exactly like spacc seeds, so the
#' across-curve confidence band and the [compare()] permutation, bootstrap, and
#' AUC tests carry over unchanged. A `specaccum` object without a permutation
#' matrix (for example `method = "exact"` or `"collector"`) contributes its
#' single mean curve.
#'
#' Coordinates cannot be recovered from a `specaccum` object, so the result
#' carries no `coords`: it is a curve container for plotting and comparison, not
#' a re-runnable spatial object.
#'
#' @return An object of class `spacc` with `print()`, `summary()`, `plot()`,
#'   `as.data.frame()`, and [compare()] support.
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
    ". Supply a 'specaccum' object from vegan::specaccum().",
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
  storage.mode(curves) <- "double"

  n_sites <- ncol(curves)
  if (n_sites < 1L) {
    stop("specaccum object produced an empty curve.", call. = FALSE)
  }
  n_species <- as.integer(round(max(curves[, n_sites])))

  structure(
    list(
      curves = curves,
      coords = NULL,
      incidence_freq = NULL,
      n_seeds = nrow(curves),
      n_sites = n_sites,
      n_species = n_species,
      method = if (!is.null(x$method)) as.character(x$method) else "vegan",
      distance = NA_character_,
      backend = NA_character_,
      call = match.call()
    ),
    class = "spacc"
  )
}
