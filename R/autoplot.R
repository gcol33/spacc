#' Autoplot Methods for spacc Objects
#'
#' `ggplot2::autoplot()` methods for spacc objects. These are thin wrappers
#' around the corresponding `plot()` methods, provided for ggplot2 integration.
#'
#' @name autoplot.spacc
#' @param object A spacc object.
#' @param ... Additional arguments passed to the corresponding `plot()` method.
#' @return A ggplot2 object.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   coords <- data.frame(x = runif(30), y = runif(30))
#'   species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)
#'   sac <- spacc(species, coords, n_seeds = 5, progress = FALSE)
#'   ggplot2::autoplot(sac)
#' }
#' }
NULL


#' @rawNamespace S3method(ggplot2::autoplot, spacc)
autoplot.spacc <- function(object, ...) {
  check_suggests("ggplot2")
  plot(object, ...)
}


#' @rawNamespace S3method(ggplot2::autoplot, spacc_hill)
autoplot.spacc_hill <- function(object, ...) {
  check_suggests("ggplot2")
  plot(object, ...)
}


#' @rawNamespace S3method(ggplot2::autoplot, spacc_beta)
autoplot.spacc_beta <- function(object, ...) {
  check_suggests("ggplot2")
  plot(object, ...)
}


#' @rawNamespace S3method(ggplot2::autoplot, spacc_coverage)
autoplot.spacc_coverage <- function(object, ...) {
  check_suggests("ggplot2")
  plot(object, ...)
}


#' @rawNamespace S3method(ggplot2::autoplot, spacc_metrics)
autoplot.spacc_metrics <- function(object, ...) {
  check_suggests("ggplot2")
  plot(object, ...)
}


#' @rawNamespace S3method(ggplot2::autoplot, spacc_endemism)
autoplot.spacc_endemism <- function(object, ...) {
  check_suggests("ggplot2")
  plot(object, ...)
}


#' @rawNamespace S3method(ggplot2::autoplot, spacc_ses)
autoplot.spacc_ses <- function(object, ...) {
  check_suggests("ggplot2")
  plot(object, ...)
}


#' @rawNamespace S3method(ggplot2::autoplot, spacc_fit)
autoplot.spacc_fit <- function(object, ...) {
  check_suggests("ggplot2")
  plot(object, ...)
}


#' @rawNamespace S3method(ggplot2::autoplot, spacc_decay)
autoplot.spacc_decay <- function(object, ...) {
  check_suggests("ggplot2")
  plot(object, ...)
}


#' @rawNamespace S3method(ggplot2::autoplot, spacc_rare)
autoplot.spacc_rare <- function(object, ...) {
  check_suggests("ggplot2")
  plot(object, ...)
}


#' @rawNamespace S3method(ggplot2::autoplot, spacc_completeness)
autoplot.spacc_completeness <- function(object, ...) {
  check_suggests("ggplot2")
  plot(object, ...)
}


#' @rawNamespace S3method(ggplot2::autoplot, spacc_zeta)
autoplot.spacc_zeta <- function(object, ...) {
  check_suggests("ggplot2")
  plot(object, ...)
}


#' @rawNamespace S3method(ggplot2::autoplot, spacc_mem)
autoplot.spacc_mem <- function(object, ...) {
  check_suggests("ggplot2")
  plot(object, ...)
}


#' @rawNamespace S3method(ggplot2::autoplot, spacc_mem_partition)
autoplot.spacc_mem_partition <- function(object, ...) {
  check_suggests("ggplot2")
  plot(object, ...)
}


#' @rawNamespace S3method(ggplot2::autoplot, spacc_model_compare)
autoplot.spacc_model_compare <- function(object, ...) {
  check_suggests("ggplot2")
  plot(object, ...)
}
