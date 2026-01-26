#' @keywords internal
"_PACKAGE"

#' @useDynLib spacc, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom stats quantile sd coef confint predict nls nls.control AIC dist
#' @importFrom parallel detectCores
NULL


# Check for suggested packages
check_suggests <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      sprintf("Package '%s' is required for this function. Install with: install.packages('%s')", pkg, pkg),
      call. = FALSE
    )
  }
}


# Safe cli messaging (falls back to cat if cli not available)
cli_info <- function(msg) {
  if (requireNamespace("cli", quietly = TRUE)) {
    cli::cli_alert_info(msg)
  } else {
    cat("i", msg, "\n")
  }
}

cli_success <- function(msg) {
  if (requireNamespace("cli", quietly = TRUE)) {
    cli::cli_alert_success(msg)
  } else {
    cat("v", msg, "\n")
  }
}
