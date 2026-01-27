#' Alpha Diversity (Per-Site)
#'
#' Compute Hill numbers for each site individually.
#'
#' @param x A site-by-species matrix (abundance data).
#' @param q Numeric vector. Orders of diversity. Default `c(0, 1, 2)`.
#'
#' @return A data.frame with columns for each q value, one row per site.
#'
#' @details
#' Alpha diversity represents local (within-site) diversity. For Hill numbers:
#' - q=0: Species richness
#' - q=1: Exponential of Shannon entropy
#' - q=2: Inverse Simpson concentration
#'
#' @references
#' Jost, L. (2007). Partitioning diversity into independent alpha and beta
#' components. Ecology, 88, 2427-2439.
#'
#' @seealso [gammaDiversity()] for regional diversity, [diversityPartition()]
#'   for full alpha-beta-gamma decomposition
#'
#' @examples
#' species <- matrix(rpois(50 * 30, 2), nrow = 50)
#' alpha <- alphaDiversity(species, q = c(0, 1, 2))
#' head(alpha)
#'
#' # Mean alpha diversity
#' colMeans(alpha)
#'
#' @export
alphaDiversity <- function(x, q = c(0, 1, 2)) {
  x <- as.matrix(x)
  n_sites <- nrow(x)
  result <- matrix(NA, nrow = n_sites, ncol = length(q))
  colnames(result) <- paste0("q", q)
 for (i in seq_len(n_sites)) {
    abundances <- as.numeric(x[i, ])
    for (j in seq_along(q)) {
      result[i, j] <- calc_hill_number(abundances, q[j])
    }
  }
  as.data.frame(result)
}
#' Gamma Diversity (Regional)
#'
#' Compute Hill numbers for the pooled community across all sites.
#'
#' @param x A site-by-species matrix (abundance data).
#' @param q Numeric vector. Orders of diversity. Default `c(0, 1, 2)`.
#'
#' @return A named numeric vector with gamma diversity for each q.
#'
#' @details
#' Gamma diversity represents regional (total) diversity across all sites.
#' It is computed by pooling abundances across all sites and calculating
#' Hill numbers on the combined community.
#'
#' @references
#' Jost, L. (2007). Partitioning diversity into independent alpha and beta
#' components. Ecology, 88, 2427-2439.
#'
#' @seealso [alphaDiversity()] for local diversity, [diversityPartition()]
#'   for full alpha-beta-gamma decomposition
#'
#' @examples
#' species <- matrix(rpois(50 * 30, 2), nrow = 50)
#' gammaDiversity(species, q = c(0, 1, 2))
#'
#' @export
gammaDiversity <- function(x, q = c(0, 1, 2)) {
  x <- as.matrix(x)
  # Pool abundances across all sites
  pooled <- colSums(x)
  result <- sapply(q, function(qi) calc_hill_number(pooled, qi))
  names(result) <- paste0("q", q)
  result
}
#' Alpha-Beta-Gamma Diversity Partitioning
#'
#' Decompose regional (gamma) diversity into local (alpha) and turnover (beta)
#' components using multiplicative partitioning of Hill numbers.
#'
#' @param x A site-by-species matrix (abundance data).
#' @param q Numeric vector. Orders of diversity. Default `c(0, 1, 2)`.
#' @param weights Character or numeric. Site weights for alpha calculation:
#'   - `"equal"`: Equal weights (default)
#'   - `"proportional"`: Weights proportional to site abundance
#'   - Numeric vector of custom weights
#'
#' @return An object of class `spacc_partition` containing:
#'   \item{alpha}{Mean alpha diversity (effective number of species per site)}
#'   \item{beta}{Beta diversity (effective number of communities)}
#'   \item{gamma}{Gamma diversity (regional species pool)}
#'   \item{q}{Orders of diversity}
#'   \item{n_sites}{Number of sites}
#'   \item{n_species}{Total species count}
#'
#' @details
#' This function implements Jost (2007) multiplicative partitioning:
#'
#' \deqn{\gamma = \alpha \times \beta}
#'
#' Where:
#' - **Alpha**: Mean effective number of species per site
#' - **Beta**: Effective number of distinct communities (1 = all identical,
#'   n_sites = all completely different)
#' - **Gamma**: Total effective number of species in the region
#'
#' Beta diversity is interpreted as the effective number of communities:
#' - Beta = 1: All sites have identical composition
#' - Beta = n_sites: Sites share no species
#'
#' @references
#' Jost, L. (2007). Partitioning diversity into independent alpha and beta
#' components. Ecology, 88, 2427-2439.
#'
#' Chao, A., Chiu, C.H. & Jost, L. (2014). Unifying species diversity, phylogenetic
#' diversity, functional diversity, and related similarity and differentiation
#' measures through Hill numbers. Annual Review of Ecology, Evolution, and
#' Systematics, 45, 297-324.
#'
#' @seealso [alphaDiversity()], [gammaDiversity()], [spaccBeta()] for
#'   spatial beta diversity accumulation
#'
#' @examples
#' # Simulated community data
#' species <- matrix(rpois(50 * 30, 2), nrow = 50)
#'
#' # Partition diversity
#' part <- diversityPartition(species, q = c(0, 1, 2))
#' print(part)
#'
#' # Beta near 1 = homogeneous region
#' # Beta near n_sites = heterogeneous region
#'
#' @export
diversityPartition <- function(x, q = c(0, 1, 2), weights = "equal") {
  x <- as.matrix(x)
  n_sites <- nrow(x)
  n_species <- ncol(x)
  # Compute per-site alpha
  alpha_per_site <- alphaDiversity(x, q)
  # Handle weights
  if (is.character(weights) && weights == "equal") {
    w <- rep(1 / n_sites, n_sites)
  } else if (is.character(weights) && weights == "proportional") {
    site_totals <- rowSums(x)
    w <- site_totals / sum(site_totals)
  } else if (is.numeric(weights)) {
    stopifnot(length(weights) == n_sites)
    w <- weights / sum(weights)
  } else {
    stop("weights must be 'equal', 'proportional', or a numeric vector")
  }
  # Weighted mean alpha (using generalized mean for Hill numbers)
  alpha <- sapply(seq_along(q), function(j) {
    if (q[j] == 1) {
      # Geometric mean for q=1
      exp(sum(w * log(alpha_per_site[[j]] + 1e-10)))
    } else {
      # Power mean for other q
      (sum(w * alpha_per_site[[j]]^(1 - q[j])))^(1 / (1 - q[j]))
    }
  })
  names(alpha) <- paste0("q", q)
  # Gamma diversity
  gamma <- gammaDiversity(x, q)
  # Beta = gamma / alpha (multiplicative partitioning)
  beta <- gamma / alpha
  names(beta) <- paste0("q", q)
  structure(
    list(
      alpha = alpha,
      beta = beta,
      gamma = gamma,
      alpha_per_site = alpha_per_site,
      q = q,
      weights = w,
      n_sites = n_sites,
      n_species = n_species,
      call = match.call()
    ),
    class = "spacc_partition"
  )
}
#' @export
print.spacc_partition <- function(x, ...) {
  cat("Alpha-Beta-Gamma Diversity Partitioning\n")
  cat(sprintf("%d sites, %d species\n\n", x$n_sites, x$n_species))
  df <- data.frame(
    q = x$q,
    alpha = round(x$alpha, 2),
    beta = round(x$beta, 2),
    gamma = round(x$gamma, 2)
  )
  print(df, row.names = FALSE)
  cat("\nInterpretation:\n")
  cat("  Alpha = mean effective species per site\n")
  cat("  Beta  = effective number of communities (1 to n_sites)\n")
  cat("  Gamma = regional effective species (gamma = alpha x beta)\n")
  invisible(x)
}
#' @export
summary.spacc_partition <- function(object, ...) {
  data.frame(
    q = object$q,
    alpha = object$alpha,
    beta = object$beta,
    gamma = object$gamma,
    n_sites = object$n_sites,
    beta_normalized = (object$beta - 1) / (object$n_sites - 1)
  )
}
