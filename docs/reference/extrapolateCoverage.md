# Extrapolate Diversity Beyond Observed Coverage

Hill numbers at target sample coverages along a coverage-tracking
accumulation, interpolated within the observed coverage and extrapolated
beyond it (Chao et al. 2014).

## Usage

``` r
extrapolateCoverage(x, target_coverage = c(0.9, 0.95, 0.99), q = 0)
```

## Arguments

- x:

  A `spacc_coverage` object from
  [`spaccCoverage()`](https://gillescolling.com/spacc/reference/spaccCoverage.md)
  (richness, `q = 0` only) or a `spacc_hill_coverage` object from
  [`spaccHillCoverage()`](https://gillescolling.com/spacc/reference/spaccHillCoverage.md)
  (any of its orders).

- target_coverage:

  Numeric vector of target coverage levels (0 to 1). Can exceed observed
  coverage for extrapolation. Default `c(0.90, 0.95, 0.99)`.

- q:

  Numeric. Diversity order: 0 (richness, default), or for a
  `spacc_hill_coverage` object one of its orders.

## Value

An object of class `spacc_coverage_ext` containing:

- richness:

  Matrix of interpolated/extrapolated Hill numbers of order `q` (n_seeds
  x n_targets)

- target_coverage:

  Target coverage levels

- q:

  Diversity order used

- observed_coverage:

  Mean observed final coverage

- observed_richness:

  Mean observed final Hill number of order `q`

## Details

A target at or below a seed's final coverage is read off that seed's
spatial accumulation curve by linear interpolation in coverage. Above
it, every seed has accumulated the same pooled sample, so the value is
the coverage-based extrapolation of the pooled abundances (Chao et al.
2014): the sample size reaching the target coverage is found from the
Chao-Jost coverage deficit, and the Hill number at that size approaches
the Chao1 (q = 0), Chao et al. (2013) Shannon (q = 1) or unbiased
Simpson (q = 2) asymptote; `target_coverage = 1` gives the asymptote
itself. Extrapolation is defined for integer `q` and for the
individual-based coverage (`coverage = "chao"`).

## References

Chao, A. & Jost, L. (2012). Coverage-based rarefaction and
extrapolation: standardizing samples by completeness rather than size.
Ecology, 93, 2533-2547.

Chao, A., Gotelli, N.J., Hsieh, T.C., et al. (2014). Rarefaction and
extrapolation with Hill numbers: a framework for sampling and estimation
in species diversity studies. Ecological Monographs, 84, 45-67.

## See also

[`spaccCoverage()`](https://gillescolling.com/spacc/reference/spaccCoverage.md),
[`interpolateCoverage()`](https://gillescolling.com/spacc/reference/interpolateCoverage.md)

## Examples

``` r
# \donttest{
coords <- data.frame(x = runif(50), y = runif(50))
species <- matrix(rpois(50 * 30, 2), nrow = 50)
cov <- spaccCoverage(species, coords)
ext <- extrapolateCoverage(cov, target_coverage = c(0.95, 0.99))
print(ext)
# }
```
