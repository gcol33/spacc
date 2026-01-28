# Extrapolate Richness Beyond Observed Coverage

Predict species richness at coverage levels beyond the empirical
maximum, following the Chao et al. (2014) framework. Provides seamless
interpolation and extrapolation as a function of sample coverage.

## Usage

``` r
extrapolateCoverage(x, target_coverage = c(0.9, 0.95, 0.99), q = 0)
```

## Arguments

- x:

  A `spacc_coverage` object from
  [`spaccCoverage()`](https://gillescolling.com/spacc/reference/spaccCoverage.md).

- target_coverage:

  Numeric vector of target coverage levels (0 to 1). Can exceed observed
  coverage for extrapolation. Default `c(0.90, 0.95, 0.99)`.

- q:

  Numeric. Diversity order for extrapolation: 0 (richness, default), 1
  (Shannon), or 2 (Simpson).

## Value

An object of class `spacc_coverage_ext` containing:

- richness:

  Matrix of interpolated/extrapolated richness (n_seeds x n_targets)

- target_coverage:

  Target coverage levels

- q:

  Diversity order used

- observed_coverage:

  Mean observed final coverage

- observed_richness:

  Mean observed final richness

## Details

For targets within observed coverage, linear interpolation is used. For
targets beyond observed coverage, asymptotic estimators are applied:

- **q = 0**: Chao1 estimator: S_est = S_obs + f1^2 / (2 \* f2), where
  f1/f2 are singleton/doubleton counts. Extrapolation via coverage
  deficit.

- **q = 1**: Shannon extrapolation based on the Good-Turing frequency
  formula.

- **q = 2**: Simpson extrapolation using the unbiased estimator.

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
if (FALSE) { # \dontrun{
cov <- spaccCoverage(species, coords)
ext <- extrapolateCoverage(cov, target_coverage = c(0.95, 0.99, 1.0))
print(ext)
plot(ext)
} # }
```
