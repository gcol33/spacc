# Extrapolate Total Species Richness

Fit an asymptotic model to the spatial accumulation curve and estimate
total species richness beyond the observed sampling effort.

## Usage

``` r
extrapolate(
  object,
  model = c("michaelis-menten", "lomolino", "asymptotic", "weibull", "logistic", "evt"),
  interval = c("bootstrap", "profile", "none"),
  R = 200L,
  level = 0.95,
  compare = TRUE,
  warn_ratio = 2,
  ...
)
```

## Arguments

- object:

  A `spacc` object.

- model:

  Character. Model to fit: `"michaelis-menten"` (default), `"lomolino"`,
  `"asymptotic"`, `"weibull"`, `"logistic"`, or `"evt"` (Extreme Value
  Theory, Borda-de-Agua et al. 2025).

- interval:

  Character. How to build the asymptote confidence interval:
  `"bootstrap"` (default) refits the model across resampled seed curves
  and takes percentile bounds; `"profile"` uses the (over-confident)
  `nls` profile interval on the mean-curve fit; `"none"` skips it.

- R:

  Integer. Number of bootstrap refits when `interval = "bootstrap"`.
  Default 200.

- level:

  Numeric. Confidence level for the interval. Default 0.95.

- compare:

  Logical. If `TRUE` (default) and the object carries incidence
  frequencies, compare the asymptote to the nonparametric
  [`chao2()`](https://gillescolling.com/spacc/reference/chao2.md) /
  [`iChao2()`](https://gillescolling.com/spacc/reference/iChao2.md)
  estimates and flag large disagreement.

- warn_ratio:

  Numeric. Warn when the fitted asymptote exceeds the observed richness
  by more than this factor. Default 2. Set to `Inf` to silence.

- ...:

  Additional arguments passed to
  [`stats::nls()`](https://rdrr.io/r/stats/nls.html).

## Value

An object of class `spacc_fit` containing:

- asymptote:

  Estimated total species richness (asymptote of the model)

- asymptote_ci:

  Confidence interval for the asymptote

- model:

  Model name

- interval:

  Interval method used

- fit:

  The `nls` fit object

- aic:

  AIC of the model

- gof:

  List with residual `rmse` and relative `rmse_rel` over the observed
  range

- compare:

  Nonparametric `chao2` / `iChao2` estimates, if available

- boot:

  Bootstrap coefficient draws, if `interval = "bootstrap"`

## Bias caveat

The asymptote is a *parametric extrapolation* of the accumulation curve.
On clustered or strongly under-sampled presence-absence data it tends to
overestimate true richness, sometimes substantially, because the curve
is far from saturation. The bootstrap interval quantifies the
uncertainty of the fitted asymptote (curve-fit and across-seed
variability) and is much wider than the `nls` profile interval, but it
is *not* a calibrated interval for true total richness: it is centred on
a possibly biased point estimate. For calibrated total-richness
estimates prefer the nonparametric estimators
[`chao2()`](https://gillescolling.com/spacc/reference/chao2.md) /
[`iChao2()`](https://gillescolling.com/spacc/reference/iChao2.md), which
are unbiased on the same data. The printout shows their values alongside
the asymptote for comparison.

## References

Lomolino, M.V. (2000). Ecology's most general, yet protean pattern: the
species-area relationship. Journal of Biogeography, 27, 17-26.

Flather, C.H. (1996). Fitting species-accumulation functions and
assessing regional land use impacts on avian diversity. Journal of
Biogeography, 23, 155-168.

Borda-de-Agua, L., Whittaker, R.J., Cardoso, P., et al. (2025). Extreme
value theory explains the topography and scaling of the species-area
relationship. Nature Communications, 16, 5346.

## See also

[`extrapolateArea()`](https://gillescolling.com/spacc/reference/extrapolateArea.md)
for area-based extrapolation to a region larger than the one sampled;
[`chao2()`](https://gillescolling.com/spacc/reference/chao2.md),
[`iChao2()`](https://gillescolling.com/spacc/reference/iChao2.md) for
nonparametric richness.

## Examples

``` r
# \donttest{
coords <- data.frame(x = runif(50), y = runif(50))
species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)
sac <- spacc(species, coords)
fit <- extrapolate(sac)
print(fit)
# }
```
