# Extrapolate Total Species Richness

Fit an asymptotic model to estimate total species richness beyond the
observed sampling effort.

## Usage

``` r
extrapolate(
  object,
  model = c("michaelis-menten", "lomolino", "asymptotic", "weibull", "logistic", "evt"),
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

- ...:

  Additional arguments passed to
  [`stats::nls()`](https://rdrr.io/r/stats/nls.html).

## Value

An object of class `spacc_fit` containing:

- asymptote:

  Estimated total species richness

- asymptote_ci:

  Confidence interval for asymptote

- model:

  Model name

- fit:

  The nls fit object

- aic:

  AIC of the model

## References

Lomolino, M.V. (2000). Ecology's most general, yet protean pattern: the
species-area relationship. Journal of Biogeography, 27, 17-26.

Flather, C.H. (1996). Fitting species-accumulation functions and
assessing regional land use impacts on avian diversity. Journal of
Biogeography, 23, 155-168.

Borda-de-Agua, L., Whittaker, R.J., Cardoso, P., et al. (2025). Extreme
value theory explains the topography and scaling of the species-area
relationship. Nature Communications, 16, 5346.

## Examples

``` r
if (FALSE) { # \dontrun{
sac <- spacc(species, coords)
fit <- extrapolate(sac, model = "lomolino")

print(fit)
plot(fit)
predict(fit, n = c(100, 500, 1000))
} # }
```
