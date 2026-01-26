# Extrapolate Total Species Richness

Fit an asymptotic model to estimate total species richness beyond the
observed sampling effort.

## Usage

``` r
extrapolate(
  object,
  model = c("michaelis-menten", "lomolino", "asymptotic", "weibull", "logistic"),
  ...
)
```

## Arguments

- object:

  A `spacc` object.

- model:

  Character. Model to fit: `"michaelis-menten"` (default), `"lomolino"`,
  `"asymptotic"`, `"weibull"`, or `"logistic"`.

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
