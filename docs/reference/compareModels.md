# Compare Multiple SAR Models

Fit all (or a subset of) asymptotic species-area models and compare them
using AIC, BIC, delta-AIC, and Akaike weights.

## Usage

``` r
compareModels(
  object,
  models = c("michaelis-menten", "lomolino", "asymptotic", "weibull", "logistic", "evt"),
  ...
)
```

## Arguments

- object:

  A `spacc` object.

- models:

  Character vector of models to fit. Defaults to all six:
  `"michaelis-menten"`, `"lomolino"`, `"asymptotic"`, `"weibull"`,
  `"logistic"`, `"evt"`.

- ...:

  Additional arguments passed to
  [`extrapolate()`](https://gillescolling.com/spacc/reference/extrapolate.md).

## Value

An object of class `spacc_model_compare` containing:

- table:

  Data frame with model comparison statistics

- fits:

  Named list of `spacc_fit` objects

- best_model:

  Name of the best model by AIC

- data:

  Mean-curve data frame used for fitting

- spacc:

  Original spacc object

## Examples

``` r
# \donttest{
coords <- data.frame(x = runif(50), y = runif(50))
species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)
sac <- spacc(species, coords)
cm <- compareModels(sac)
print(cm)
# }
```
