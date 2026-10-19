# Autoplot Methods for spacc Objects

[`ggplot2::autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
methods for spacc objects. Each one calls the corresponding
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) method, which
already returns a ggplot2 object.

## Usage

``` r
# S3 method for class 'spacc'
autoplot(object, ...)
```

## Arguments

- object:

  A spacc object.

- ...:

  Additional arguments passed to the corresponding
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method.

## Value

A ggplot2 object.

## Examples

``` r
# \donttest{
if (requireNamespace("ggplot2", quietly = TRUE)) {
  coords <- data.frame(x = runif(30), y = runif(30))
  species <- matrix(rbinom(30 * 15, 1, 0.3), nrow = 30)
  sac <- spacc(species, coords, n_seeds = 5, progress = FALSE)
  ggplot2::autoplot(sac)
}
# }
```
