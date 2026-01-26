# Plot Spatial SAC

Create a ggplot2 visualization of species accumulation curves.

## Usage

``` r
# S3 method for class 'spacc'
plot(
  x,
  ci = TRUE,
  ci_level = 0.95,
  ci_alpha = 0.3,
  show_seeds = FALSE,
  saturation = FALSE,
  saturation_level = 0.9,
  ...
)
```

## Arguments

- x:

  A `spacc`, `summary.spacc`, or `spacc_multi` object.

- ci:

  Logical. Show confidence interval ribbon? Default `TRUE`.

- ci_level:

  Numeric. Confidence level for interval. Default 0.95.

- ci_alpha:

  Numeric. Transparency of CI ribbon. Default 0.3.

- show_seeds:

  Logical. Show individual seed curves? Default `FALSE`.

- saturation:

  Logical. Mark saturation point? Default `FALSE`.

- saturation_level:

  Numeric. Proportion for saturation. Default 0.9.

- ...:

  Additional arguments (ignored).

## Value

A ggplot2 object.

## Examples

``` r
if (FALSE) { # \dontrun{
sac <- spacc(species, coords)
plot(sac)

# Customize
plot(sac, ci_alpha = 0.5, saturation = TRUE) +
  ggplot2::theme_dark()
} # }
```
