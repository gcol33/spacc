# Plot Spatial SAC

Create a ggplot2 visualization of species accumulation curves. For
grouped `spacc` objects, curves are overlaid with different colors.

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
  facet = FALSE,
  ...
)
```

## Arguments

- x:

  A `spacc` or `summary.spacc` object.

- ci:

  Logical. Show confidence interval ribbon? Default `TRUE`.

- ci_level:

  Numeric. Confidence level for interval. Default 0.95.

- ci_alpha:

  Numeric. Transparency of CI ribbon. Default 0.3.

- show_seeds:

  Logical. Show individual seed curves? Default `FALSE`. Only available
  for ungrouped objects.

- saturation:

  Logical. Mark saturation point? Default `FALSE`. Only available for
  ungrouped objects.

- saturation_level:

  Numeric. Proportion for saturation. Default 0.9.

- facet:

  Logical. Use faceted panels for grouped objects? Default `FALSE`
  (overlay with colors). Ignored for ungrouped objects.

- ...:

  Additional arguments (ignored).

## Value

A ggplot2 object.

## Examples

``` r
if (FALSE) { # \dontrun{
sac <- spacc(species, coords)
plot(sac)

# Grouped
sac_grouped <- spacc(species, coords, groups = status)
plot(sac_grouped)
plot(sac_grouped, facet = TRUE)
} # }
```
