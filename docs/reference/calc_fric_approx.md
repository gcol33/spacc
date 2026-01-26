# Calculate Functional Richness (convex hull volume approximation)

Uses simplified volume estimation based on trait ranges. For true FRic,
use R's geometry::convhulln.

## Usage

``` r
calc_fric_approx(traits, species_present)
```

## Arguments

- traits:

  Numeric matrix (species x traits)

- species_present:

  Logical vector

## Value

Approximate functional volume
