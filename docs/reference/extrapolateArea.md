# Extrapolate Richness to a Larger Area (Ugland T-S Curve)

Estimate species richness for a spatial extent larger than the one
sampled, using the total-species (T-S) curve of Ugland, Gray & Ellingsen
(2003). Sites are partitioned into spatial subareas; the T-S curve is
the expected total richness of random combinations of subareas plotted
against their (convex-hull) area; an asymptotic model is fitted to that
curve and extrapolated to a target area, with a site-bootstrap
confidence interval.

## Usage

``` r
extrapolateArea(
  x,
  coords,
  target_area = NULL,
  n_subareas = 10L,
  model = c("michaelis-menten", "asymptotic", "weibull", "lomolino", "logistic"),
  n_combos = 30L,
  R = 200L,
  level = 0.95,
  seed = NULL,
  progress = TRUE
)
```

## Arguments

- x:

  A site-by-species matrix (presence/absence or abundance; binarized).

- coords:

  A two-column matrix or data frame of site coordinates (x, y).

- target_area:

  Numeric. Area to extrapolate to, in the squared units of `coords`.
  Defaults to twice the observed (convex-hull) area.

- n_subareas:

  Integer. Number of spatial subareas to partition sites into (k-means
  on coordinates). Default 10; must be between 4 and the number of
  sites.

- model:

  Character. Area model to fit: `"michaelis-menten"` (default),
  `"asymptotic"`, `"weibull"`, `"lomolino"`, or `"logistic"`.

- n_combos:

  Integer. Number of random subarea combinations averaged per T-S point
  (all combinations are used when there are fewer). Default 30.

- R:

  Integer. Number of site-bootstrap replicates for the confidence
  interval. Default 200.

- level:

  Numeric. Confidence level. Default 0.95.

- seed:

  Optional integer for reproducibility.

- progress:

  Logical. Show a progress message. Default `TRUE`.

## Value

An object of class `spacc_area` with components:

- ts_curve:

  Data frame of the T-S curve (`n_sub`, `area`, `richness`,
  `richness_sd`)

- fit:

  The fitted `nls` area model

- observed_area, observed_richness:

  Convex-hull area and richness of the full sample

- target_area:

  Area extrapolated to

- predicted_richness, predicted_ci:

  Extrapolated richness and its bootstrap interval

- asymptote, asymptote_ci:

  Model asymptote and its bootstrap interval

## Details

Unlike
[`extrapolate()`](https://gillescolling.com/spacc/reference/extrapolate.md),
which extends the curve in *sample count*, this answers "how many
species in a region of a given *area*", including areas beyond the
sampled extent.

## References

Ugland, K.I., Gray, J.S. & Ellingsen, K.E. (2003). The
species-accumulation curve and estimation of species richness. Journal
of Animal Ecology, 72, 888-897.
[doi:10.1046/j.1365-2656.2003.00748.x](https://doi.org/10.1046/j.1365-2656.2003.00748.x)

## See also

[`extrapolate()`](https://gillescolling.com/spacc/reference/extrapolate.md)
for sample-count extrapolation;
[`chao2()`](https://gillescolling.com/spacc/reference/chao2.md) for
nonparametric richness.

## Examples

``` r
# \donttest{
coords <- data.frame(x = runif(150), y = runif(150))
species <- matrix(rbinom(150 * 60, 1, 0.15), nrow = 150)
area_fit <- extrapolateArea(species, coords, n_subareas = 8, progress = FALSE)
print(area_fit)
# }
```
