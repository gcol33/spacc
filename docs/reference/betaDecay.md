# Beta Distance-Decay

Compute the distance-decay of community similarity by fitting models to
pairwise beta dissimilarity as a function of geographic distance.

## Usage

``` r
betaDecay(
  x,
  coords,
  index = c("sorensen", "jaccard"),
  model = c("all", "exponential", "power", "linear"),
  distance = c("euclidean", "haversine"),
  n_bins = 50L,
  progress = TRUE
)
```

## Arguments

- x:

  A site-by-species matrix (presence/absence or abundance), or a `spacc`
  object (will use its stored coordinates).

- coords:

  A data.frame with columns `x` and `y`, or a `spacc_dist` object.

- index:

  Character. Dissimilarity index: `"sorensen"` (default) or `"jaccard"`.

- model:

  Character. Decay model to fit: `"exponential"`, `"power"`, `"linear"`,
  or `"all"` (default). When `"all"`, fits all three and selects the
  best by AIC.

- distance:

  Character. Distance method: `"euclidean"` (default) or `"haversine"`.

- n_bins:

  Integer. Number of distance bins for binned means in plots. Default
  50.

- progress:

  Logical. Show progress messages? Default `TRUE`.

## Value

An object of class `spacc_beta_decay` containing:

- pairs:

  Data.frame with columns `distance`, `dissimilarity`, `site_i`,
  `site_j`

- fits:

  Named list of fitted model objects

- best_model:

  Name of best model by AIC

- half_life:

  Distance at which similarity halves (exponential model only)

- coefficients:

  Data.frame of model coefficients and AIC

- index:

  Dissimilarity index used

- n_sites:

  Number of sites

- n_pairs:

  Number of pairwise comparisons

## Details

Three decay models are available:

- **Exponential**: \\\beta = 1 - a \cdot e^{-b \cdot d}\\

- **Power**: \\\beta = a \cdot d^b\\

- **Linear**: \\\beta = a + b \cdot d\\

The half-life (exponential model) is the distance at which similarity
decays to half its initial value: \\d\_{1/2} = \ln(2) / b\\.

## References

Nekola, J.C. & White, P.S. (1999). The distance decay of similarity in
biogeography and ecology. Journal of Biogeography, 26, 867-878.

Morlon, H., Chuyong, G., Condit, R., et al. (2008). A general framework
for the distance-decay of similarity in ecological communities. Ecology
Letters, 11, 904-917.

## See also

[`spaccBeta()`](https://gillescolling.com/spacc/reference/spaccBeta.md)
for spatial beta diversity accumulation,
[`distanceDecay()`](https://gillescolling.com/spacc/reference/distanceDecay.md)
for species similarity decay

## Examples

``` r
# \donttest{
coords <- data.frame(x = runif(30), y = runif(30))
species <- matrix(rbinom(30 * 20, 1, 0.3), nrow = 30)

bd <- betaDecay(species, coords)
print(bd)
plot(bd)
# }
```
