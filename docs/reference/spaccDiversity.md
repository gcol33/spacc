# Spatial Accumulation of a Custom Diversity Metric

Accumulate any user-supplied diversity index along a spatial ordering of
sites. At each accumulation step the cumulative community is passed to
`fun`, which returns a single number. This is the general escape hatch
behind the built-in metric functions: use it for indices that spacc does
not implement directly.

## Usage

``` r
spaccDiversity(
  x,
  coords,
  fun,
  ...,
  method = c("knn", "kncn", "random", "radius", "collector"),
  incidence = FALSE,
  n_seeds = 50L,
  distance = c("euclidean", "haversine"),
  progress = TRUE,
  seed = NULL
)
```

## Arguments

- x:

  A site-by-species matrix (rows = sites, cols = species), abundance or
  presence/absence.

- coords:

  A data.frame with columns `x` and `y`, or a `spacc_dist` object.

- fun:

  A function applied to the cumulative community at each step. It
  receives a named numeric vector of length `ncol(x)` (cumulative summed
  abundances, or 0/1 incidences when `incidence = TRUE`) plus any
  arguments passed through `...`, and must return a single numeric
  value.

- ...:

  Additional arguments passed to `fun`.

- method:

  Character. Spatial ordering of sites: `"knn"` (default), `"kncn"`,
  `"random"`, `"radius"`, or `"collector"`.

- incidence:

  Logical. If `TRUE`, `fun` receives 0/1 incidences instead of summed
  abundances. Default `FALSE`.

- n_seeds:

  Integer. Number of random starting points / orderings. Ignored for
  `"collector"` (a single data-order curve). Default 50.

- distance:

  Character. `"euclidean"` or `"haversine"`.

- progress:

  Logical. Show progress? Default `TRUE`.

- seed:

  Integer. Random seed for reproducibility.

## Value

An object of class `spacc_diversity` that inherits from `spacc`, so the
standard [`summary()`](https://rdrr.io/r/base/summary.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html),
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) and
[`predict()`](https://rdrr.io/r/stats/predict.html) methods apply.
`curves` is an `n_seeds x n_sites` matrix of the metric along the
accumulation.

## Details

The site ordering reuses the same spatial traversals as the built-in
methods (nearest-neighbour, nearest-centroid, random, distance-rank, or
data order), then evaluates `fun` on the accumulating community. Because
the index is an arbitrary R function, this trades the speed of the
compiled metrics for full flexibility.

## See also

[`spaccHill()`](https://gillescolling.com/spacc/reference/spaccHill.md),
[`spaccPhylo()`](https://gillescolling.com/spacc/reference/spaccPhylo.md),
[`spaccFunc()`](https://gillescolling.com/spacc/reference/spaccFunc.md)
for built-in metrics.

## Examples

``` r
# \donttest{
coords <- data.frame(x = runif(40), y = runif(40))
species <- matrix(rpois(40 * 20, 2), nrow = 40)

# Shannon entropy along the accumulation
shannon <- function(comm) {
  p <- comm[comm > 0] / sum(comm)
  -sum(p * log(p))
}
div <- spaccDiversity(species, coords, shannon, n_seeds = 20)
plot(div)
# }
```
