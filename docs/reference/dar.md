# Diversity-Area Relationship (DAR)

Diversity against cumulative area along a spatial accumulation, for Hill
numbers (Ma 2018) or for any curve spacc computes: species richness,
phylogenetic and functional metrics, redundancy.

## Usage

``` r
dar(x, ...)

# Default S3 method
dar(
  x,
  coords,
  q = c(0, 1, 2),
  n_seeds = 50L,
  method = "knn",
  area_method = c("convex_hull", "voronoi", "count"),
  distance = c("euclidean", "haversine"),
  parallel = TRUE,
  n_cores = NULL,
  progress = TRUE,
  seed = NULL,
  focal_points = NULL,
  focal_domain = NULL,
  ...
)

# S3 method for class 'spacc'
dar(x, area_method = c("convex_hull", "voronoi", "count"), ...)

# S3 method for class 'spacc_hill'
dar(x, area_method = c("convex_hull", "voronoi", "count"), ...)

# S3 method for class 'spacc_phylo'
dar(x, area_method = c("convex_hull", "voronoi", "count"), ...)

# S3 method for class 'spacc_func'
dar(x, area_method = c("convex_hull", "voronoi", "count"), ...)

# S3 method for class 'spacc_redundancy'
dar(x, area_method = c("convex_hull", "voronoi", "count"), ...)
```

## Arguments

- x:

  A site-by-species matrix (abundance data recommended), or a curve
  object: `spacc`, `spacc_hill`, `spacc_phylo`, `spacc_func` or
  `spacc_redundancy`. A matrix is first passed to
  [`spaccHill()`](https://gillescolling.com/spacc/reference/spaccHill.md)
  with `q`.

- ...:

  Arguments passed to methods.

- coords:

  A data.frame with columns `x` and `y`, or a `spacc_dist` object
  (matrix input only).

- q:

  Numeric vector. Diversity orders (matrix input only). Default
  `c(0, 1, 2)`.

- n_seeds:

  Integer. Number of random starting points. Default 50.

- method:

  Character. Accumulation method: `"knn"` or `"nn_walk"`.

- area_method:

  Character. How to measure the area covered by the first k sites of an
  order: `"convex_hull"` (area of their convex hull, 0 until three sites
  span a polygon), `"voronoi"` (summed area of their Voronoi cells
  within the bounding box of all sites) or `"count"` (number of sites).
  Areas are in squared coordinate units. Default `"convex_hull"`.

- distance:

  Character. Distance method. Default `"euclidean"`.

- parallel:

  Logical. Use parallel processing? Default `TRUE`.

- n_cores:

  Integer. Number of cores.

- progress:

  Logical. Show progress? Default `TRUE`.

- seed:

  Integer. Random seed.

- focal_points:

  Optional focal points for `method = "knn"`. See
  [`spacc()`](https://gillescolling.com/spacc/reference/spacc.md).

- focal_domain:

  Optional polygonal focal domain. See
  [`spacc()`](https://gillescolling.com/spacc/reference/spacc.md).

## Value

An object of class `spacc_dar` containing:

- curves:

  Named list of diversity curves (n_seeds x n_sites), one per metric or
  order

- metric:

  Names of the curves

- area:

  Matrix of cumulative areas (n_seeds x n_sites)

- power_law:

  Data frame with the power-law fit \\D = c A^z\\ of each mean curve
  (`metric`, `c`, `z`)

- q:

  Diversity orders for Hill-number input, else `NULL`

- area_method:

  Method used for area estimation

- source:

  Class of the curve object the DAR was built from

- object:

  That curve object

## Details

The DAR (Ma, 2018) generalizes the SAR by replacing species richness
(q=0) with Hill numbers of any order. This reveals how different facets
of diversity scale with area:

- q=0 (richness) recovers the classic SAR

- q=1 (Shannon) shows how common species diversity scales

- q=2 (Simpson) shows how dominant species diversity scales

The area is accumulated along the same site orders as the diversity
curve, so any curve object gives its own diversity-area relationship:
functional richness, dispersion and other facets can scale with area
differently from richness (Dias et al. 2025). The power-law exponent `z`
is the slope of log mean diversity against log mean area over the steps
where both are positive.

## References

Ma, Z.S. (2018). DAR (diversity-area relationship): extending classic
SAR for biodiversity and biogeography analyses. Ecology and Evolution,
8, 10023-10038.

Arrhenius, O. (1921). Species and area. Journal of Ecology, 9, 95-99.

Dias, R.A., Bastazini, V.A.G., Gianuca, A.T. & Matthews, T.J. (2025).
The spatial scaling of multiple dimensions of functional diversity in
habitat islands. Functional Ecology.
[doi:10.1111/1365-2435.70158](https://doi.org/10.1111/1365-2435.70158)

## See also

[`spaccHill()`](https://gillescolling.com/spacc/reference/spaccHill.md),
[`spaccFunc()`](https://gillescolling.com/spacc/reference/spaccFunc.md),
[`spaccPhylo()`](https://gillescolling.com/spacc/reference/spaccPhylo.md),
[`extrapolate()`](https://gillescolling.com/spacc/reference/extrapolate.md)

## Examples

``` r
# \donttest{
coords <- data.frame(x = runif(50), y = runif(50))
species <- matrix(rpois(50 * 30, 2), nrow = 50)

d <- dar(species, coords, q = c(0, 1, 2))
plot(d)

traits <- matrix(rnorm(30 * 2), nrow = 30)
func <- spaccFunc(species, coords, traits, metric = c("fdis", "fric"))
dar(func)
# }
```
