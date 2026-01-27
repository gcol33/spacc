# Spatial Species Accumulation Curves

Compute species accumulation curves using various spatial sampling
methods with C++ backend for performance.

## Usage

``` r
spacc(
  x,
  coords,
  n_seeds = 50L,
  method = c("knn", "kncn", "random", "radius", "gaussian", "cone", "collector"),
  distance = c("euclidean", "haversine"),
  backend = c("auto", "exact", "kdtree"),
  support = NULL,
  include_halo = TRUE,
  sigma = NULL,
  cone_width = pi/4,
  parallel = TRUE,
  n_cores = NULL,
  progress = TRUE,
  groups = NULL,
  seed = NULL
)
```

## Arguments

- x:

  A site-by-species matrix (rows = sites, cols = species) with
  presence/absence (0/1) or abundance data. Can also be a data.frame.

- coords:

  A data.frame with columns `x` and `y` containing site coordinates, or
  a `spacc_dist` object from
  [`distances()`](https://gillescolling.com/spacc/reference/distances.md).

- n_seeds:

  Integer. Number of random starting points for uncertainty
  quantification. Default 50.

- method:

  Character. Accumulation method:

  - `"knn"`: k-Nearest Neighbor (always visit closest unvisited)

  - `"kncn"`: k-Nearest Centroid Neighbor (visit closest to centroid)

  - `"random"`: Random order (null model)

  - `"radius"`: Expand by distance from seed

  - `"gaussian"`: Probabilistic selection weighted by distance

  - `"cone"`: Directional expansion within angular constraint

  - `"collector"`: Sites in data order (no randomization, single curve)

- distance:

  Character. Distance method: `"euclidean"` or `"haversine"`.

- backend:

  Character. Nearest-neighbor backend for `knn` and `kncn`:

  - `"auto"` (default): Uses exact (brute-force) for ≤500 sites, spatial
    tree for \>500 sites.

  - `"exact"`: Always use brute-force with precomputed distance matrix.

  - `"kdtree"`: Always use spatial tree. Uses k-d tree (nanoflann) for
    Euclidean distances and ball tree for haversine distances. Faster
    for large datasets, no distance matrix needed.

- support:

  Optional. Spatial support for core/halo classification via
  [`areaOfEffect::aoe()`](https://gcol33.github.io/areaOfEffect/reference/aoe.html).
  Can be:

  - Country name or ISO code: `"France"`, `"FR"`, `"FRA"`

  - Vector of countries: `c("France", "Germany")`

  - An `sf` polygon object

  - An `aoe_result` object (pre-computed) When provided, seeds are
    sampled only from "core" sites (inside support), while accumulation
    can expand into "halo" sites (buffer zone).

- include_halo:

  Logical. When `support` is provided, should halo sites be included in
  accumulation? Default `TRUE` (ecological boundary). Set to `FALSE` for
  political/hard boundary.

- sigma:

  Numeric. Bandwidth for Gaussian method. Default auto-calculated.

- cone_width:

  Numeric. Half-width in radians for cone method. Default pi/4.

- parallel:

  Logical. Use parallel processing? Default `TRUE`.

- n_cores:

  Integer. Number of cores. Default `NULL` uses `detectCores() - 1`.

- progress:

  Logical. Show progress bar? Default `TRUE`.

- groups:

  Optional. A factor, character, or integer vector of length `ncol(x)`
  assigning each species (column) to a group. When provided, separate
  accumulation curves are computed for each group using the **same
  spatial site ordering**, and a `spacc_multi` object is returned.
  Useful for comparing native vs alien species, families, or any
  categorical split. Default `NULL` (no grouping).

- seed:

  Integer. Random seed for reproducibility. Default `NULL`.

## Value

When `groups = NULL`, an object of class `spacc` containing:

- curves:

  Matrix of cumulative species counts (n_seeds x n_sites)

- coords:

  Original coordinates

- n_seeds:

  Number of seeds used

- method:

  Method used

- n_species:

  Total species in dataset

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage
sac <- spacc(species, coords)
plot(sac)

# Different methods
sac_knn <- spacc(species, coords, method = "knn")
sac_gauss <- spacc(species, coords, method = "gaussian", sigma = 10)
sac_cone <- spacc(species, coords, method = "cone", cone_width = pi/6)

# Compare to null model
sac_rand <- spacc(species, coords, method = "random")
comp <- compare(sac_knn, sac_rand)

# With spatial support (seeds from France, accumulate into neighbors)
sac_france <- spacc(species, coords, support = "France")

# Hard boundary (France only, no halo)
sac_france_only <- spacc(species, coords, support = "France", include_halo = FALSE)

# Grouped accumulation (e.g., native vs alien)
status <- ifelse(grepl("alien", colnames(species)), "alien", "native")
sac_grouped <- spacc(species, coords, groups = status, seed = 42)
plot(sac_grouped)  # Overlaid curves per group
} # }
```
