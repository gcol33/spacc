# Spatial Species Accumulation Curves

Compute species accumulation curves using various spatial sampling
methods with C++ backend for performance.

## Usage

``` r
spacc(
  x,
  coords,
  n_seeds = 50L,
  method = c("knn", "kncn", "nn_walk", "random", "gaussian", "cone", "collector"),
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
  time = NULL,
  w_space = 1,
  w_time = 1,
  seed = NULL,
  focal_points = NULL,
  focal_domain = NULL,
  order = NULL
)
```

## Arguments

- x:

  A site-by-species matrix (rows = sites, cols = species) with
  presence/absence (0/1) or abundance data. Can also be a data.frame.

- coords:

  Site coordinates. Can be:

  - A data.frame with columns `x` and `y`

  - An `sf` object with point geometries (CRS is preserved)

  - A `spacc_dist` object from
    [`distances()`](https://gillescolling.com/spacc/reference/distances.md)

- n_seeds:

  Integer. Number of random focal points or starting sites for
  uncertainty quantification. Default 50.

- method:

  Character. Accumulation method:

  - `"knn"`: fixed-focus spatially constrained rarefaction

  - `"kncn"`: k-Nearest Centroid Neighbor (visit closest to centroid)

  - `"nn_walk"`: nearest-neighbour walk (move from the current site to
    the closest unvisited site)

  - `"random"`: Random order (null model)

  - `"gaussian"`: Probabilistic selection weighted by distance

  - `"cone"`: Directional expansion within angular constraint

  - `"collector"`: Sites in data order (no randomization, single curve)

- distance:

  Character. Distance method: `"euclidean"` or `"haversine"`.

- backend:

  Character. Nearest-neighbor backend for `nn_walk` and `kncn`:

  - `"auto"` (default): Uses exact (brute-force) for \<=500 sites,
    spatial tree for \>500 sites.

  - `"exact"`: Always use brute-force with precomputed distance matrix.

  - `"kdtree"`: Always use spatial tree. Uses k-d tree (nanoflann) for
    Euclidean distances and ball tree for haversine distances. Faster
    for large datasets, no distance matrix needed.

- support:

  Optional. Spatial support for core/halo classification via
  [`areaOfEffect::aoe()`](https://gcol33.github.io/areaOfEffect/reference/aoe.html).
  Can be:

  - `"auto"`: Auto-detect countries and run per-country accumulation,
    returning a grouped `spacc` object with one curve per country

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
  spatial site ordering**, and a grouped `spacc` object is returned.
  Useful for comparing native vs alien species, families, or any
  categorical split. Default `NULL` (no grouping).

- time:

  Optional. Numeric vector of length `nrow(x)` giving a temporal
  coordinate for each site. When provided, a combined spatiotemporal
  distance matrix is computed as
  `w_space * d_spatial + w_time * d_temporal` and used for accumulation.
  Forces exact (brute-force) backend since spatial trees cannot handle
  composite distances. Only supported for methods that use a distance
  matrix: `"nn_walk"`, `"gaussian"`.

- w_space:

  Numeric. Weight for spatial distance when `time` is provided. Default
  1.

- w_time:

  Numeric. Weight for temporal distance when `time` is provided. Default
  1.

- seed:

  Integer. Random seed for reproducibility. Default `NULL`.

- focal_points:

  Optional data frame or `sf` point object with `x` and `y` coordinates.
  For `method = "knn"`, each row is a fixed focal point and sites are
  accumulated by increasing distance from it. When supplied, `n_seeds`
  is set to the number of focal points. Default `NULL` samples
  continuous focal points uniformly from the spatial domain.

- focal_domain:

  Optional `sf` or `sfc` polygon defining the domain from which
  continuous `knn` focal points are sampled. When omitted, the convex
  hull of the eligible site coordinates is used. For haversine
  distances, hull sampling uses a local equal-area projection. Supplying
  the study polygon preserves concavities, holes, and disconnected
  components in irregular sampling domains. Its coordinates must use the
  same coordinate reference system as `coords`.

- order:

  Optional user-defined accumulation order(s). When supplied, `method`
  is ignored and sites are accumulated in the given sequence(s) (no
  distance computation or seed sampling). Can be:

  - A numeric vector of length `nrow(x)`: a single ordering (one curve).

  - A list of such vectors, or a matrix with one ordering per row: each
    produces a curve, treated like seeds for uncertainty quantification.
    Each ordering must be a permutation of `seq_len(nrow(x))`. Cannot be
    combined with `support` or `time`. Default `NULL`.

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

- focal_points:

  Continuous focal points used by `method = "knn"`

## Details

The canonical `knn` method follows Chiarucci et al. (2009): each curve
uses one fixed focal point and accumulates sites by increasing distance
from that point. The `nn_walk` method is a greedy traversal whose
reference point moves to the most recently selected site. The two
methods represent different spatial sampling models.

## References

Arrhenius, O. (1921). Species and area. Journal of Ecology, 9, 95-99.

Scheiner, S.M. (2003). Six types of species-area curves. Global Ecology
and Biogeography, 12, 441-447.

Chiarucci, A., Bacaro, G., Rocchini, D., Ricotta, C., Palmer, M.W. &
Scheiner, S.M. (2009). Spatially constrained rarefaction: incorporating
the autocorrelated structure of biological communities into sample-based
rarefaction. Community Ecology, 10, 209-214.

## Examples

``` r
# \donttest{
coords <- data.frame(x = runif(50), y = runif(50))
species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)

# Basic usage
sac <- spacc(species, coords)
plot(sac)

# Different methods
sac_knn <- spacc(species, coords, method = "knn")
sac_rand <- spacc(species, coords, method = "random")
comp <- compare(sac_knn, sac_rand)
# }
```
