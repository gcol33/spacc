# Spatial Functional Diversity Accumulation

Compute spatial accumulation of functional diversity metrics based on
traits.

## Usage

``` r
spaccFunc(
  x,
  coords,
  traits,
  metric = c("fdis", "fric"),
  n_seeds = 50L,
  method = c("knn", "nn_walk"),
  distance = c("euclidean", "haversine"),
  parallel = TRUE,
  n_cores = NULL,
  progress = TRUE,
  seed = NULL,
  map = FALSE,
  focal_points = NULL,
  focal_domain = NULL,
  axis_length_k = NULL,
  itv = FALSE
)
```

## Arguments

- x:

  A site-by-species matrix (abundance data recommended).

- coords:

  A data.frame with columns `x` and `y`, or a `spacc_dist` object.

- traits:

  A species-by-traits matrix. Row names should match species (columns of
  x).

- metric:

  Character vector. Metrics to compute:

  - `"fdis"`: Functional Dispersion (mean distance to centroid)

  - `"fric"`: Functional Richness (convex hull volume approximation)

  - `"rao"`: Rao's quadratic entropy (abundance-weighted mean pairwise
    Euclidean trait distance)

  - `"K"`: aggregated functional diversity (Wojcik et al. 2024),
    reported with its four facets as `k_richness`, `k_biomass_evenness`,
    `k_trait_evenness`, `k_dispersion` and `K`

- n_seeds:

  Integer. Number of random starting points. Default 50.

- method:

  Character. Accumulation method: `"knn"` or `"nn_walk"`.

- distance:

  Character. Site distance method: `"euclidean"` or `"haversine"`.

- parallel:

  Logical. Use parallel processing? Default `TRUE`.

- n_cores:

  Integer. Number of cores.

- progress:

  Logical. Show progress? Default `TRUE`.

- seed:

  Integer. Random seed.

- map:

  Logical. If `TRUE`, run accumulation from every site as seed and store
  per-site final values for spatial mapping. Enables
  [`as_sf()`](https://gillescolling.com/spacc/reference/as_sf.md) and
  `plot(type = "map")`. Default `FALSE`.

- focal_points:

  Optional focal points for `method = "knn"`. See
  [`spacc()`](https://gillescolling.com/spacc/reference/spacc.md).

- focal_domain:

  Optional polygonal focal domain. See
  [`spacc()`](https://gillescolling.com/spacc/reference/spacc.md).

- axis_length_k:

  Numeric vector, one value per trait. Length of each trait axis for the
  richness facet of `"K"`. Default `NULL` uses the range of each trait
  over all species in `traits`.

- itv:

  Logical. Use intraspecific trait variability: `traits` is then a long
  data frame with one row per individual or population, a `species`
  column, an optional `site` column (matched like the sites of
  [`spaccNetwork()`](https://gillescolling.com/spacc/reference/spaccNetwork.md))
  and numeric trait columns. Default `FALSE`.

## Value

An object of class `spacc_func` containing:

- curves:

  Named list of matrices, one per metric (n_seeds x n_sites)

- metric:

  Metrics computed

- coords, n_seeds, n_sites, method:

  Parameters used

## Details

Functional diversity metrics quantify trait space occupation:

- **FDis (Functional Dispersion)**: Abundance-weighted mean distance
  from the community centroid in trait space. Captures functional
  divergence.

- **FRic (Functional Richness)**: Volume of trait space occupied (convex
  hull). Requires more species than traits to compute.

- **K (aggregated functional diversity)**: the geometric mean of four
  facets between 0 and 1, each equal to 1 when species are spread
  uniformly along the trait axis with equal abundances (Wojcik et al.
  2024). Per trait, with the \\n\\ present species sorted by trait value
  and relative abundances \\b_i\\: functional richness is the trait
  range over the axis length; biomass evenness compares the \\b_i\\ with
  \\1/n\\; trait evenness compares the gaps between neighbouring trait
  values with an even spacing; dispersion is the abundance-weighted
  trait variance over that of a uniform distribution on the same range,
  capped at 1. Several traits combine as the geometric mean over traits,
  and each reported facet is its geometric mean over traits. Trait
  evenness needs three species: with two, K is the cube root of the
  other three facets and `k_trait_evenness` is `NA`. With fewer than two
  species K is 0. With presence/absence data biomass evenness is 1.

**Intraspecific trait variability** (`itv = TRUE`; Palacio et al. 2025).
Each trait row is a unit of the community in place of a species: at each
step the metrics are evaluated on the rows describing the accumulated
sites, each weighted by its species' abundance at its site, shared
equally among the rows describing that species there. Rows without a
site (no `site` column, or `NA`) stand for the species at every site
where it has no site-specific rows, so one row per species reproduces
the species-mean curves. A widespread species whose traits change along
a gradient then adds trait range as the accumulation reaches new parts
of the gradient, which its mean trait value hides.

## References

Laliberté, E. & Legendre, P. (2010). A distance-based framework for
measuring functional diversity from multiple traits. Ecology, 91,
299-305.

Palacio, F.X., Ottaviani, G., Mammola, S., Graco-Roza, C., de Bello, F.
& Carmona, C.P. (2025). Integrating intraspecific trait variability in
functional diversity: an overview of methods and a guide for ecologists.
Ecological Monographs, 95, e70024.
[doi:10.1002/ecm.70024](https://doi.org/10.1002/ecm.70024)

Wojcik, L.A., Gaedke, U., van Velzen, E. & Klauschies, T. (2024).
Measuring overall functional diversity by aggregating its multiple
facets: functional richness, biomass evenness, trait evenness and
dispersion. Methods in Ecology and Evolution, 16, 215-227.
[doi:10.1111/2041-210X.14470](https://doi.org/10.1111/2041-210X.14470)

## See also

`FD::dbFD()` for comprehensive functional diversity analysis

## Examples

``` r
# \donttest{
coords <- data.frame(x = runif(50), y = runif(50))
species <- matrix(rpois(50 * 20, 2), nrow = 50)

# Trait matrix (species x traits)
traits <- matrix(rnorm(20 * 3), nrow = 20)
rownames(traits) <- paste0("sp", 1:20)
colnames(species) <- rownames(traits)

func <- spaccFunc(species, coords, traits, metric = c("fdis", "fric"))
plot(func)
# }
```
