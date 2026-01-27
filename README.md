# spacc

[![R-CMD-check](https://github.com/gcol33/spacc/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/gcol33/spacc/actions/workflows/R-CMD-check.yml)
[![codecov](https://codecov.io/gh/gcol33/spacc/branch/main/graph/badge.svg)](https://app.codecov.io/gh/gcol33/spacc)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Fast Spatial Species Accumulation Curves with C++ Performance**

The `spacc` package computes spatially-explicit species accumulation curves and diversity metrics. Unlike traditional SACs that assume random sampling order, `spacc` respects the spatial arrangement of sites using nearest-neighbor algorithms with a C++ backend (Rcpp/RcppParallel) for speed.

## Quick Start

```r
library(spacc)

coords <- data.frame(x = runif(100), y = runif(100))
species <- matrix(rpois(100 * 50, 2), nrow = 100)

# Spatial accumulation curve
sac <- spacc(species, coords, method = "knn", n_seeds = 50)
plot(sac)

# Diversity partitioning
part <- diversityPartition(species, q = c(0, 1, 2))
print(part)
```

## Statement of Need

Species accumulation curves are fundamental tools for understanding sampling completeness and comparing biodiversity. Standard SACs assume random sampling order, but real ecological sampling is spatially structured. `spacc` addresses this by computing SACs that respect spatial arrangement, with applications in:

- biodiversity assessment for spatially structured surveys,
- community ecology and scale-dependent diversity patterns,
- invasion modeling with wavefront and directional methods,
- conservation planning to identify high-diversity areas.

## Performance

`spacc` automatically selects the fastest backend for nearest-neighbor queries:

| Backend | Distance | Method | When |
|---------|----------|--------|------|
| **k-d tree** ([nanoflann](https://github.com/jlblancoc/nanoflann)) | Euclidean | O(log n) per query | >500 sites (auto) |
| **Ball tree** (custom) | Haversine | O(log n) per query | >500 sites + `distance = "haversine"` |
| **Exact** (brute-force) | Any | O(n) per query | ≤500 sites (auto) |

Override with `backend = "kdtree"` or `backend = "exact"`. The ball tree uses haversine distance natively — no coordinate projection needed for geographic data.

```r
# Geographic coordinates with haversine distance
sac <- spacc(species, coords_lonlat, distance = "haversine", n_seeds = 50)

# Force tree backend for smaller datasets
sac <- spacc(species, coords, backend = "kdtree", n_seeds = 50)
```

## Features

### Accumulation Curves

- **`spacc()`**: Spatial accumulation with 7 methods (`knn`, `kncn`, `random`, `radius`, `gaussian`, `cone`, `collector`)
- **`spaccHill()`**: Hill number accumulation (q = 0, 1, 2)
- **`spaccBeta()`**: Beta diversity with turnover/nestedness partitioning (Baselga 2010)
- **`spaccCoverage()`**: Coverage-based rarefaction (Chao & Jost 2012)
- **`spaccPhylo()`**: Phylogenetic diversity accumulation (MPD, MNTD)
- **`spaccFunc()`**: Functional diversity accumulation (FDis, FRic)
- **`wavefront()`**: Expanding radius accumulation for invasion modeling
- **`distanceDecay()`**: Species richness vs distance from focal points

### Diversity

- **`alphaDiversity()`**: Per-site Hill numbers (local diversity)
- **`gammaDiversity()`**: Pooled Hill numbers (regional diversity)
- **`diversityPartition()`**: Alpha-beta-gamma decomposition (Jost 2007)
- **`spaccMetrics()`**: Per-site accumulation metrics (slope, AUC, half-richness)

All diversity and accumulation functions support spatial mapping via `map = TRUE`, `plot(x, type = "map")`, and `as_sf()` for sf/GIS integration.

### Analysis

- **`extrapolate()`**: Fit asymptotic models (5 options)
- **`compare()`**: Statistical comparison between curves
- **`rarefy()`**: Rarefaction to common sampling effort
- **`subsample()`**: Spatial subsampling of sites
- **`coleman()`**, **`mao_tau()`**, **`spatialRarefaction()`**: Analytical (non-simulation) methods
- **`distances()`**: Pre-compute distance matrices

## Installation

```r
# install.packages("pak")
pak::pak("gcol33/spacc")
```

## Usage Examples

### Spatial Accumulation

```r
coords <- data.frame(x = runif(100), y = runif(100))
species <- matrix(rbinom(100 * 50, 1, 0.3), nrow = 100)

sac <- spacc(species, coords, method = "knn", n_seeds = 100)
plot(sac)

# Compare spatial vs random
sac_random <- spacc(species, coords, method = "random", n_seeds = 100)
comp <- compare(sac, sac_random)
plot(comp)
```

### Diversity Partitioning

```r
species <- matrix(rpois(100 * 50, 2), nrow = 100)

part <- diversityPartition(species, q = c(0, 1, 2))
part
#> Alpha-Beta-Gamma Diversity Partitioning
#> 100 sites, 50 species
#>
#>  q alpha  beta gamma
#>  0 12.45  4.01    50
#>  1  8.23  3.12 25.68
#>  2  6.54  2.87 18.77

# Per-site alpha with spatial mapping
alpha <- alphaDiversity(species, q = c(0, 1, 2), coords = coords)
plot(alpha, type = "map", q = 0)
```

### Hill Number & Beta Diversity Accumulation

```r
hill <- spaccHill(species, coords, q = c(0, 1, 2), n_seeds = 50, map = TRUE)
plot(hill)                         # Accumulation curves
plot(hill, type = "map", q = 0)    # Spatial map

beta <- spaccBeta(species, coords, index = "sorensen", n_seeds = 50, map = TRUE)
plot(beta, partition = TRUE)
plot(beta, type = "map")
```

### Coverage-Based Rarefaction

```r
cov <- spaccCoverage(species, coords, n_seeds = 50, map = TRUE)
plot(cov, xaxis = "coverage")
plot(cov, type = "map", metric = "final_coverage")
```

### Phylogenetic & Functional Diversity

```r
library(ape)
tree <- rtree(50)
colnames(species) <- tree$tip.label

phylo <- spaccPhylo(species, coords, tree,
                    metric = c("mpd", "mntd"), map = TRUE)
plot(phylo, type = "map", metric = "mpd")

traits <- matrix(rnorm(50 * 3), nrow = 50)
rownames(traits) <- colnames(species)
func <- spaccFunc(species, coords, traits,
                  metric = c("fdis", "fric"), map = TRUE)
plot(func, type = "map", metric = "fdis")
```

### Spatial Mapping with sf

```r
# Convert any diversity measure to sf for GIS integration
metrics <- spaccMetrics(species, coords,
                        metrics = c("slope_10", "half_richness", "auc"))
metrics_sf <- as_sf(metrics, crs = 32631)

# Works for all diversity types
hill_sf <- as_sf(hill, crs = 32631)
beta_sf <- as_sf(beta, crs = 32631)
```

## Documentation

- [Function Reference](https://gillescolling.com/spacc/reference/)

## Support

> "Software is like sex: it's better when it's free." — Linus Torvalds

I'm a PhD student who builds R packages in my free time because I believe good tools should be free and open. I started these projects for my own work and figured others might find them useful too.

If this package saved you some time, buying me a coffee is a nice way to say thanks. It helps with my coffee addiction.

[![Buy Me A Coffee](https://img.shields.io/badge/-Buy%20me%20a%20coffee-FFDD00?logo=buymeacoffee&logoColor=black)](https://buymeacoffee.com/gcol33)

## License

MIT (see the LICENSE.md file)

## Citation

```bibtex
@software{spacc,
  author = {Colling, Gilles},
  title = {spacc: Fast Spatial Species Accumulation Curves},
  year = {2025},
  url = {https://github.com/gcol33/spacc}
}
```
