# spacc

> Small exact engines for scientific computing in R.

*species accumulating in the neighbourhood*

[![R-CMD-check](https://github.com/gcol33/spacc/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/gcol33/spacc/actions/workflows/R-CMD-check.yml)
[![codecov](https://codecov.io/gh/gcol33/spacc/branch/main/graph/badge.svg)](https://app.codecov.io/gh/gcol33/spacc)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Spatial species accumulation curves built by nearest-neighbour expansion in C++.**

Pick a starting site, then walk outward to its nearest neighbours, counting new species
as you go. `spacc` does that from many starting points (Rcpp/RcppParallel backend), so the
curve carries a confidence ribbon showing how richness depends on *where* you begin.
The classical curve in `vegan::specaccum()` shuffles sites in random order and never sees
the map. `spacc` reproduces that random-order curve and adds the spatial one, so you can
read the gap between them.

```r
library(spacc)

coords  <- data.frame(x = runif(100), y = runif(100))
species <- matrix(rbinom(100 * 50, 1, 0.3), nrow = 100)

# spatial accumulation: walk outward by nearest neighbour
sac <- spacc(species, coords, method = "knn", n_seeds = 100)
plot(sac)
```

## Spatial order, not random order

`vegan::specaccum()` accumulates sites in random or collector order, so geography drops out.
`spacc()` keeps it: each curve follows nearest-neighbour expansion from a real seed site, and
the spread across seeds is the spatial signal. Run both and compare them directly.

```r
sac        <- spacc(species, coords, method = "knn",    n_seeds = 100)
sac_random <- spacc(species, coords, method = "random", n_seeds = 100)

comp <- compare(sac, sac_random)   # permutation, bootstrap, or AUC test
plot(comp)
```

Seven expansion methods are available: `knn` and `kncn` (nearest-neighbour and
nearest-centroid walks), `radius`, `gaussian`, and `cone` (distance- and direction-weighted
growth), `random` (the classical null), and `collector` (data order). Geographic coordinates
use `distance = "haversine"` with no projection step.

## Backend selection

Nearest-neighbour queries pick a backend by site count, overridable with `backend`.

| Backend | Distance | Cost per query | When |
|---|---|---|---|
| Exact (brute force) | any | O(n) | <= 500 sites (auto) |
| k-d tree ([nanoflann](https://github.com/jlblancoc/nanoflann)) | Euclidean | O(log n) | > 500 sites (auto) |
| Ball tree | Haversine | O(log n) | > 500 sites, `distance = "haversine"` |

```r
sac <- spacc(species, coords_lonlat, distance = "haversine", n_seeds = 50)
sac <- spacc(species, coords, backend = "kdtree", n_seeds = 50)
```

## Diversity beyond richness

Accumulation works for more than species counts. Each function carries the spatial ordering
through to a different diversity measure, with `map = TRUE` and `as_sf()` for GIS output.

- **`spaccHill()`**: Hill numbers (q = 0, 1, 2) along the accumulation, extending the iNEXT
  framework to spatial order.
- **`spaccBeta()`**: beta diversity with turnover and nestedness components (Baselga 2010).
- **`spaccCoverage()`**: coverage-based rarefaction via the Good-Turing estimator (Chao & Jost 2012).
- **`spaccPhylo()` / `spaccFunc()`**: phylogenetic (MPD, MNTD) and functional (FDis, FRic) accumulation.
- **`alphaDiversity()` / `gammaDiversity()` / `diversityPartition()`**: local, regional, and
  alpha-beta-gamma decomposition (Jost 2007).

```r
hill <- spaccHill(species, coords, q = c(0, 1, 2), n_seeds = 50, map = TRUE)
plot(hill)                       # accumulation curves
plot(hill, type = "map", q = 0)  # per-site richness map

beta <- spaccBeta(species, coords, index = "sorensen", n_seeds = 50)
plot(beta, partition = TRUE)     # turnover vs nestedness
```

## Grouped and spatiotemporal accumulation

Split species into groups (native vs alien, families, any factor) and accumulate each under
the **same** spatial site ordering, or add a temporal axis to weight space against time.

```r
# same spatial walk, separate curve per group
sac_grouped <- spacc(species, coords, groups = status, seed = 42)
plot(sac_grouped, facet = TRUE)

# sites sampled across years: weighted space-time distance
sac_st <- spacc(species, coords, method = "knn",
                time = site_years, w_space = 1, w_time = 0.5)
```

## Fitting and standardizing

- **`extrapolate()`**: asymptotic richness models (Michaelis-Menten, Lomolino, Weibull,
  logistic, EVT).
- **`compare()`**: permutation, bootstrap, or AUC tests between curves.
- **`rarefy()` / `subsample()`**: rarefaction to common effort, spatial subsampling.
- **`coleman()` / `mao_tau()` / `spatialRarefaction()`**: analytical (non-simulation) curves.
- **`as_spacc()`**: bring an existing `vegan::specaccum()` object into the same plotting and
  comparison machinery.

## Installation

```r
install.packages("pak")
pak::pak("gcol33/spacc")
```

## Documentation

- [Getting Started](https://gillescolling.com/spacc/articles/quickstart.html)
- [Diversity Accumulation](https://gillescolling.com/spacc/articles/diversity.html)
- [Extrapolation](https://gillescolling.com/spacc/articles/extrapolation.html)
- [Spatial Analysis](https://gillescolling.com/spacc/articles/spatial-analysis.html)
- [Community Assembly](https://gillescolling.com/spacc/articles/community-assembly.html)
- [Richness Estimation](https://gillescolling.com/spacc/articles/richness-estimation.html)
- [Rarefaction](https://gillescolling.com/spacc/articles/rarefaction-standardization.html)
- [Function Reference](https://gillescolling.com/spacc/reference/)

## Support

> "Software is like sex: it's better when it's free." — Linus Torvalds

I'm a PhD student who builds R packages in my free time because I believe good tools
should be free and open. I started these projects for my own work and figured others
might find them useful too.

If this package saved you some time, buying me a coffee is a nice way to say thanks.
It helps with my coffee addiction.

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
