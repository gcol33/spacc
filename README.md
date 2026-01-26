<!-- badges: start -->
[![R-CMD-check](https://github.com/gcol33/spacc/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/gcol33/spacc/actions/workflows/R-CMD-check.yml)
[![codecov](https://codecov.io/gh/gcol33/spacc/branch/main/graph/badge.svg)](https://app.codecov.io/gh/gcol33/spacc)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

# spacc

**Fast spatial species accumulation curves with C++ performance.**

## Quick Start

```r
library(spacc)

# Site coordinates and species matrix
coords <- data.frame(x = runif(100), y = runif(100))
species <- matrix(rbinom(100 * 50, 1, 0.3), nrow = 100)

sac <- spacc(species, coords, method = "knn", n_seeds = 50)
plot(sac)
```

## Statement of Need

Species accumulation curves (SACs) are fundamental in ecology for understanding sampling completeness and comparing biodiversity across sites. Standard SACs assume random sampling, but real sampling is often spatially structured.

**spacc** computes spatial SACs that respect the spatial arrangement of sampling sites, providing more realistic estimates for:

- Assessing sampling completeness in ecological surveys
- Comparing biodiversity across regions with different spatial configurations
- Modeling invasion spread or range expansion
- Understanding distance-decay relationships in species composition

## Features

### Spatial Accumulation Methods

- `spacc()` - Main function with 7 spatial methods:
  - `knn`: k-Nearest Neighbor (always visit closest unvisited)
  - `kncn`: k-Nearest Centroid Neighbor (visit closest to growing centroid)
  - `random`: Random order (null model for comparison)
  - `radius`: Expand by distance from seed
  - `gaussian`: Probabilistic distance-weighted selection
  - `cone`: Directional expansion within angular constraint
  - `collector`: Sites in data order (no randomization)

### Additional Functions

- `wavefront()` - Expanding radius accumulation (invasion modeling)
- `distanceDecay()` - Species richness vs distance from focal points
- `distances()` - Pre-compute distance matrices (Euclidean or Haversine)

### Analytical Methods (No Simulation)

- `coleman()` - Coleman expected accumulation
- `mao_tau()` - Mao Tau (exact) expected curve
- `spatialRarefaction()` - Spatially-constrained rarefaction

### Analysis

- `extrapolate()` - Fit asymptotic models (Michaelis-Menten, Lomolino, etc.)
- `compare()` - Statistical comparison between curves
- `rarefy()` - Rarefaction to common sampling effort
- `subsample()` - Spatial subsampling of sites

## Installation

Development version:

```r
# install.packages("pak")
pak::pak("gcol33/spacc")
```

## Usage Examples

### Basic Spatial SAC

```r
library(spacc)

# Generate example data
set.seed(42)
coords <- data.frame(x = runif(100), y = runif(100))
species <- matrix(rbinom(100 * 50, 1, 0.3), nrow = 100)

# Compute spatial SAC with k-nearest neighbor method
sac <- spacc(species, coords, method = "knn", n_seeds = 100)
sac
#> spacc: 100 sites, 50 species, 100 seeds (knn)

plot(sac)
```

### Compare Spatial vs Random

```r
# Spatial accumulation
sac_spatial <- spacc(species, coords, method = "knn", n_seeds = 100)

# Random accumulation (null model)
sac_random <- spacc(species, coords, method = "random", n_seeds = 100)

# Statistical comparison
comp <- compare(sac_spatial, sac_random)
plot(comp)
```

### Extrapolate to Asymptote

```r
sac <- spacc(species, coords, method = "knn", n_seeds = 50)

# Fit asymptotic model
fit <- extrapolate(sac, model = "michaelis-menten")
summary(fit)
#> Asymptotic richness: 52.3 (95% CI: 48.1-56.5)

plot(fit)
```

### Geographic Coordinates

```r
# Use Haversine distance for lat/lon data
coords_geo <- data.frame(x = runif(100, -5, 5), y = runif(100, 45, 50))
sac <- spacc(species, coords_geo, distance = "haversine")
```

## Performance

spacc uses a C++ backend via Rcpp and RcppParallel for high performance:

| Sites | Species | Seeds | Time (spacc) | Time (R loop) |
|-------|---------|-------|--------------|---------------|
| 100   | 50      | 100   | 0.02s        | 2.1s          |
| 500   | 100     | 100   | 0.3s         | 45s           |
| 1000  | 200     | 100   | 1.2s         | 5+ min        |

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
  year = {2026},
  url = {https://CRAN.R-project.org/package=spacc},
  doi = {10.32614/CRAN.package.spacc}
}
```
