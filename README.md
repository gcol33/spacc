# spacc

> Spatial species accumulation curves. Fast.

High-performance spatial species accumulation curves using nearest-neighbor algorithms with a C++ backend.

## Installation

```r
# Install from GitHub
remotes::install_github("gcol33/spacc")
```

## Usage

```r
library(spacc)

# Your species occurrence data
# coords: data.frame with x, y columns
# species: data.frame with plot_id, species columns

sac <- fast_spatial_sac(
  species = species_data,
  coords = coord_data,
  n_seeds = 50,
  method = "knn"
)

# Summary statistics
sac_summary(sac)

# Plot
plot(sac)
```

## Why spacc?

Existing R packages for spatial species accumulation (e.g., mobr) use pure R loops. spacc uses:

- **C++ core** via Rcpp for 10-100x speedup
- **Parallel processing** via RcppParallel for multi-seed runs
- **Spatial indexing** for efficient nearest-neighbor search

## Methods

- `knn`: k-Nearest Neighbor - always visit the closest unvisited plot
- `kncn`: k-Nearest Centroid Neighbor - recalculate centroid after each step

## License

MIT
