# spacc Implementation Plan

## Overview

Fast spatial species accumulation curves with C++ backend.

## Architecture

    spacc/
    ├── R/
    │   ├── spacc-package.R      # Package docs, Rcpp/RcppParallel imports
    │   ├── fast_spatial_sac.R   # Main user-facing function
    │   ├── sac_summary.R        # Summary methods
    │   ├── plot.R               # ggplot2 visualization
    │   └── utils.R              # Input validation, helpers
    ├── src/
    │   ├── distance.cpp         # Fast distance matrix computation
    │   ├── knn.cpp              # k-Nearest Neighbor accumulation
    │   ├── kncn.cpp             # k-Nearest Centroid Neighbor accumulation
    │   ├── parallel.cpp         # RcppParallel workers for multi-seed
    │   └── RcppExports.cpp      # Auto-generated
    └── tests/
        └── testthat/
            ├── test-distance.R
            ├── test-knn.R
            └── test-benchmark.R

## Implementation Phases

### Phase 1: Core C++ Distance Computation

**Goal:** Fast pairwise distance matrix

**Files:** - `src/distance.cpp` - `R/fast_distance_matrix.R`

**Implementation:**

``` cpp
// Euclidean distance, O(n²) but vectorized
NumericMatrix fast_dist(NumericVector x, NumericVector y);
```

**Benchmark target:** 10x faster than
[`dist()`](https://rdrr.io/r/stats/dist.html) for n \> 1000

------------------------------------------------------------------------

### Phase 2: kNN Accumulation (Single Seed)

**Goal:** Single nearest-neighbor traversal

**Files:** - `src/knn.cpp` - `R/fast_spatial_sac.R` (partial)

**Algorithm:**

    1. Start at seed plot
    2. Mark visited, record species
    3. Find nearest unvisited (use pre-computed distances)
    4. Repeat until all visited
    5. Return cumulative species vector

**Key optimization:** Pre-compute distance matrix once, reuse for all
lookups

------------------------------------------------------------------------

### Phase 3: Multi-Seed Parallelization

**Goal:** Run n_seeds iterations in parallel

**Files:** - `src/parallel.cpp` - `R/fast_spatial_sac.R` (complete)

**Implementation:**

``` cpp
// RcppParallel worker
struct SacWorker : public Worker {
  // Input: distance matrix, species matrix, seed indices
  // Output: matrix of accumulation curves (n_seeds x n_plots)
};
```

**Benchmark target:** Near-linear scaling with cores

------------------------------------------------------------------------

### Phase 4: Summary Statistics

**Goal:** Mean curve, CI, saturation metrics

**Files:** - `R/sac_summary.R`

**Outputs:** - Mean accumulation curve - 2.5% and 97.5% quantiles (95%
CI) - Saturation point (where curve reaches X% of max) - Coefficient of
variation at midpoint

------------------------------------------------------------------------

### Phase 5: Visualization

**Goal:** ggplot2 method for sac objects

**Files:** - `R/plot.R`

**Features:** - Mean line with CI ribbon - Optional: individual seed
curves - Optional: saturation point marker

------------------------------------------------------------------------

### Phase 6: kNCN Method (Optional)

**Goal:** Centroid-recalculating variant

**Files:** - `src/kncn.cpp`

**Difference from kNN:** After each step, recalculate centroid of
visited plots, find nearest unvisited to centroid (not to current plot)

------------------------------------------------------------------------

## Data Structures

### Input

``` r

# species: data.frame or matrix
#   - Rows: plots
#   - Cols: species (presence/absence or abundance)
# OR
#   - Long format: plot_id, species_id columns

# coords: data.frame
#   - x: longitude or easting
#   - y: latitude or northing
```

### Output (S3 class “spatial_sac”)

``` r

list(
  curves = matrix,      # n_seeds x n_plots, cumulative species
  coords = data.frame,  # original coordinates
  n_seeds = integer,
  method = character,   # "knn" or "kncn"
  call = call
)
```

------------------------------------------------------------------------

## Benchmarks to Beat

| Dataset size | mobr (est.) | spacc target |
|--------------|-------------|--------------|
| 100 plots    | ~1s         | \<0.1s       |
| 1,000 plots  | ~30s        | \<1s         |
| 10,000 plots | ~30min      | \<30s        |

------------------------------------------------------------------------

## Dependencies

**Runtime:** - Rcpp (C++ interface) - RcppParallel (parallel workers) -
ggplot2 (plotting)

**Build:** - C++11 compiler - GNU make

------------------------------------------------------------------------

## Testing Strategy

1.  **Correctness:** Compare output to mobr on small datasets
2.  **Performance:** Benchmark against mobr, pure R implementation
3.  **Edge cases:** Single plot, single species, empty data

------------------------------------------------------------------------

## Milestones

Phase 1: Distance matrix (1 day)

Phase 2: kNN single seed (1 day)

Phase 3: Parallel multi-seed (1 day)

Phase 4: Summary stats (0.5 day)

Phase 5: Plotting (0.5 day)

Phase 6: kNCN (optional, 1 day)

Documentation & vignette (1 day)

CRAN submission prep (1 day)
