# spacc 0.8.3

## Documentation

- Expanded all seven vignettes from brief overviews into full worked
  guides (quickstart, diversity, rarefaction/standardization, community
  assembly, spatial analysis, extrapolation, richness estimation). Each
  now covers the underlying model, simulation with known ground truth,
  fitting, uncertainty, prediction/comparison, and practical guidance,
  and exercises the full exported API (including `evenness()`,
  `diversityProfile()`, `spatialEigenvectors()`/`spatialPartition()`,
  `wavefront()`, and `compareModels()`).


# spacc 0.8.2

## Improvements

- `spaccDiversity()` objects gain a dedicated `plot()` method with a
  metric-neutral default y-axis label ("Cumulative diversity") and a
  `ylab`/`title` argument, so custom-metric curves are no longer labelled
  as species counts.


# spacc 0.8.1

## Improvements

- `spaccPhylo()` and `spaccFunc()` now accept non-integer abundances (e.g.
  cover or biomass) for Rao's Q and FDis weighting. The phylogenetic and
  functional accumulation backends use double-precision community matrices
  instead of truncating to integers. Presence-based metrics (MPD/MNTD/PD/FRic)
  are unaffected.


# spacc 0.8.0

## New Features

### Rao's Quadratic Entropy (v0.8.0)
- `spaccPhylo()` and `spaccFunc()` gain a `"rao"` metric: abundance-weighted
  mean pairwise distance (phylogenetic cophenetic distance, or Euclidean trait
  distance), accumulated along the spatial curve.
- `spaccPhylo()` no longer binarises its input, so abundance data now weights
  Rao; presence metrics (MPD/MNTD/PD) are unchanged.
- Exact-math recovery tests against the Rao definition `sum_i sum_j p_i p_j d_ij`.

### Custom Diversity Metrics (v0.8.0)
- `spaccDiversity()` accumulates any user-supplied index along a spatial
  ordering: at each step the cumulative community is passed to a function that
  returns a scalar. Supports `knn`, `kncn`, `random`, `radius`, and
  `collector` orderings, abundance or incidence input, and extra arguments.
- Returns a `spacc_diversity` object inheriting `spacc`, so the standard
  `summary()`, `plot()`, `as.data.frame()`, and `predict()` methods apply.

### Arbitrary-Order Rarefaction (v0.8.0)
- `rarefy()` now accepts any `q >= 0` (q = 0, 1, 2 keep their exact estimators;
  other orders report the Hill number of order `q`) instead of silently
  falling back to richness.


# spacc 0.7.1

## New Features

### User-Defined Accumulation Order (v0.7.1)
- `spacc()` gains an `order` argument for supplying an explicit accumulation
  sequence, bypassing distance computation and seed sampling.
  - Accepts a single ordering vector, a list of vectors, or a matrix with one
    ordering per row (each row treated like a seed for uncertainty bounds).
  - Each ordering must be a permutation of `seq_len(nrow(x))`.
  - Backed by `cpp_order_parallel()`, which reuses the random-accumulation
    worker with caller-supplied orderings.


# spacc 0.7.0

## New Features

### Hill Numbers (v0.2.0)
- `spaccHill()` - Spatial accumulation with Hill numbers (q = 0, 1, 2)
  - q = 0: Species richness
  - q = 1: Exponential Shannon entropy (effective common species)
  - q = 2: Inverse Simpson (effective dominant species)
- Extends iNEXT framework to spatial accumulation

### Spatial Beta Diversity (v0.3.0)
- `spaccBeta()` - Beta diversity accumulation with partitioning
  - Total beta diversity (Sorensen or Jaccard)
  - Turnover component (species replacement)
  - Nestedness component (species loss)
- Based on Baselga (2010) framework

### Coverage-Based Rarefaction (v0.4.0)
- `spaccCoverage()` - Track sample coverage during accumulation
- `interpolateCoverage()` - Interpolate richness at target coverage levels
- Implements Chao & Jost (2012) Good-Turing coverage estimator

### Phylogenetic/Functional Diversity (v0.5.0)
- `spaccPhylo()` - Phylogenetic diversity accumulation
  - MPD: Mean Pairwise Distance
  - MNTD: Mean Nearest Taxon Distance
- `spaccFunc()` - Functional diversity accumulation
  - FDis: Functional Dispersion
  - FRic: Functional Richness (approximation)

### Per-Site Metrics & Heatmaps (v0.6.0)
- `spaccMetrics()` - Extract accumulation metrics per site
  - `slope_10`, `slope_25`: Initial accumulation slopes
  - `half_richness`, `richness_Npct`: Sites to reach richness thresholds
  - `auc`: Area under accumulation curve
- `as_sf()` - Convert metrics to sf for spatial analysis
- Heatmap plotting via `plot(type = "heatmap")`

### Spatial Support Integration (v0.7.0)
- New `support` parameter in `spacc()` for areaOfEffect integration
  - Accepts country names: `spacc(species, coords, support = "France")`
  - Accepts sf polygons or pre-computed `aoe_result` objects
- Seeds sampled from **core** sites only (inside support)
- Accumulation expands into **halo** sites (buffer zone) by default
- `include_halo = FALSE` for hard/political boundaries
- Eliminates edge effects at arbitrary administrative boundaries

## Internal Changes
- New C++ implementations: hill.cpp, beta.cpp, coverage.cpp, phylo.cpp, metrics.cpp
- `cpp_knn_parallel_seeds()` - kNN with explicit seed indices
- Added sf, areaOfEffect, ape, iNEXT, betapart to Suggests

---

# spacc 0.1.0

Initial release.

## Features

### Core Spatial Accumulation Methods
- `spacc()` - Main function with multiple spatial sampling methods:
  - `knn`: k-Nearest Neighbor (always visit closest unvisited)
  - `kncn`: k-Nearest Centroid Neighbor (visit closest to centroid)
  - `random`: Random order (null model)
  - `radius`: Expand by distance from seed
  - `gaussian`: Probabilistic selection weighted by distance
  - `cone`: Directional expansion within angular constraint
  - `collector`: Sites in data order

### Additional Accumulation Methods
- `wavefront()` - Expanding radius accumulation
- `distanceDecay()` - Distance-decay relationships

### Analytical Methods (No Simulation)
- `coleman()` - Coleman expected accumulation
- `mao_tau()` - Mao Tau (exact) expected accumulation
- `collector()` - Collector's curve (data order)
- `spatialRarefaction()` - Spatially-constrained rarefaction

### Analysis Functions
- `extrapolate()` - Fit asymptotic models (Michaelis-Menten, Lomolino, etc.)
- `compare()` - Statistical comparison between curves
- `rarefy()` - Rarefaction to common effort
- `subsample()` - Subsample sites spatially

### Infrastructure
- C++ backend via Rcpp for performance
- Parallel processing via RcppParallel
- Haversine distance support for geographic coordinates
