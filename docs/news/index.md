# Changelog

## spacc 0.5.0 (development)

### New Features

#### Hill Numbers (v0.2.0)

- [`spaccHill()`](https://gillescolling.com/spacc/reference/spaccHill.md) -
  Spatial accumulation with Hill numbers (q = 0, 1, 2)
  - q = 0: Species richness
  - q = 1: Exponential Shannon entropy (effective common species)
  - q = 2: Inverse Simpson (effective dominant species)
- Extends iNEXT framework to spatial accumulation

#### Spatial Beta Diversity (v0.3.0)

- [`spaccBeta()`](https://gillescolling.com/spacc/reference/spaccBeta.md) -
  Beta diversity accumulation with partitioning
  - Total beta diversity (Sorensen or Jaccard)
  - Turnover component (species replacement)
  - Nestedness component (species loss)
- Based on Baselga (2010) framework

#### Coverage-Based Rarefaction (v0.4.0)

- [`spaccCoverage()`](https://gillescolling.com/spacc/reference/spaccCoverage.md) -
  Track sample coverage during accumulation
- [`interpolateCoverage()`](https://gillescolling.com/spacc/reference/interpolateCoverage.md) -
  Interpolate richness at target coverage levels
- Implements Chao & Jost (2012) Good-Turing coverage estimator

#### Phylogenetic/Functional Diversity (v0.5.0)

- [`spaccPhylo()`](https://gillescolling.com/spacc/reference/spaccPhylo.md) -
  Phylogenetic diversity accumulation
  - MPD: Mean Pairwise Distance
  - MNTD: Mean Nearest Taxon Distance
- [`spaccFunc()`](https://gillescolling.com/spacc/reference/spaccFunc.md) -
  Functional diversity accumulation
  - FDis: Functional Dispersion
  - FRic: Functional Richness (approximation)

#### Per-Site Metrics & Heatmaps (v0.6.0)

- [`spaccMetrics()`](https://gillescolling.com/spacc/reference/spaccMetrics.md) -
  Extract accumulation metrics per site
  - `slope_10`, `slope_25`: Initial accumulation slopes
  - `half_richness`, `richness_Npct`: Sites to reach richness thresholds
  - `auc`: Area under accumulation curve
- [`as_sf()`](https://gillescolling.com/spacc/reference/as_sf.md) -
  Convert metrics to sf for spatial analysis
- Heatmap plotting via `plot(type = "heatmap")`

#### Spatial Support Integration (v0.7.0)

- New `support` parameter in
  [`spacc()`](https://gillescolling.com/spacc/reference/spacc.md) for
  areaOfEffect integration
  - Accepts country names: `spacc(species, coords, support = "France")`
  - Accepts sf polygons or pre-computed `aoe_result` objects
- Seeds sampled from **core** sites only (inside support)
- Accumulation expands into **halo** sites (buffer zone) by default
- `include_halo = FALSE` for hard/political boundaries
- Eliminates edge effects at arbitrary administrative boundaries

### Internal Changes

- New C++ implementations: hill.cpp, beta.cpp, coverage.cpp, phylo.cpp,
  metrics.cpp
- [`cpp_knn_parallel_seeds()`](https://gillescolling.com/spacc/reference/cpp_knn_parallel_seeds.md) -
  kNN with explicit seed indices
- Added sf, areaOfEffect, ape, iNEXT, betapart to Suggests

------------------------------------------------------------------------

## spacc 0.1.0

Initial release.

### Features

#### Core Spatial Accumulation Methods

- [`spacc()`](https://gillescolling.com/spacc/reference/spacc.md) - Main
  function with multiple spatial sampling methods:
  - `knn`: k-Nearest Neighbor (always visit closest unvisited)
  - `kncn`: k-Nearest Centroid Neighbor (visit closest to centroid)
  - `random`: Random order (null model)
  - `radius`: Expand by distance from seed
  - `gaussian`: Probabilistic selection weighted by distance
  - `cone`: Directional expansion within angular constraint
  - `collector`: Sites in data order

#### Additional Accumulation Methods

- [`wavefront()`](https://gillescolling.com/spacc/reference/wavefront.md) -
  Expanding radius accumulation
- [`distanceDecay()`](https://gillescolling.com/spacc/reference/distanceDecay.md) -
  Distance-decay relationships

#### Analytical Methods (No Simulation)

- [`coleman()`](https://gillescolling.com/spacc/reference/coleman.md) -
  Coleman expected accumulation
- [`mao_tau()`](https://gillescolling.com/spacc/reference/mao_tau.md) -
  Mao Tau (exact) expected accumulation
- [`collector()`](https://gillescolling.com/spacc/reference/collector.md) -
  Collector’s curve (data order)
- [`spatialRarefaction()`](https://gillescolling.com/spacc/reference/spatialRarefaction.md) -
  Spatially-constrained rarefaction

#### Analysis Functions

- [`extrapolate()`](https://gillescolling.com/spacc/reference/extrapolate.md) -
  Fit asymptotic models (Michaelis-Menten, Lomolino, etc.)
- [`compare()`](https://gillescolling.com/spacc/reference/compare.md) -
  Statistical comparison between curves
- [`rarefy()`](https://gillescolling.com/spacc/reference/rarefy.md) -
  Rarefaction to common effort
- [`subsample()`](https://gillescolling.com/spacc/reference/subsample.md) -
  Subsample sites spatially

#### Infrastructure

- C++ backend via Rcpp for performance
- Parallel processing via RcppParallel
- Haversine distance support for geographic coordinates
