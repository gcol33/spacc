# spacc 0.5.0 (development)

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
