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
