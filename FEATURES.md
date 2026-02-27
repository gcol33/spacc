# spacc Feature Plan

## Implemented (v0.2.0--v0.8.0)

All features from the original plan are now implemented:

| # | Feature | Status | Version |
|---|---------|--------|---------|
| 1 | Faith's PD accumulation | Done | v0.5.0 |
| 2 | Coverage-based extrapolation | Done | v0.6.0 |
| 3 | Diversity-Area Relationship (DAR) | Done | v0.5.0 |
| 4 | EVT SAR model | Done | v0.7.0 |
| 5 | Functional/phylogenetic beta | Done | v0.6.0 |
| 6 | SESARS | Done | v0.6.0 |
| 7 | SFAR | Done | v0.6.0 |
| 8 | Endemism-area | Done | v0.6.0 |

---

## Planned Features (v0.9.0+)

Based on a literature review of methods published 2023--2026. Each feature improves an existing module or adds a new analytical capability.

---

### 1. Sample-Based Coverage Estimator

**Status**: Not implemented

**Builds on**: `spaccCoverage()`

**What**: The current Good-Turing coverage estimator assumes independent random sampling, which is inappropriate for spatially structured data (kNN neighborhoods, plot-based surveys). Chiu (2023) derived a corrected estimator for sample-based (non-independently sampled) abundance data that accounts for spatial autocorrelation.

**Implementation**:
- Add `method = "sample-based"` option to `spaccCoverage()` alongside the current `"individual-based"` default
- Implement Chiu's modified Good-Turing formula that adjusts for non-independence between sampling units
- The correction factor depends on the number of sampling units and the frequency distribution of species across units (not individuals)
- Returns the same `spacc_coverage` object; the adjusted coverage values feed into downstream `interpolateCoverage()` and `extrapolateCoverage()`

**References**:
- Chiu, C.H. (2023). Sample coverage for sample-based abundance data. *Ecology*, 104, e4099.

**Complexity**: Low (formula change within existing function)

---

### 2. Coverage-Standardized Beta Diversity

**Status**: Not implemented

**Builds on**: `spaccBeta()`, `spaccBetaFunc()`, `spaccBetaPhylo()`

**What**: `spaccBeta()` computes beta diversity at fixed spatial scales (site counts). Coverage-based standardization removes sampling intensity artifacts, yielding pure compositional differentiation. The iNEXT.beta3D framework unifies taxonomic, phylogenetic, and functional beta diversity under Hill numbers, with coverage-based rarefaction and extrapolation.

**Implementation**:
- Add `standardize = "coverage"` argument to `spaccBeta()` (default `"none"` preserves current behavior)
- At each accumulation step, standardize alpha and gamma diversity to equal coverage before computing beta
- Compute four dissimilarity indices (Jaccard-type and Sorensen-type turnover, and their complements) as functions of coverage
- Extend to `spaccBetaFunc()` and `spaccBetaPhylo()` with the same `standardize` argument
- New plot panel showing beta diversity vs. coverage rather than vs. site count

**References**:
- Chao, A., Chiu, C.H., Villeger, S., et al. (2023). Rarefaction and extrapolation with beta diversity under a framework of Hill numbers: the iNEXT.beta3D standardization. *Ecological Monographs*, 93, e1588.

**CRAN packages**: [iNEXT.beta3D](https://cran.r-project.org/package=iNEXT.beta3D)

**Complexity**: High (requires coverage computation at each accumulation step for each beta component)

---

### 3. Aggregated Functional Diversity Index (K)

**Status**: Not implemented

**Builds on**: `spaccFunc()`

**What**: A single composite functional diversity metric that is the geometric mean of four independent FD facets: functional richness, biomass evenness, trait evenness, and dispersion. Each facet and K range from 0 to 1 (uniform distribution = maximum diversity). Overcomes limitations of existing indices that fail for simple communities or lack clear interpretation.

**Implementation**:
- Add `metric = "K"` option to `spaccFunc()` alongside existing `"FDis"` and `"FRic"`
- At each accumulation step, compute the four facets from the trait matrix of accumulated species
- Return the geometric mean as a single normalized trajectory
- The [0,1] normalization makes cross-site and cross-dataset comparison straightforward
- New `spacc_func` column or separate output for the composite index

**References**:
- Wojcik, R., et al. (2025). An aggregated overall functional diversity index. *Methods in Ecology and Evolution*, 16, 215--227.

**Complexity**: Medium (new metric computation, reuses existing trait handling)

---

### 4. EvoHeritage (Generalized Phylogenetic Diversity)

**Status**: Not implemented

**Builds on**: `spaccPhylo()`

**What**: Generalizes Faith's PD by accounting for both accumulation and attrition of evolutionary features along phylogenetic branches. Unlike PD (which assumes features only accumulate), EvoHeritage incorporates gradual loss of features through processes other than extinction (convergent evolution, trait lability). The metric smoothly interpolates between species richness (high attrition rate) and Faith's PD (zero attrition) via a single parameter.

**Implementation**:
- Add `metric = "evoheritage"` option to `spaccPhylo()` alongside `"PD"`, `"MPD"`, `"MNTD"`
- Accept an attrition rate parameter `phi` (0 = Faith's PD, Inf = species richness)
- At each accumulation step, compute the EvoHeritage score by weighting each branch by its expected feature retention
- Default `phi` could be estimated from trait data if available
- Returns EvoHeritage accumulation curve in the `spacc_phylo` object

**References**:
- Rosauer, D.F., et al. (2024). EvoHeritage: a generalization of phylogenetic diversity. *Systematic Biology*, 73, 158--174.

**Depends on**: ape (already in Suggests)

**Complexity**: Medium (new metric, requires branch-level weighting)

---

### 5. Multi-Resolution Beta Partitioning

**Status**: Not implemented

**Builds on**: `spaccBeta()`

**What**: Partition beta diversity (turnover vs. nestedness) at two spatial resolutions simultaneously: site-level pairwise and habitat-level pooled assemblages. Reveals whether observed turnover arises from fine-scale species replacement or coarse-scale biotic homogenization -- patterns invisible at a single spatial grain.

**Implementation**:
- Add `scale = c("site", "habitat")` argument to `spaccBeta()`, or a new `multiscaleBeta()` function
- Accept a `habitat` grouping vector that assigns sites to broader spatial units (e.g., landscape patches, regions)
- Compute beta decomposition at both resolutions and return a paired result
- Plot method shows turnover/nestedness at both scales side by side
- New class `spacc_beta_multiscale` or extend `spacc_beta` with a `scale` attribute

**References**:
- Jones, M.M., et al. (2025). Multi-resolution beta diversity partitioning. *Diversity and Distributions*, 31, e70080.

**Complexity**: Medium (two-level aggregation of existing beta computation)

---

### 6. Shift-and-Rotate Null Models

**Status**: Not implemented

**Builds on**: `sesars()`, `ses()`

**What**: Current null models for SES computation break species-environment associations by permuting species across sites, which destroys spatial autocorrelation. Shift-and-Rotate (S&R) generates spatially realistic null expectations by randomly translating and rotating the sampling area within environmental layers. Works with any sampling design (unlike torus translation which requires regular grids) and preserves the spatial covariance structure.

**Implementation**:
- Add `null_model = "shift_rotate"` option to `ses()` and `sesars()` alongside the current `"swap"` and `"frequency"` methods
- Accept an `env` argument (raster or matrix of environmental predictors)
- For each null iteration: randomly translate (dx, dy) and rotate (theta) the sampling coordinates within the environmental layer extent, then extract environmental values at the shifted positions
- Compare observed diversity metrics against the distribution of shifted-environment expectations
- Falls back to standard permutation if no environmental data provided

**References**:
- Ridder, J., Hardy, O.J. & Ovaskainen, O. (2024). Shift-and-Rotate: a novel null model approach for testing joint species-environment associations. *Methods in Ecology and Evolution*, 15, 1885--1898.

**Complexity**: Medium (coordinate transformation + raster extraction per null iteration)

---

### 7. Spatial Interaction Accumulation

**Status**: Not implemented

**Builds on**: New module

**What**: Extend spatial accumulation from species to species *interactions*. Track how network complexity (number of unique interactions, connectance, modularity) builds as the spatial neighborhood expands. Relevant for pollination networks, food webs, and host-parasite systems where interaction diversity matters alongside species diversity.

**Implementation**:
- New function `spaccNetwork(x, coords, interactions, ...)` where `interactions` is an edge list or interaction matrix
- At each accumulation step, compute network metrics for the subnetwork formed by species present in accumulated sites
- Metrics: interaction richness, connectance, links per species, modularity (optional)
- Coverage-based rarefaction and extrapolation for network data following the iNEXT.link framework
- New class `spacc_network` with print/summary/plot methods
- Plot: Interaction accumulation curve alongside species accumulation for comparison

**References**:
- Chiu, C.H., Chao, A., Vogel, S., Kriegel, P. & Thorn, S. (2023). Quantifying and estimating ecological network diversity based on incomplete sampling data. *Philosophical Transactions of the Royal Society B*, 378, 20220183.

**CRAN packages**: [iNEXT.link](https://github.com/AnneChao/iNEXT.link), [bipartite](https://cran.r-project.org/package=bipartite)

**Complexity**: High (new data type, network metrics at each step)

---

### 8. Functional Diversity-Area Relationships (FDARs)

**Status**: Not implemented

**Builds on**: `dar()`, `spaccFunc()`

**What**: Extend the DAR framework from Hill numbers to multiple functional diversity facets (FRic, FDiv, FEve, FDis) simultaneously. Shows how different FD facets scale differently with area -- FRic tracks species richness closely, but standardized FD metrics reveal independent scaling behaviors driven by trait filtering vs. competitive exclusion.

**Implementation**:
- New function `fdar(object, traits, q = c(0, 1, 2))` or extend `dar()` with `type = "functional"`
- At each accumulation step, compute FRic, FDiv, FEve from the trait matrix of accumulated species
- Area estimation via Voronoi tessellation (reuse from `dar()`)
- Plot: Multi-panel FD-area curves for each facet, optionally with log-log transformation
- Power-law fitting for each facet separately

**References**:
- Dias, R.A., et al. (2025). Spatial scaling of multiple dimensions of functional diversity. *Functional Ecology*, 39, 1065--1078.

**Complexity**: Medium (extends existing DAR machinery with new metrics)

---

### 9. Intraspecific Trait Variability in Functional Diversity

**Status**: Not implemented

**Builds on**: `spaccFunc()`, `spaccBetaFunc()`

**What**: Current functional diversity computation uses species-level mean traits. Incorporating intraspecific trait variability (ITV) makes functional accumulation curves more realistic, especially for widespread species whose traits vary across the spatial extent of a study. Uses probabilistic hypervolumes rather than point representations in trait space.

**Implementation**:
- Add `itv = TRUE` option to `spaccFunc()` and `spaccBetaFunc()`
- Accept a traits object with per-individual or per-population measurements (not just species means)
- Compute kernel density hypervolumes at each accumulation step using the accumulated individuals' trait values
- FRic becomes the hypervolume volume; FDis becomes the mean distance to centroid in probability space
- Falls back to species means when `itv = FALSE` (default, current behavior)

**References**:
- Palacio, F.X., et al. (2025). Integrating intraspecific trait variability in functional diversity. *Ecological Monographs*, 95, e70024.

**Depends on**: [hypervolume](https://cran.r-project.org/package=hypervolume) (Suggests)

**Complexity**: High (per-individual trait tracking, hypervolume computation at each step)

---

### 10. Phylogenetic/Functional Redundancy

**Status**: Not implemented

**Builds on**: `spaccPhylo()`, `spaccFunc()`

**What**: Track whether spatial accumulation adds genuinely novel evolutionary lineages or functionally distinct species, versus redundant ones. Redundancy = ratio of taxonomic Hill diversity to phylogenetic/functional Hill diversity. High redundancy means new species are closely related or functionally similar to those already accumulated.

**Implementation**:
- New function `redundancy(spacc_hill, spacc_phylo)` or `redundancy(spacc_hill, spacc_func)`
- At each accumulation step, compute R = D_tax / D_phylo (or D_func)
- R > 1 indicates phylogenetic/functional clustering (redundancy); R < 1 indicates overdispersion (complementarity)
- Plot: Redundancy trajectory overlaid on the accumulation curve
- Alternatively, integrate as a `metric = "redundancy"` option in `spaccPhylo()` and `spaccFunc()`

**References**:
- Alberdi, A. (2024). hilldiv2: Hill diversity with phylogenetic and functional redundancy. GitHub package.
- Cadotte, M.W. & Davies, T.J. (2016). *Phylogenies in Ecology*. Princeton University Press.

**Complexity**: Low (ratio of existing outputs)

---

## Implementation Priority

| # | Feature | Complexity | Builds on | Priority |
|---|---------|-----------|-----------|----------|
| 1 | Sample-based coverage estimator | Low | `spaccCoverage()` | High |
| 10 | Phylo/functional redundancy | Low | `spaccPhylo()`, `spaccFunc()` | High |
| 3 | Aggregated FD index (K) | Medium | `spaccFunc()` | High |
| 4 | EvoHeritage | Medium | `spaccPhylo()` | High |
| 5 | Multi-resolution beta | Medium | `spaccBeta()` | Medium |
| 6 | Shift-and-Rotate null models | Medium | `ses()`, `sesars()` | Medium |
| 8 | Functional DAR | Medium | `dar()`, `spaccFunc()` | Medium |
| 2 | Coverage-standardized beta | High | `spaccBeta()` | Medium |
| 9 | Intraspecific trait variability | High | `spaccFunc()` | Low |
| 7 | Spatial interaction accumulation | High | New module | Low |
