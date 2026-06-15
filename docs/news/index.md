# Changelog

## spacc 0.9.0

### New Features

- [`extrapolateArea()`](https://gillescolling.com/spacc/reference/extrapolateArea.md)
  extrapolates richness to a spatial extent larger than the one sampled,
  using the total-species (T-S) curve of Ugland, Gray & Ellingsen
  (2003): sites are partitioned into spatial subareas, the expected
  total richness of random subarea combinations is plotted against their
  convex-hull area, and an asymptotic model is fitted and extended to a
  target area with a median-bias-corrected site-bootstrap interval. New
  `spacc_area` class with
  [`print()`](https://rdrr.io/r/base/print.html)/[`summary()`](https://rdrr.io/r/base/summary.html)/[`plot()`](https://rdrr.io/r/graphics/plot.default.html)/[`predict()`](https://rdrr.io/r/stats/predict.html)/[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  methods.

- [`extrapolate()`](https://gillescolling.com/spacc/reference/extrapolate.md)
  gains a calibrated `interval = "bootstrap"` (now the default) that
  refits the chosen model across resampled seed curves and returns
  percentile bounds, replacing the over-confident `nls` profile interval
  as the default. `predict(fit, interval = "bootstrap")` returns a
  prediction band, and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws it.
  `interval = "profile"` recovers the old behaviour.

- [`extrapolate()`](https://gillescolling.com/spacc/reference/extrapolate.md)
  now reports goodness of fit (residual RMSE over the observed range),
  an extrapolation-range note, and the nonparametric
  [`chao2()`](https://gillescolling.com/spacc/reference/chao2.md)/[`iChao2()`](https://gillescolling.com/spacc/reference/iChao2.md)
  estimates alongside the asymptote. It warns when the asymptote exceeds
  the observed richness by more than `warn_ratio` (default 2), when it
  disagrees with
  [`chao2()`](https://gillescolling.com/spacc/reference/chao2.md) by
  more than 50%, and when
  [`predict()`](https://rdrr.io/r/stats/predict.html) is evaluated
  beyond ~2.5x the sampled effort.

### Bug Fixes

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a
  `spacc_endemism` object no longer lists “Endemic species” twice in the
  legend. The confidence-band fill was drawn as its own guide because
  the ribbon covered only the endemism series; that fill guide is now
  suppressed so a single line-colour legend remains.

- Removed a stray blank line beneath the y-axis title in
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) for
  `spacc_beta` objects.

- [`extrapolate()`](https://gillescolling.com/spacc/reference/extrapolate.md)
  asymptote confidence intervals are no longer over-confident. The
  previous [`confint()`](https://rdrr.io/r/stats/confint.html) used an
  `nls` profile interval on the (ultra-smooth) mean curve, ignoring
  across-seed variability; it could exclude the truth in every replicate
  of a recovery simulation.
  [`confint()`](https://rdrr.io/r/stats/confint.html) now returns the
  calibrated bootstrap interval by default (`method = "profile"` is
  still available). The documentation notes that the parametric
  asymptote can be biased high on clustered or under-sampled data and
  points to the nonparametric estimators for calibrated total-richness
  intervals ([\#1](https://github.com/gcol33/spacc/issues/1),
  [\#3](https://github.com/gcol33/spacc/issues/3),
  [\#5](https://github.com/gcol33/spacc/issues/5)).

### Testing

- Added parameter-recovery and CI-coverage tests on data simulated with
  a known truth: richness estimators recover the true pool size
  (jackknife intervals at ~nominal coverage; chao/iChao near-unbiased),
  [`extrapolate()`](https://gillescolling.com/spacc/reference/extrapolate.md)
  recovers truth on near-saturation data,
  [`extrapolateArea()`](https://gillescolling.com/spacc/reference/extrapolateArea.md)
  interpolates the observed area, and the Hill, phylogenetic,
  functional, and coverage accumulation endpoints match their defining
  formulae ([\#2](https://github.com/gcol33/spacc/issues/2)).

## spacc 0.8.3

### Documentation

- Expanded all seven vignettes from brief overviews into full worked
  guides (quickstart, diversity, rarefaction/standardization, community
  assembly, spatial analysis, extrapolation, richness estimation). Each
  now covers the underlying model, simulation with known ground truth,
  fitting, uncertainty, prediction/comparison, and practical guidance,
  and exercises the full exported API (including
  [`evenness()`](https://gillescolling.com/spacc/reference/evenness.md),
  [`diversityProfile()`](https://gillescolling.com/spacc/reference/diversityProfile.md),
  [`spatialEigenvectors()`](https://gillescolling.com/spacc/reference/spatialEigenvectors.md)/[`spatialPartition()`](https://gillescolling.com/spacc/reference/spatialPartition.md),
  [`wavefront()`](https://gillescolling.com/spacc/reference/wavefront.md),
  and
  [`compareModels()`](https://gillescolling.com/spacc/reference/compareModels.md)).

## spacc 0.8.2

### Improvements

- [`spaccDiversity()`](https://gillescolling.com/spacc/reference/spaccDiversity.md)
  objects gain a dedicated
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method with a
  metric-neutral default y-axis label (“Cumulative diversity”) and a
  `ylab`/`title` argument, so custom-metric curves are no longer
  labelled as species counts.

## spacc 0.8.1

### Improvements

- [`spaccPhylo()`](https://gillescolling.com/spacc/reference/spaccPhylo.md)
  and
  [`spaccFunc()`](https://gillescolling.com/spacc/reference/spaccFunc.md)
  now accept non-integer abundances (e.g. cover or biomass) for Rao’s Q
  and FDis weighting. The phylogenetic and functional accumulation
  backends use double-precision community matrices instead of truncating
  to integers. Presence-based metrics (MPD/MNTD/PD/FRic) are unaffected.

## spacc 0.8.0

### New Features

#### Rao’s Quadratic Entropy (v0.8.0)

- [`spaccPhylo()`](https://gillescolling.com/spacc/reference/spaccPhylo.md)
  and
  [`spaccFunc()`](https://gillescolling.com/spacc/reference/spaccFunc.md)
  gain a `"rao"` metric: abundance-weighted mean pairwise distance
  (phylogenetic cophenetic distance, or Euclidean trait distance),
  accumulated along the spatial curve.
- [`spaccPhylo()`](https://gillescolling.com/spacc/reference/spaccPhylo.md)
  no longer binarises its input, so abundance data now weights Rao;
  presence metrics (MPD/MNTD/PD) are unchanged.
- Exact-math recovery tests against the Rao definition
  `sum_i sum_j p_i p_j d_ij`.

#### Custom Diversity Metrics (v0.8.0)

- [`spaccDiversity()`](https://gillescolling.com/spacc/reference/spaccDiversity.md)
  accumulates any user-supplied index along a spatial ordering: at each
  step the cumulative community is passed to a function that returns a
  scalar. Supports `knn`, `kncn`, `random`, `radius`, and `collector`
  orderings, abundance or incidence input, and extra arguments.
- Returns a `spacc_diversity` object inheriting `spacc`, so the standard
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), and
  [`predict()`](https://rdrr.io/r/stats/predict.html) methods apply.

#### Arbitrary-Order Rarefaction (v0.8.0)

- [`rarefy()`](https://gillescolling.com/spacc/reference/rarefy.md) now
  accepts any `q >= 0` (q = 0, 1, 2 keep their exact estimators; other
  orders report the Hill number of order `q`) instead of silently
  falling back to richness.

## spacc 0.7.1

### New Features

#### User-Defined Accumulation Order (v0.7.1)

- [`spacc()`](https://gillescolling.com/spacc/reference/spacc.md) gains
  an `order` argument for supplying an explicit accumulation sequence,
  bypassing distance computation and seed sampling.
  - Accepts a single ordering vector, a list of vectors, or a matrix
    with one ordering per row (each row treated like a seed for
    uncertainty bounds).
  - Each ordering must be a permutation of `seq_len(nrow(x))`.
  - Backed by `cpp_order_parallel()`, which reuses the
    random-accumulation worker with caller-supplied orderings.

## spacc 0.7.0

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
- `cpp_knn_parallel_seeds()` - kNN with explicit seed indices
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
