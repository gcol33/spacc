# Changelog

## spacc 0.11.0

CRAN release: 2026-10-11

### New features

- `spaccPhylo(metric = "evoheritage", rho = )` accumulates EvoHeritage
  (Rosindell et al. 2024), phylogenetic diversity with attrition of
  evolutionary heritage at rate `rho`. `rho = 0` reproduces Faith’s PD;
  a vector `rho` gives one curve per rate
  ([\#7](https://github.com/gcol33/spacc/issues/7)).
- `chao2(x, n_units = )` estimates richness for sites drawn without
  replacement from a finite set of `n_units` units (MoRE, Chiu 2023),
  which reaches the observed richness when every unit is sampled instead
  of overestimating it. `extrapolate(n_units = )` compares the asymptote
  with MoRE and uses it as the reference for the disagreement warning
  ([\#14](https://github.com/gcol33/spacc/issues/14)).
- [`arc()`](https://gillescolling.com/spacc/reference/arc.md) estimates
  total richness with the accumulation rate curve estimator of
  Shestopaloff (2025): log-normal, log-logistic and Weibull fits to the
  rate at which new species appear, combined with capped weights, with a
  parametric bootstrap interval from a penalised Poisson mixture of the
  frequency counts. `data = "incidence"` runs it on site occupancy. The
  weights are the minimum of the paper’s loss reached by gradient flow
  from zero weights, the one the author’s reference code reaches
  ([\#13](https://github.com/gcol33/spacc/issues/13)).
- `spaccCoverage(coverage = "chiu", n_units = )` gives the matching
  sample coverage: the undetected abundance uses the MoRE estimate of
  the undetected species, so the coverage reaches 1 when every unit is
  sampled ([\#14](https://github.com/gcol33/spacc/issues/14)).
- `spaccHillBeta(standardize = "coverage", coverage = )` standardises
  gamma and alpha to a common sample coverage at every step
  (iNEXT.beta3D, Chao et al. 2023) and returns the four dissimilarities
  `1-CqN`, `1-UqN`, `1-VqN`, `1-SqN`. Final-step values match
  [`iNEXT.beta3D::iNEXTbeta3D()`](https://rdrr.io/pkg/iNEXT.beta3D/man/iNEXTbeta3D.html)
  ([\#6](https://github.com/gcol33/spacc/issues/6)).
- `spaccFunc(metric = "K")` accumulates the aggregated functional
  diversity K of Wojcik et al. (2024) with its four facets
  (`k_richness`, `k_biomass_evenness`, `k_trait_evenness`,
  `k_dispersion`); K is defined from two species on
  ([\#8](https://github.com/gcol33/spacc/issues/8)).
- [`spaccRedundancy()`](https://gillescolling.com/spacc/reference/spaccRedundancy.md)
  tracks phylogenetic or functional redundancy along the accumulation:
  the ratio of the taxonomic Hill number to the phylogenetic (Chao et
  al. 2010) or similarity-sensitive functional (Leinster & Cobbold
  2012. Hill number of the accumulated community
        ([\#9](https://github.com/gcol33/spacc/issues/9)).
- `spaccBeta(pool = )` repeats the turnover/nestedness partition between
  pooled assemblages (habitats, regions) along the same site orders,
  next to the site-level partition (Jones et al. 2025);
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) take
  `resolution = "site"`, `"pooled"` or (plot) `"both"`
  ([\#11](https://github.com/gcol33/spacc/issues/11)).
- `ses(null_model = "shift_rotate")` keeps the community and moves a
  fixed sampling design (user-supplied `focal_points` of a `knn`
  analysis) by random rotations and shifts (Ridder et al. 2024), testing
  the design’s position rather than the community’s spatial structure
  ([\#12](https://github.com/gcol33/spacc/issues/12)).
- [`as_spacc()`](https://gillescolling.com/spacc/reference/as_spacc.md)
  converts a numeric curve or a data frame column, including the output
  of
  [`Rarefy::directionalSAC()`](https://rdrr.io/pkg/Rarefy/man/directionalSAC.html)
  and
  [`Rarefy::rare_alpha()`](https://rdrr.io/pkg/Rarefy/man/rare_alpha.html).
  `spacc(method = "knn", focal_points = coords)` is checked to reproduce
  Rarefy’s directional curve exactly
  ([\#16](https://github.com/gcol33/spacc/issues/16)).
- [`spaccNetwork()`](https://gillescolling.com/spacc/reference/spaccNetwork.md)
  accumulates an interaction network over space: interaction richness
  and Hill numbers, sample coverage, the species curves of both levels
  and connectance along the same orders (Chiu et al. 2023)
  ([\#15](https://github.com/gcol33/spacc/issues/15)).
- `spaccFunc(itv = TRUE)` takes traits per individual or population (a
  long table with `species`, optional `site` and trait columns) and
  evaluates FDis, FRic, Rao’s Q and K on the trait rows of the
  accumulated sites; one row per species reproduces the species-mean
  curves ([\#17](https://github.com/gcol33/spacc/issues/17)).
- [`dar()`](https://gillescolling.com/spacc/reference/dar.md) is an S3
  generic: besides a site-by-species matrix it takes a `spacc`,
  `spacc_hill`, `spacc_phylo`, `spacc_func` or `spacc_redundancy` object
  and returns the diversity-area relationship of each of its curves,
  along the same site orders, with a power-law fit per curve
  ([`coef()`](https://rdrr.io/r/stats/coef.html)). The area kernels are
  C++ (incremental convex hull, Voronoi cells), so
  [`dar()`](https://gillescolling.com/spacc/reference/dar.md) no longer
  needs `sf`. The result carries `curves`, `metric` and the curve
  `object` it was built from (formerly `hill`)
  ([\#10](https://github.com/gcol33/spacc/issues/10)).
- New C++ engine for size- and coverage-based rarefaction and
  extrapolation of Hill numbers (`src/core/inext_core.h`), matching
  iNEXT.3D.

### Bug fixes

- [`ses()`](https://gillescolling.com/spacc/reference/ses.md) on
  `spacc_phylo` and `spacc_func` objects always stopped with “Too few
  valid permutations”: their null analyses were never rerun. The objects
  now carry their tree, traits and metric settings, and
  [`ses()`](https://gillescolling.com/spacc/reference/ses.md) reruns
  them.
- [`ses()`](https://gillescolling.com/spacc/reference/ses.md) reran
  `knn` analyses with newly drawn focal points; the nulls now reuse the
  observed design, and `spacc_beta` nulls keep the distance method.
- [`extrapolate()`](https://gillescolling.com/spacc/reference/extrapolate.md)
  and
  [`compareModels()`](https://gillescolling.com/spacc/reference/compareModels.md)
  passed unknown arguments (such as `progress = FALSE`) on to
  [`nls()`](https://rdrr.io/r/stats/nls.html), which treated them as
  data, so every model “failed to converge”. Unknown arguments are now
  an error. The extrapolation vignette’s model comparison had been all
  failures for this reason.
- `dar(area_method = "voronoi")` computed bounding-box areas; it now
  sums the Voronoi cells of the accumulated sites. The convex-hull area
  of one or two sites is 0 instead of their bounding box.
- `diversityProfile(tree = )` scaled branch lengths by the maximum
  root-to-tip depth. It now uses the abundance-weighted mean depth of
  Chao et al. (2010), which is the same for an ultrametric tree and
  keeps the profile equal to the taxonomic one on a star tree otherwise.
  Phylogenetic and functional profiles run on the C++ kernels shared
  with
  [`spaccRedundancy()`](https://gillescolling.com/spacc/reference/spaccRedundancy.md).
- [`rarefy()`](https://gillescolling.com/spacc/reference/rarefy.md) with
  `q = 1` or `q = 2` returned the Hill number of the full sample at
  every size (q = 1 plugged the sample proportions into Shannon, q = 2
  ignored the subsample size). It now returns the expected Hill number
  of a random subsample at every order, matching iNEXT.3D, and
  extrapolates sizes above the sample size.
- [`extrapolateCoverage()`](https://gillescolling.com/spacc/reference/extrapolateCoverage.md)
  with `q = 1` or `q = 2` interpolated the richness curve within the
  observed coverage and applied ad hoc corrections beyond it. It now
  reads the Hill curve of the requested order (a
  [`spaccHillCoverage()`](https://gillescolling.com/spacc/reference/spaccHillCoverage.md)
  object; a `spacc_coverage` object carries `q = 0` only) and
  extrapolates the pooled sample with the Chao et al. (2014) estimators,
  which also replace the reconstructed singleton counts used for
  `q = 0`.
- [`spaccPhylo()`](https://gillescolling.com/spacc/reference/spaccPhylo.md)
  PD matched species to tree tips by column position, so it was wrong
  when `colnames(x)` were not in `tree$tip.label` order and read past
  the community when the tree held more tips than `x` had species.
  Species are now matched to tips by name.

### Internals

- Every accumulation method is now an order generator
  (`src/orders.cpp`), and one kernel accumulates species along those
  orders (`src/accumulate.cpp`). The per-method serial/parallel/tree
  copies of the traversal (`knn.cpp`, `kncn.cpp`, `random.cpp`,
  `methods.cpp`, `metrics.cpp`) are gone;
  [`spacc()`](https://gillescolling.com/spacc/reference/spacc.md),
  [`spaccMetrics()`](https://gillescolling.com/spacc/reference/spaccMetrics.md),
  [`spaccDiversity()`](https://gillescolling.com/spacc/reference/spaccDiversity.md)
  and the metric front doors all take their site orders from
  `.accumulation_orders()`.
- The Hill, coverage, beta, phylogenetic and functional accumulations
  run on one templated worker (`OrderMetricWorker` in
  `src/parallel_utils.h`); each metric is a small struct with an
  `at_step()` method. The pure kernels live in `src/core/`
  (`phylo_core.h` and `func_core.h` are new) and are the single
  implementation behind both the exported single-community functions
  (`calc_mpd()`, `calc_fdis()`, …) and the accumulation workers.
- Faith’s PD precomputes the parent of every tree node once per call
  instead of scanning the edge list at every step.
- The S3 methods of the curve classes share `.curve_stats()`,
  `.curve_frame()`, `.plot_curves_by()`, `.plot_site_map()` and
  `.as_sf_site_values()`; the
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  methods are one delegating function registered per class.
- `spacc` objects carry the site `orders` their curves were accumulated
  along.

### Behaviour changes

- `method = "gaussian"` draws from R’s random number generator, so
  `seed` (or [`set.seed()`](https://rdrr.io/r/base/Random.html)) now
  reproduces its curves. It previously used an unseeded C++ generator.
- `method = "kncn"` with `distance = "haversine"` measures centroid
  distances with the haversine formula on both backends; the exact
  backend used planar distances on longitude/latitude.
- `.accumulation_orders()` computes a distance matrix only for the
  methods that need one, so
  [`spaccHill()`](https://gillescolling.com/spacc/reference/spaccHill.md),
  [`spaccPhylo()`](https://gillescolling.com/spacc/reference/spaccPhylo.md),
  [`spaccFunc()`](https://gillescolling.com/spacc/reference/spaccFunc.md),
  [`spaccCoverage()`](https://gillescolling.com/spacc/reference/spaccCoverage.md)
  and
  [`spaccEndemism()`](https://gillescolling.com/spacc/reference/spaccEndemism.md)
  no longer build an `n x n` matrix for the fixed-focus `knn` method.
- [`as_sf()`](https://gillescolling.com/spacc/reference/as_sf.md) and
  `plot(type = "map")` on an object computed without `map = TRUE` fail
  with one message for every class: “No map data. Rerun () with map =
  TRUE.”

## spacc 0.10.2

### API Changes

- `method = "knn"` now implements fixed-focus spatially constrained
  rarefaction: continuous focal points are sampled from the spatial
  domain and all sites are ranked by distance from each focus. Supply
  `focal_domain` for a known polygonal study boundary or `focal_points`
  for exact focal coordinates. The former greedy current-site traversal
  is now `method = "nn_walk"`. `method = "radius"` was removed because
  its observed-site fixed-focus ordering duplicated the corrected `knn`
  concept.
  [`spaccWavefront()`](https://gillescolling.com/spacc/reference/spaccWavefront.md)
  remains the radius-indexed accumulation front door.

- Diversity, coverage, beta, phylogenetic, functional, endemism, metric,
  and diversity-area calculations now consume shared site-order
  matrices. This keeps the selected spatial rule identical across
  downstream quantities.

- Spatial-accumulation curves use the
  [`spaccWavefront()`](https://gillescolling.com/spacc/reference/spaccWavefront.md)
  front door. The distance- and area-relationship functions keep their
  established names:
  [`distanceDecay()`](https://gillescolling.com/spacc/reference/distanceDecay.md),
  [`betaDecay()`](https://gillescolling.com/spacc/reference/betaDecay.md),
  [`zetaDiversity()`](https://gillescolling.com/spacc/reference/zetaDiversity.md),
  and [`dar()`](https://gillescolling.com/spacc/reference/dar.md).

- [`spaccBeta()`](https://gillescolling.com/spacc/reference/spaccBeta.md)
  gains `traits` and `tree` arguments and computes taxonomic,
  functional, or phylogenetic beta diversity from a single front door.
  Supplying `traits` gives the trait-weighted Baselga partition;
  supplying `tree` gives the branch-length-weighted (PhyloSor)
  partition.

- [`diversityProfile()`](https://gillescolling.com/spacc/reference/diversityProfile.md)
  gains `traits` and `tree` arguments and computes taxonomic, functional
  (Leinster-Cobbold), or phylogenetic (Chao et al.) Hill profiles from a
  single front door.

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

CRAN release: 2026-06-20

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
  `wavefront()`, and
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
  scalar. Supports `knn`, `kncn`, `nn_walk`, `random`, and `collector`
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
- `cpp_knn_parallel_seeds()` - legacy nearest-neighbour walk with
  explicit seed indices
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

- `wavefront()` - Expanding radius accumulation
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
