# Package index

## Core Functions

Main spatial accumulation functions

- [`spacc()`](https://gillescolling.com/spacc/reference/spacc.md) :
  Spatial Species Accumulation Curves
- [`wavefront()`](https://gillescolling.com/spacc/reference/wavefront.md)
  : Wavefront Expansion Accumulation
- [`distanceDecay()`](https://gillescolling.com/spacc/reference/distanceDecay.md)
  : Distance-Decay Analysis
- [`distances()`](https://gillescolling.com/spacc/reference/distances.md)
  : Compute Distance Matrix

## Hill Numbers

Diversity accumulation with Hill numbers (q=0,1,2)

- [`spaccHill()`](https://gillescolling.com/spacc/reference/spaccHill.md)
  : Spatial Accumulation with Hill Numbers

## Beta Diversity

Spatial beta diversity with turnover/nestedness partitioning

- [`spaccBeta()`](https://gillescolling.com/spacc/reference/spaccBeta.md)
  : Spatial Beta Diversity Accumulation

## Coverage-Based Rarefaction

Standardize by sample completeness

- [`spaccCoverage()`](https://gillescolling.com/spacc/reference/spaccCoverage.md)
  : Coverage-Based Spatial Rarefaction
- [`interpolateCoverage()`](https://gillescolling.com/spacc/reference/interpolateCoverage.md)
  : Interpolate Richness at Target Coverage Levels

## Phylogenetic & Functional Diversity

Accumulation of PD and FD metrics

- [`spaccPhylo()`](https://gillescolling.com/spacc/reference/spaccPhylo.md)
  : Spatial Phylogenetic Diversity Accumulation
- [`spaccFunc()`](https://gillescolling.com/spacc/reference/spaccFunc.md)
  : Spatial Functional Diversity Accumulation

## Per-Site Metrics

Extract metrics per starting site for heatmap visualization

- [`spaccMetrics()`](https://gillescolling.com/spacc/reference/spaccMetrics.md)
  : Per-Site Accumulation Metrics
- [`as_sf()`](https://gillescolling.com/spacc/reference/as_sf.md) :
  Convert spacc_metrics to sf

## Analytical Methods

Expected curves without simulation

- [`coleman()`](https://gillescolling.com/spacc/reference/coleman.md) :
  Coleman Expected Accumulation
- [`mao_tau()`](https://gillescolling.com/spacc/reference/mao_tau.md) :
  Exact (Mao Tau) Expected Accumulation
- [`collector()`](https://gillescolling.com/spacc/reference/collector.md)
  : Collector's Curve
- [`spatialRarefaction()`](https://gillescolling.com/spacc/reference/spatialRarefaction.md)
  : Spatially-Constrained Rarefaction

## Analysis

Post-hoc analysis of accumulation curves

- [`extrapolate()`](https://gillescolling.com/spacc/reference/extrapolate.md)
  : Extrapolate Total Species Richness
- [`compare()`](https://gillescolling.com/spacc/reference/compare.md) :
  Compare Two Accumulation Curves
- [`rarefy()`](https://gillescolling.com/spacc/reference/rarefy.md) :
  Individual-Based Rarefaction
- [`subsample()`](https://gillescolling.com/spacc/reference/subsample.md)
  : Spatial Subsampling

## Internal

Internal functions (not for end users)

- [`cpp_beta_knn_parallel()`](https://gillescolling.com/spacc/reference/cpp_beta_knn_parallel.md)
  : Parallel kNN Beta Diversity Accumulation
- [`cpp_beta_knn_single()`](https://gillescolling.com/spacc/reference/cpp_beta_knn_single.md)
  : Single kNN Beta Diversity Accumulation
- [`cpp_collector_single()`](https://gillescolling.com/spacc/reference/cpp_collector_single.md)
  : Collector Accumulation (data order)
- [`cpp_cone_parallel()`](https://gillescolling.com/spacc/reference/cpp_cone_parallel.md)
  : Parallel Directional Cone Accumulation
- [`cpp_cone_single()`](https://gillescolling.com/spacc/reference/cpp_cone_single.md)
  : Single Directional Cone Accumulation
- [`cpp_distance_decay_parallel()`](https://gillescolling.com/spacc/reference/cpp_distance_decay_parallel.md)
  : Parallel Distance-Decay Accumulation
- [`cpp_distance_decay_single()`](https://gillescolling.com/spacc/reference/cpp_distance_decay_single.md)
  : Distance-Decay Accumulation
- [`cpp_distance_matrix()`](https://gillescolling.com/spacc/reference/cpp_distance_matrix.md)
  : Fast Distance Matrix
- [`cpp_func_knn_parallel()`](https://gillescolling.com/spacc/reference/cpp_func_knn_parallel.md)
  : Parallel Functional Diversity Accumulation
- [`cpp_func_knn_single()`](https://gillescolling.com/spacc/reference/cpp_func_knn_single.md)
  : Single kNN Functional Diversity Accumulation
- [`cpp_gaussian_parallel()`](https://gillescolling.com/spacc/reference/cpp_gaussian_parallel.md)
  : Parallel Gaussian-weighted Accumulation
- [`cpp_gaussian_single()`](https://gillescolling.com/spacc/reference/cpp_gaussian_single.md)
  : Single Gaussian-weighted Accumulation
- [`cpp_kncn_metrics_parallel()`](https://gillescolling.com/spacc/reference/cpp_kncn_metrics_parallel.md)
  : Parallel kNCN Metrics Accumulation
- [`cpp_kncn_parallel()`](https://gillescolling.com/spacc/reference/cpp_kncn_parallel.md)
  : Parallel kNCN Accumulation
- [`cpp_kncn_single()`](https://gillescolling.com/spacc/reference/cpp_kncn_single.md)
  : Single kNCN Accumulation Curve
- [`cpp_knn_coverage_parallel()`](https://gillescolling.com/spacc/reference/cpp_knn_coverage_parallel.md)
  : Parallel kNN Accumulation with Coverage
- [`cpp_knn_coverage_single()`](https://gillescolling.com/spacc/reference/cpp_knn_coverage_single.md)
  : Single kNN Accumulation with Coverage Tracking
- [`cpp_knn_hill_parallel()`](https://gillescolling.com/spacc/reference/cpp_knn_hill_parallel.md)
  : Parallel kNN Accumulation with Hill Numbers
- [`cpp_knn_hill_single()`](https://gillescolling.com/spacc/reference/cpp_knn_hill_single.md)
  : Single kNN Accumulation with Hill Numbers
- [`cpp_knn_metrics_parallel()`](https://gillescolling.com/spacc/reference/cpp_knn_metrics_parallel.md)
  : Parallel kNN Metrics Accumulation
- [`cpp_knn_parallel()`](https://gillescolling.com/spacc/reference/cpp_knn_parallel.md)
  : Parallel kNN Accumulation
- [`cpp_knn_parallel_seeds()`](https://gillescolling.com/spacc/reference/cpp_knn_parallel_seeds.md)
  : Parallel kNN Accumulation with Explicit Seeds
- [`cpp_knn_single()`](https://gillescolling.com/spacc/reference/cpp_knn_single.md)
  : Single kNN Accumulation Curve
- [`cpp_phylo_knn_parallel()`](https://gillescolling.com/spacc/reference/cpp_phylo_knn_parallel.md)
  : Parallel Phylogenetic Diversity Accumulation
- [`cpp_phylo_knn_single()`](https://gillescolling.com/spacc/reference/cpp_phylo_knn_single.md)
  : Single kNN Phylogenetic Diversity Accumulation
- [`cpp_radius_parallel()`](https://gillescolling.com/spacc/reference/cpp_radius_parallel.md)
  : Parallel Radius-Order Accumulation
- [`cpp_radius_single()`](https://gillescolling.com/spacc/reference/cpp_radius_single.md)
  : Single Radius-Order Accumulation
- [`cpp_random_parallel()`](https://gillescolling.com/spacc/reference/cpp_random_parallel.md)
  : Parallel Random Accumulation
- [`cpp_random_single()`](https://gillescolling.com/spacc/reference/cpp_random_single.md)
  : Single Random Accumulation Curve
- [`cpp_wavefront_parallel()`](https://gillescolling.com/spacc/reference/cpp_wavefront_parallel.md)
  : Parallel Wavefront Accumulation
- [`cpp_wavefront_single()`](https://gillescolling.com/spacc/reference/cpp_wavefront_single.md)
  : Single Wavefront Accumulation
- [`calc_coverage()`](https://gillescolling.com/spacc/reference/calc_coverage.md)
  : Calculate Good-Turing Coverage Estimate
- [`calc_faith_pd()`](https://gillescolling.com/spacc/reference/calc_faith_pd.md)
  : Calculate Faith's PD from branch lengths
- [`calc_fdis()`](https://gillescolling.com/spacc/reference/calc_fdis.md)
  : Calculate Functional Dispersion (FDis)
- [`calc_fric_approx()`](https://gillescolling.com/spacc/reference/calc_fric_approx.md)
  : Calculate Functional Richness (convex hull volume approximation)
- [`calc_hill_number()`](https://gillescolling.com/spacc/reference/calc_hill_number.md)
  : Calculate Hill Number for a vector of abundances
- [`calc_mntd()`](https://gillescolling.com/spacc/reference/calc_mntd.md)
  : Calculate Mean Nearest Taxon Distance (MNTD)
- [`calc_mpd()`](https://gillescolling.com/spacc/reference/calc_mpd.md)
  : Calculate Mean Pairwise Distance (MPD)
- [`interpolate_at_coverage()`](https://gillescolling.com/spacc/reference/interpolate_at_coverage.md)
  : Interpolate Richness at Target Coverage Levels
- [`analytical`](https://gillescolling.com/spacc/reference/analytical.md)
  : Analytical Accumulation Methods
- [`c(`*`<spacc>`*`)`](https://gillescolling.com/spacc/reference/c.spacc.md)
  : Combine spacc Objects
- [`plot(`*`<spacc>`*`)`](https://gillescolling.com/spacc/reference/plot.spacc.md)
  : Plot Spatial SAC
- [`plot(`*`<spacc_metrics>`*`)`](https://gillescolling.com/spacc/reference/plot.spacc_metrics.md)
  : Plot spacc_metrics
