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

## Diversity

Alpha, beta, and gamma diversity partitioning

- [`alphaDiversity()`](https://gillescolling.com/spacc/reference/alphaDiversity.md)
  : Alpha Diversity (Per-Site)
- [`gammaDiversity()`](https://gillescolling.com/spacc/reference/gammaDiversity.md)
  : Gamma Diversity (Regional)
- [`diversityPartition()`](https://gillescolling.com/spacc/reference/diversityPartition.md)
  : Alpha-Beta-Gamma Diversity Partitioning

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

- [`analytical`](https://gillescolling.com/spacc/reference/analytical.md)
  : Analytical Accumulation Methods
- [`c(`*`<spacc>`*`)`](https://gillescolling.com/spacc/reference/c.spacc.md)
  : Combine spacc Objects
- [`plot(`*`<spacc>`*`)`](https://gillescolling.com/spacc/reference/plot.spacc.md)
  : Plot Spatial SAC
- [`plot(`*`<spacc_metrics>`*`)`](https://gillescolling.com/spacc/reference/plot.spacc_metrics.md)
  : Plot spacc_metrics
