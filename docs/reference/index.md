# Package index

## Core Functions

Main spatial accumulation functions

- [`spacc()`](https://gillescolling.com/spacc/reference/spacc.md) :
  Spatial Species Accumulation Curves
- [`spaccDiversity()`](https://gillescolling.com/spacc/reference/spaccDiversity.md)
  : Spatial Accumulation of a Custom Diversity Metric
- [`wavefront()`](https://gillescolling.com/spacc/reference/wavefront.md)
  : Wavefront Expansion Accumulation
- [`distanceDecay()`](https://gillescolling.com/spacc/reference/distanceDecay.md)
  : Distance-Decay Analysis
- [`distances()`](https://gillescolling.com/spacc/reference/distances.md)
  : Compute Distance Matrix

## Diversity

Alpha, beta, and gamma diversity partitioning

- [`alphaDiversity()`](https://gillescolling.com/spacc/reference/alphaDiversity.md)
  : Alpha Diversity (Per-Site)
- [`gammaDiversity()`](https://gillescolling.com/spacc/reference/gammaDiversity.md)
  : Gamma Diversity (Regional)
- [`diversityPartition()`](https://gillescolling.com/spacc/reference/diversityPartition.md)
  : Alpha-Beta-Gamma Diversity Partitioning
- [`diversityProfile()`](https://gillescolling.com/spacc/reference/diversityProfile.md)
  : Diversity Profile
- [`diversityProfilePhylo()`](https://gillescolling.com/spacc/reference/diversityProfilePhylo.md)
  : Phylogenetic Diversity Profile
- [`diversityProfileFunc()`](https://gillescolling.com/spacc/reference/diversityProfileFunc.md)
  : Functional Diversity Profile
- [`evenness()`](https://gillescolling.com/spacc/reference/evenness.md)
  : Evenness Profiles

## Hill Numbers

Diversity accumulation with Hill numbers (q=0,1,2)

- [`spaccHill()`](https://gillescolling.com/spacc/reference/spaccHill.md)
  : Spatial Accumulation with Hill Numbers
- [`spaccHillCoverage()`](https://gillescolling.com/spacc/reference/spaccHillCoverage.md)
  : Spatial Hill Numbers at Standardized Coverage
- [`spaccHillBeta()`](https://gillescolling.com/spacc/reference/spaccHillBeta.md)
  : Spatial Hill Number Beta Diversity

## Beta Diversity

Spatial beta diversity with turnover/nestedness partitioning

- [`spaccBeta()`](https://gillescolling.com/spacc/reference/spaccBeta.md)
  : Spatial Beta Diversity Accumulation
- [`spaccBetaFunc()`](https://gillescolling.com/spacc/reference/spaccBetaFunc.md)
  : Functional Beta Diversity Accumulation
- [`spaccBetaPhylo()`](https://gillescolling.com/spacc/reference/spaccBetaPhylo.md)
  : Phylogenetic Beta Diversity Accumulation

## Community Turnover

Distance-decay, zeta diversity, and null models

- [`betaDecay()`](https://gillescolling.com/spacc/reference/betaDecay.md)
  : Beta Distance-Decay
- [`zetaDiversity()`](https://gillescolling.com/spacc/reference/zetaDiversity.md)
  : Zeta Diversity
- [`ses()`](https://gillescolling.com/spacc/reference/ses.md) :
  Standardized Effect Size (SES) via Null Models

## Phylogenetic & Functional Diversity

Accumulation of PD and FD metrics

- [`spaccPhylo()`](https://gillescolling.com/spacc/reference/spaccPhylo.md)
  : Spatial Phylogenetic Diversity Accumulation
- [`spaccFunc()`](https://gillescolling.com/spacc/reference/spaccFunc.md)
  : Spatial Functional Diversity Accumulation

## Coverage-Based Rarefaction

Standardize by sample completeness

- [`spaccCoverage()`](https://gillescolling.com/spacc/reference/spaccCoverage.md)
  : Coverage-Based Spatial Rarefaction
- [`interpolateCoverage()`](https://gillescolling.com/spacc/reference/interpolateCoverage.md)
  : Interpolate Richness at Target Coverage Levels
- [`extrapolateCoverage()`](https://gillescolling.com/spacc/reference/extrapolateCoverage.md)
  : Extrapolate Richness Beyond Observed Coverage

## Per-Site Metrics

Extract metrics per starting site for heatmap visualization

- [`spaccMetrics()`](https://gillescolling.com/spacc/reference/spaccMetrics.md)
  : Per-Site Accumulation Metrics
- [`as_sf()`](https://gillescolling.com/spacc/reference/as_sf.md) :
  Convert spacc_metrics to sf

## Richness Estimators

Non-parametric richness estimation

- [`chao1()`](https://gillescolling.com/spacc/reference/chao1.md) :
  Chao1 Richness Estimator
- [`chao2()`](https://gillescolling.com/spacc/reference/chao2.md) :
  Chao2 Richness Estimator
- [`iChao1()`](https://gillescolling.com/spacc/reference/iChao1.md) :
  Improved Chao1 (iChao1) Richness Estimator
- [`iChao2()`](https://gillescolling.com/spacc/reference/iChao2.md) :
  Improved Chao2 (iChao2) Richness Estimator
- [`ace()`](https://gillescolling.com/spacc/reference/ace.md) : ACE
  Richness Estimator
- [`jackknife()`](https://gillescolling.com/spacc/reference/jackknife.md)
  : Jackknife Richness Estimator
- [`bootstrap_richness()`](https://gillescolling.com/spacc/reference/bootstrap_richness.md)
  : Bootstrap Richness Estimator
- [`completenessProfile()`](https://gillescolling.com/spacc/reference/completenessProfile.md)
  : Sample Completeness Profile

## Spatial Structure

Spatial eigenvectors and variance partitioning

- [`spatialEigenvectors()`](https://gillescolling.com/spacc/reference/spatialEigenvectors.md)
  : Spatial Eigenvectors (PCNM/dbMEM)
- [`spatialPartition()`](https://gillescolling.com/spacc/reference/spatialPartition.md)
  : Spatial Variance Partitioning with MEMs

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

## Species-Area Relationships

Diversity-area relationships, endemism, fragmentation

- [`dar()`](https://gillescolling.com/spacc/reference/dar.md) :
  Diversity-Area Relationship (DAR)
- [`spaccEndemism()`](https://gillescolling.com/spacc/reference/spaccEndemism.md)
  : Spatial Endemism Accumulation
- [`sesars()`](https://gillescolling.com/spacc/reference/sesars.md) :
  Sampling Effort Species-Area Relationship (SESARS)
- [`sfar()`](https://gillescolling.com/spacc/reference/sfar.md) :
  Species-Fragmented Area Relationship (SFAR)

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
- [`c(`*`<spacc_beta>`*`)`](https://gillescolling.com/spacc/reference/c.spacc_beta.md)
  : Combine spacc_beta Objects
- [`c(`*`<spacc_coverage>`*`)`](https://gillescolling.com/spacc/reference/c.spacc_coverage.md)
  : Combine spacc_coverage Objects
- [`c(`*`<spacc_hill>`*`)`](https://gillescolling.com/spacc/reference/c.spacc_hill.md)
  : Combine spacc_hill Objects
- [`compareModels()`](https://gillescolling.com/spacc/reference/compareModels.md)
  : Compare Multiple SAR Models
- [`predict(`*`<spacc>`*`)`](https://gillescolling.com/spacc/reference/predict.spacc.md)
  : Predict from Empirical Accumulation Curve
- [`plot(`*`<spacc>`*`)`](https://gillescolling.com/spacc/reference/plot.spacc.md)
  : Plot Spatial SAC
- [`plot(`*`<spacc_metrics>`*`)`](https://gillescolling.com/spacc/reference/plot.spacc_metrics.md)
  : Plot spacc_metrics
- [`autoplot.spacc`](https://gillescolling.com/spacc/reference/autoplot.spacc.md)
  : Autoplot Methods for spacc Objects
