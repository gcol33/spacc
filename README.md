# spacc
[![R-CMD-check](https://github.com/gcol33/spacc/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/gcol33/spacc/actions/workflows/R-CMD-check.yml)
[![codecov](https://codecov.io/gh/gcol33/spacc/branch/main/graph/badge.svg)](https://app.codecov.io/gh/gcol33/spacc)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
**Fast Spatial Species Accumulation Curves with C++ Performance**
The `spacc` package provides high-performance functions for computing **spatially-explicit species accumulation curves** and related diversity metrics. Unlike traditional accumulation curves that assume random sampling, `spacc` respects the spatial arrangement of sampling sites, providing ecologically meaningful estimates for biodiversity assessment, invasion modeling, and community ecology.
## Quick Start
```r
library(spacc)
coords <- data.frame(x = runif(100), y = runif(100))
species <- matrix(rpois(100 * 50, 2), nrow = 100)
# Spatial species accumulation curve
sac <- spacc(species, coords, method = "knn", n_seeds = 50)
plot(sac)
# Alpha-beta-gamma diversity partitioning
part <- diversityPartition(species, q = c(0, 1, 2))
print(part)
```
## Statement of Need
Species accumulation curves (SACs) are fundamental tools in ecology for understanding sampling completeness and comparing biodiversity. Standard SACs assume random sampling order, but real sampling is spatially structured.
**spacc** addresses this by computing SACs that respect spatial arrangement:
- **Biodiversity assessment**: Realistic sampling completeness for spatially structured surveys
- **Community ecology**: Understand how diversity accumulates across space
- **Invasion modeling**: Model spread patterns with wavefront and directional methods
- **Conservation planning**: Identify areas with high diversity accumulation rates
## Features
### Core Accumulation Functions
| Function | Description |
|----------|-------------|
| `spacc()` | Main spatial accumulation with 7 methods |
| `wavefront()` | Expanding radius accumulation (invasion modeling) |
| `distanceDecay()` | Species richness vs distance from focal points |
**Spatial methods:** `knn`, `kncn`, `random`, `radius`, `gaussian`, `cone`, `collector`
### Alpha-Beta-Gamma Diversity
| Function | Description |
|----------|-------------|
| `alphaDiversity()` | Per-site Hill numbers (local diversity) |
| `gammaDiversity()` | Pooled Hill numbers (regional diversity) |
| `diversityPartition()` | Full alpha-beta-gamma decomposition (Jost 2007) |
### Spatial Diversity Accumulation
| Function | Description |
|----------|-------------|
| `spaccHill()` | Hill numbers across spatial accumulation |
| `spaccBeta()` | Beta diversity with turnover/nestedness partitioning |
| `spaccCoverage()` | Coverage-based rarefaction (Chao & Jost 2012) |
| `spaccPhylo()` | Phylogenetic diversity (MPD, MNTD) |
| `spaccFunc()` | Functional diversity (FDis, FRic) |
| `spaccMetrics()` | Per-site metrics with spatial mapping |
### Analytical Methods (No Simulation)
| Function | Description |
|----------|-------------|
| `coleman()` | Coleman expected accumulation |
| `mao_tau()` | Mao Tau expected curve |
| `spatialRarefaction()` | Spatially-constrained rarefaction |
### Analysis & Post-Processing
| Function | Description |
|----------|-------------|
| `extrapolate()` | Fit asymptotic models (5 options) |
| `compare()` | Statistical comparison between curves |
| `rarefy()` | Rarefaction to common sampling effort |
| `subsample()` | Spatial subsampling of sites |
| `distances()` | Pre-compute distance matrices |
## Installation
```r
# install.packages("pak")
pak::pak("gcol33/spacc")
```
## Usage Examples
### Alpha-Beta-Gamma Partitioning
```r
species <- matrix(rpois(100 * 50, 2), nrow = 100)
# Full diversity partitioning (Jost 2007)
part <- diversityPartition(species, q = c(0, 1, 2))
part
#> Alpha-Beta-Gamma Diversity Partitioning
#> 100 sites, 50 species
#>
#>  q alpha  beta gamma
#>  0 12.45  4.01    50
#>  1  8.23  3.12 25.68
#>  2  6.54  2.87 18.77
#>
#> Interpretation:
#>   Alpha = mean effective species per site
#>   Beta  = effective number of communities (1 to n_sites)
#>   Gamma = regional effective species (gamma = alpha x beta)
# Individual components
alpha <- alphaDiversity(species, q = c(0, 1, 2))  # per-site
gamma <- gammaDiversity(species, q = c(0, 1, 2))  # regional total
```
### Basic Spatial SAC
```r
coords <- data.frame(x = runif(100), y = runif(100))
species <- matrix(rbinom(100 * 50, 1, 0.3), nrow = 100)
sac <- spacc(species, coords, method = "knn", n_seeds = 100)
plot(sac)
# Compare to random accumulation
sac_random <- spacc(species, coords, method = "random", n_seeds = 100)
comp <- compare(sac, sac_random)
plot(comp)
```
### Hill Number Accumulation
```r
hill <- spaccHill(species, coords, q = c(0, 1, 2), n_seeds = 50)
plot(hill)
```
### Beta Diversity Accumulation
```r
beta <- spaccBeta(species, coords, index = "sorensen", n_seeds = 50)
plot(beta, partition = TRUE)  # Shows turnover vs nestedness
```
### Coverage-Based Rarefaction
```r
cov <- spaccCoverage(species, coords, n_seeds = 50)
interp <- interpolateCoverage(cov, target = c(0.90, 0.95))
plot(cov, xaxis = "coverage")
```
### Phylogenetic & Functional Diversity
```r
# Phylogenetic
library(ape)
tree <- rtree(50)
colnames(species) <- tree$tip.label
phylo <- spaccPhylo(species, coords, tree, metric = c("mpd", "mntd"))
# Functional
traits <- matrix(rnorm(50 * 3), nrow = 50)
rownames(traits) <- colnames(species)
func <- spaccFunc(species, coords, traits, metric = c("fdis", "fric"))
```
### Per-Site Metrics
```r
metrics <- spaccMetrics(species, coords, metrics = c("slope_10", "half_richness", "auc"))
plot(metrics, metric = "slope_10", type = "heatmap")
```
## Choosing the Right Function
| Goal | Function |
|------|----------|
| Alpha-beta-gamma partitioning | `diversityPartition()` |
| Per-site diversity | `alphaDiversity()` |
| Regional diversity | `gammaDiversity()` |
| Spatial accumulation | `spacc()` |
| Hill numbers over space | `spaccHill()` |
| Turnover patterns | `spaccBeta()` |
| Coverage standardization | `spaccCoverage()` |
| Phylogenetic diversity | `spaccPhylo()` |
| Functional diversity | `spaccFunc()` |
| Map diversity patterns | `spaccMetrics()` |
| Invasion modeling | `wavefront()` |
## Performance
C++ backend via Rcpp and RcppParallel:
| Sites | Species | Seeds | Time |
|-------|---------|-------|------|
| 100   | 50      | 100   | 0.02s |
| 500   | 100     | 100   | 0.3s |
| 1000  | 200     | 100   | 1.2s |
## Documentation
- [Function Reference](https://gillescolling.com/spacc/reference/)
## Support
> "Software is like sex: it is better when it is free." - Linus Torvalds
[![Buy Me A Coffee](https://img.shields.io/badge/-Buy%20me%20a%20coffee-FFDD00?logo=buymeacoffee&logoColor=black)](https://buymeacoffee.com/gcol33)
## License
MIT
## Citation
```bibtex
@software{spacc,
  author = {Colling, Gilles},
  title = {spacc: Fast Spatial Species Accumulation Curves},
  year = {2025},
  url = {https://github.com/gcol33/spacc}
}
```
