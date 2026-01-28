# Getting Started with spacc

## Overview

**spacc** computes spatial species accumulation curves (SACs) using
nearest-neighbour algorithms with a C++ backend for speed. Unlike
classical random SACs, spatial SACs respect the geographic arrangement
of sampling sites, revealing how species richness accumulates as you
expand outward from a focal point.

## Simulating Example Data

``` r

library(spacc)

set.seed(42)
n_sites <- 80
n_species <- 40

coords <- data.frame(
  x = runif(n_sites, 0, 100),
  y = runif(n_sites, 0, 100)
)

# Simulate presence/absence with spatial structure
species <- matrix(0L, n_sites, n_species)
for (sp in seq_len(n_species)) {
  cx <- runif(1, 20, 80)
  cy <- runif(1, 20, 80)
  prob <- exp(-0.001 * ((coords$x - cx)^2 + (coords$y - cy)^2))
  species[, sp] <- rbinom(n_sites, 1, prob)
}
colnames(species) <- paste0("sp", seq_len(n_species))
```

## Basic Spatial Accumulation Curve

The core function
[`spacc()`](https://gillescolling.com/spacc/reference/spacc.md) computes
kNN-based accumulation curves from multiple random starting sites
(seeds):

``` r

sac <- spacc(species, coords, n_seeds = 30, method = "knn", progress = FALSE)
sac
#> spacc: 80 sites, 40 species, 30 seeds (knn)
```

``` r

plot(sac)
```

![Spatial species accumulation curve with 95% confidence
ribbon.](quickstart_files/figure-html/plot-sac-1.svg)

Spatial species accumulation curve with 95% confidence ribbon.

The grey ribbon shows variability across starting points — wider ribbons
indicate that species richness depends on *where* you start sampling.

## Comparing Methods

spacc supports several accumulation methods:

- **knn**: k-Nearest Neighbour (default, fast)

- **kncn**: k-Nearest Centroid Neighbour (geometrically centred
  expansion)

- **random**: Classical random accumulation (no spatial structure)

``` r

sac_kncn <- spacc(species, coords, n_seeds = 30, method = "kncn", progress = FALSE)
sac_rand <- spacc(species, coords, n_seeds = 30, method = "random", progress = FALSE)
```

Combine curves for comparison:

``` r

combined <- c(knn = sac, kncn = sac_kncn, random = sac_rand)
plot(combined)
```

![Comparison of accumulation
methods.](quickstart_files/figure-html/combine-1.svg)

Comparison of accumulation methods.

## Extrapolation

Fit an asymptotic model to estimate total species richness:

``` r

fit <- extrapolate(sac, model = "lomolino")
#> Waiting for profiling to be done...
fit
#> Extrapolation: lomolino 
#> ------------------------------ 
#> Estimated asymptote: 42.7 species
#> 95% CI: 42.2 - 43.3
#> AIC: 165.3
#> Observed: 40.0 species (94% of estimated)
```

``` r

plot(fit)
```

![Fitted Lomolino curve with asymptote
estimate.](quickstart_files/figure-html/plot-fit-1.svg)

Fitted Lomolino curve with asymptote estimate.

## Pre-computing Distances

For repeated analyses on the same sites, pre-compute the distance matrix
once:

``` r

d <- distances(coords, method = "euclidean")

sac1 <- spacc(species[, 1:20], d, n_seeds = 30, progress = FALSE)
sac2 <- spacc(species[, 21:40], d, n_seeds = 30, progress = FALSE)
```

## Comparing Curves

Test whether two accumulation curves differ significantly:

``` r

comp <- compare(sac1, sac2, method = "permutation", n_perm = 199)
comp
#> Comparison: sac1 vs sac2 
#> ---------------------------------------- 
#> Method: permutation (n=199)
#> AUC difference: -8.6 (p = 0.573)
#> Saturation: sac1 at 23 sites, sac2 at 19 sites
```

## Next Steps

- [Diversity](https://gillescolling.com/spacc/articles/diversity.md):
  Hill numbers, beta diversity, phylogenetic and functional diversity
  accumulation
- [Extrapolation](https://gillescolling.com/spacc/articles/extrapolation.md):
  Coverage-based extrapolation, EVT models, diversity-area relationships
- [Spatial
  Analysis](https://gillescolling.com/spacc/articles/spatial-analysis.md):
  Endemism curves, fragmentation analysis, sampling-effort correction
