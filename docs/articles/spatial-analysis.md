# Spatial Analysis: Endemism, Fragmentation, and SAR

## Overview

spacc includes tools for conservation-oriented spatial analyses:
endemism-area curves, fragmentation effects on species richness (SFAR),
and sampling-effort correction for species-area relationships (SESARS).

## Data Setup

``` r

library(spacc)

set.seed(42)
n_sites <- 80
n_species <- 40

coords <- data.frame(
  x = runif(n_sites, 0, 100),
  y = runif(n_sites, 0, 100)
)

# Species with varying range sizes (some endemic, some widespread)
species <- matrix(0L, n_sites, n_species)
for (sp in seq_len(n_species)) {
  cx <- runif(1, 10, 90)
  cy <- runif(1, 10, 90)
  # First 10 species are narrow-ranged (endemics)
  spread <- if (sp <= 10) 0.005 else 0.001
  prob <- exp(-spread * ((coords$x - cx)^2 + (coords$y - cy)^2))
  species[, sp] <- rbinom(n_sites, 1, prob)
}
colnames(species) <- paste0("sp", seq_len(n_species))
```

## Endemism-Area Curves

[`spaccEndemism()`](https://gillescolling.com/spacc/reference/spaccEndemism.md)
tracks how many species are endemic (restricted) to the set of sites
visited so far. As the sampled area expands, the number of endemic
species initially rises then falls as widespread species are
encountered:

``` r

end <- spaccEndemism(species, coords, n_seeds = 20, progress = FALSE)
end
#> spacc endemism: 80 sites, 40 species, 20 seeds
#> Final mean endemism: 40.0 species (100% of total)
```

``` r

plot(end)
```

![Endemism-area curve showing endemic species as area
expands.](spatial-analysis_files/figure-html/plot-endemism-1.svg)

Endemism-area curve showing endemic species as area expands.

The endemism curve is informative for identifying irreplaceable sites —
areas where unique species would be lost if not protected.

## Species-Fragmented Area Relationship (SFAR)

[`sfar()`](https://gillescolling.com/spacc/reference/sfar.md) fits the
SFAR model (Hanski et al. 2013), which quantifies how habitat
fragmentation (number of fragments) reduces species richness beyond what
area loss alone predicts:

``` math
\log(S) = c + z \cdot \log(A) + f \cdot \log(n)
```

where $`S`$ is species richness, $`A`$ is total area, $`n`$ is number of
fragments, and $`f`$ is the fragmentation exponent.

``` r

# First create a spacc object
sac <- spacc(species, coords, n_seeds = 20, progress = FALSE)

# Define patch assignments (e.g., from k-means clustering)
set.seed(123)
patches <- kmeans(coords, centers = 5)$cluster

# Fit SFAR model
sfar_fit <- sfar(sac, patches)
sfar_fit
#> SFAR: Species-Fragmented Area Relationship
#> -------------------------------------------- 
#> Fragments: 5
#> R-squared: 0.985
#> 
#> Model: S = 7.74 * A^0.300 * n^(--0.248)
#> Fragmentation effect (f): -0.248
#>   No additional fragmentation penalty detected
```

``` r

plot(sfar_fit)
```

![Species-Fragmented Area
Relationship.](spatial-analysis_files/figure-html/plot-sfar-1.svg)

Species-Fragmented Area Relationship.

A negative fragmentation exponent $`f`$ indicates that splitting habitat
into more fragments reduces species richness beyond the area effect.

## Sampling-Effort Corrected SAR (SESARS)

[`sesars()`](https://gillescolling.com/spacc/reference/sesars.md) fits a
species-area relationship corrected for variation in sampling effort
across sites, avoiding the bias that arises when species-poor areas are
also under-sampled:

``` math
\log(S) = c + z \cdot \log(A) + e \cdot \log(E)
```

``` r

# Simulate varying sampling effort per site
effort <- rpois(n_sites, 10) + 1

# Fit SESARS model using the spacc object from above
sesars_fit <- sesars(sac, effort)
sesars_fit
#> SESARS: Sampling Effort Species-Area Relationship
#> ------------------------------------------------ 
#> Model: power
#> R-squared: 0.996
#> 
#> Coefficients:
#>   log_c       z       w 
#>  4.2099  2.2586 -1.5486 
#> 
#> Interpretation: S = 67.35 * A^2.259 * E^-1.549
```

``` r

plot(sesars_fit)
```

![Sampling-effort corrected
SAR.](spatial-analysis_files/figure-html/plot-sesars-1.svg)

Sampling-effort corrected SAR.

## Per-Site Metrics for Spatial Prioritisation

[`spaccMetrics()`](https://gillescolling.com/spacc/reference/spaccMetrics.md)
computes accumulation metrics from every site as a starting point. This
reveals spatial variation in diversity accumulation rates, useful for
conservation planning:

``` r

met <- spaccMetrics(species, coords,
                    metrics = c("slope_10", "half_richness", "auc"),
                    progress = FALSE)
```

``` r

plot(met, metric = "slope_10", type = "heatmap")
#> Warning in viridisLite::viridis(n, alpha, begin, end, direction, option):
#> Option 'v' does not exist. Defaulting to 'viridis'.
```

![Spatial heatmap of initial accumulation
slope.](spatial-analysis_files/figure-html/plot-metrics-1.svg)

Spatial heatmap of initial accumulation slope.

Sites with steep initial slopes represent areas where many species are
concentrated — potential biodiversity hotspots.

## Spatial Subsampling

[`subsample()`](https://gillescolling.com/spacc/reference/subsample.md)
reduces spatial autocorrelation when preparing data for macroecological
analyses:

``` r

# Grid-based thinning
keep <- subsample(coords, method = "grid", cell_size = 20)
length(keep)
#> [1] 24

# Compare thinned vs full SAC
sac_full <- spacc(species, coords, n_seeds = 20, progress = FALSE)
sac_thin <- spacc(species[keep, ], coords[keep, ], n_seeds = 20, progress = FALSE)
combined <- c(full = sac_full, thinned = sac_thin)
```

``` r

plot(combined)
```

![Effect of spatial thinning on
SAC.](spatial-analysis_files/figure-html/plot-subsample-1.svg)

Effect of spatial thinning on SAC.

## References

- Hanski, I., Zurita, G.A., Bellocq, M.I. & Rybicki, J. (2013).
  Species-fragmented area relationship. Proceedings of the National
  Academy of Sciences, 110, 12715-12720.
- Aiello-Lammens, M.E., Boria, R.A., Radosavljevic, A., et al. (2015).
  spThin: an R package for spatial thinning. Ecography, 38, 541-545.
