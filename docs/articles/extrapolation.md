# Extrapolation and Species-Area Models

## Overview

spacc provides tools for estimating total species richness beyond
observed sampling effort: asymptotic model fitting, coverage-based
extrapolation, and diversity-area relationships (DAR).

## Data Setup

``` r

library(spacc)

set.seed(42)
n_sites <- 100
n_species <- 50

coords <- data.frame(
  x = runif(n_sites, 0, 100),
  y = runif(n_sites, 0, 100)
)

species <- matrix(0L, n_sites, n_species)
for (sp in seq_len(n_species)) {
  cx <- runif(1, 10, 90)
  cy <- runif(1, 10, 90)
  lambda <- 4 * exp(-0.0008 * ((coords$x - cx)^2 + (coords$y - cy)^2))
  species[, sp] <- rpois(n_sites, lambda)
}
colnames(species) <- paste0("sp", seq_len(n_species))
pa <- (species > 0) * 1L
```

## Asymptotic Extrapolation

[`extrapolate()`](https://gillescolling.com/spacc/reference/extrapolate.md)
fits an asymptotic model to the mean accumulation curve to estimate
total species richness:

``` r

sac <- spacc(pa, coords, n_seeds = 30, progress = FALSE)
```

### Model Comparison

spacc supports six asymptotic models. Compare them to select the best
fit:

``` r

models <- c("michaelis-menten", "lomolino", "asymptotic", "weibull", "logistic")
fits <- lapply(models, function(m) extrapolate(sac, model = m))
#> Waiting for profiling to be done...
#> Waiting for profiling to be done...
#> Waiting for profiling to be done...
#> Warning in value[[3L]](cond): Model fitting failed: Missing value or an
#> infinity produced when evaluating the model
#> Waiting for profiling to be done...
names(fits) <- models

# Compare AIC
data.frame(
  model = models,
  asymptote = sapply(fits, function(f) round(f$asymptote, 1)),
  AIC = sapply(fits, function(f) round(f$aic, 1))
)
#>                               model asymptote   AIC
#> michaelis-menten.a michaelis-menten      51.8 271.8
#> lomolino.a                 lomolino      51.8 273.8
#> asymptotic.a             asymptotic      49.6 390.8
#> weibull                     weibull        NA    NA
#> logistic.a                 logistic      49.9 157.3
```

``` r

best <- fits[[which.min(sapply(fits, function(f) f$aic))]]
plot(best)
```

![Best-fitting asymptotic
model.](extrapolation_files/figure-html/plot-best-1.svg)

Best-fitting asymptotic model.

### EVT Model

The Extreme Value Theory model (Borda-de-Agua et al. 2025) uses a
two-component Weibull mixture to capture the triphasic pattern observed
in empirical species-area relationships:

``` r

fit_evt <- extrapolate(sac, model = "evt")
#> Warning in stats::nls(y ~ a * (w * (1 - exp(-(x/b1)^c1)) + (1 - w) * (1 - :
#> Convergence failure: singular convergence (7)
#> Waiting for profiling to be done...
#> Warning in regularize.values(x, y, ties, missing(ties), na.rm = na.rm):
#> collapsing to unique 'x' values
fit_evt
#> Extrapolation: evt 
#> ------------------------------ 
#> Estimated asymptote: 50.0 species
#> 95% CI: NA - NA
#> AIC: 50.9
#> Observed: 50.0 species (100% of estimated)
```

The EVT model requires sufficient data points (\>10) and works best with
real ecological data exhibiting a clear triphasic SAR pattern.

## Coverage-Based Extrapolation

[`extrapolateCoverage()`](https://gillescolling.com/spacc/reference/extrapolateCoverage.md)
estimates richness at coverage levels beyond the observed maximum, using
asymptotic estimators (Chao et al. 2014):

``` r

cov <- spaccCoverage(species, coords, n_seeds = 20, progress = FALSE)

ext <- extrapolateCoverage(cov, target_coverage = c(0.95, 0.99, 1.0), q = 0)
ext
#> Coverage-based extrapolation
#> -------------------------------- 
#> Diversity order: q = 0
#> Observed coverage: 100.0%
#> Observed richness: 50.0
#> 
#> Extrapolated richness:
#>   C=95%: 36.6 (+/- 5.9)
#>   C=99%: 42.4 (+/- 4.0)
#>   C=100%: 48.3 (+/- 2.5)
```

``` r

plot(ext)
```

![Coverage-based extrapolation of species
richness.](extrapolation_files/figure-html/plot-cov-extrap-1.svg)

Coverage-based extrapolation of species richness.

This is particularly useful when comparing communities with different
abundances — standardizing by coverage ensures equal completeness rather
than equal sample size.

## Diversity-Area Relationship (DAR)

[`dar()`](https://gillescolling.com/spacc/reference/dar.md) fits the
relationship between diversity (Hill numbers) and area (Ma 2018). This
generalises the classical species-area relationship to all Hill number
orders:

``` r

dar_result <- dar(species, coords, q = c(0, 1, 2), n_seeds = 20, progress = FALSE)
dar_result
#> spacc DAR: 100 sites, 20 seeds
#> Orders (q): 0, 1, 2
#> Area method: convex_hull
```

``` r

plot(dar_result)
```

![Diversity-area relationship for q = 0, 1,
2.](extrapolation_files/figure-html/plot-dar-1.svg)

Diversity-area relationship for q = 0, 1, 2.

The DAR captures how diversity scales with area. For q=0 this reduces to
the classical SAR; for higher q it quantifies how evenness scales
spatially.

## Predicting Richness at New Effort Levels

Use [`predict()`](https://rdrr.io/r/stats/predict.html) on a fitted
extrapolation model:

``` r

predict(best, n = c(50, 100, 200, 500))
#> [1] 49.87375 49.87643 49.87643 49.87643
```

## References

- Borda-de-Agua, L., Whittaker, R.J., Cardoso, P., et al. (2025).
  Extreme value theory explains the topography and scaling of the
  species-area relationship. Nature Communications, 16, 5346.
- Chao, A., Gotelli, N.J., Hsieh, T.C., et al. (2014). Rarefaction and
  extrapolation with Hill numbers. Ecological Monographs, 84, 45-67.
- Lomolino, M.V. (2000). Ecology’s most general, yet protean pattern:
  the species-area relationship. Journal of Biogeography, 27, 17-26.
- Ma, Z. (2018). Generalizing the species-area relationship with
  diversity-area relationships. Ecology and Evolution, 8, 8645-8655.
