# spacc API Design

> Spatial species accumulation curves. Fast.

## Core Philosophy

- [`spacc()`](https://gillescolling.com/spacc/reference/spacc.md) is the
  main entry point
- All objects have [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods
- Power users can access individual components
- C++ backend, R interface

------------------------------------------------------------------------

## Class Hierarchy

    spacc        <- main accumulation result
    spacc_dist   <- distance matrix with metadata
    spacc_fit    <- extrapolation model fit
    spacc_comp   <- comparison between two curves
    spacc_rare   <- rarefaction result

------------------------------------------------------------------------

## Exported Functions

### Main Entry Point

#### `spacc()`

Compute spatial species accumulation curves.

``` r

spacc(
  x,                          # species matrix OR formula
  coords,                     # coordinates (data.frame or matrix)
  data = NULL,                # optional data.frame for formula
  n_seeds = 50L,              # number of starting points
  method = "knn",             # "knn", "kncn", "random"
  distance = "euclidean",     # "euclidean", "haversine", or spacc_dist object
  parallel = TRUE,            # use parallel processing
  n_cores = NULL,             # NULL = detectCores() - 1
  progress = TRUE,            # show progress bar
  seed = NULL                 # RNG seed for reproducibility
)
```

**Returns:** Object of class `spacc`

------------------------------------------------------------------------

### Distance Matrix

#### `distances()`

Pre-compute distance matrix for reuse.

``` r

distances(
  coords,                     # coordinates
  method = "euclidean",       # "euclidean", "haversine", "custom"
  fun = NULL                  # custom distance function
)
```

**Returns:** Object of class `spacc_dist`

**Why?** Computing distances is O(n²). For large datasets, compute once
and reuse:

``` r

d <- distances(coords, method = "haversine")
sac_native <- spacc(native_species, d)
sac_alien <- spacc(alien_species, d)
```

------------------------------------------------------------------------

### Extrapolation

#### `extrapolate()`

Fit asymptotic model to estimate total species richness.

``` r

extrapolate(
  object,                     # spacc object
  model = "michaelis-menten", # "michaelis-menten", "lomolino", "asymptotic", "weibull"
  ...
)
```

**Returns:** Object of class `spacc_fit`

``` r

fit <- extrapolate(sac, model = "lomolino")
print(fit)
# Estimated asymptote: 342 species (95% CI: 318-371)
# Model: Lomolino (AIC: 234.5)

plot(fit)  # curve with extrapolation + CI
```

------------------------------------------------------------------------

### Comparison

#### `compare()`

Test if two accumulation curves differ significantly.

``` r

compare(
  x,                          # spacc object 1
  y,                          # spacc object 2
  method = "permutation",     # "permutation", "bootstrap", "auc"
  n_perm = 999L,              # permutations for p-value
  ...
)
```

**Returns:** Object of class `spacc_comp`

``` r

comp <- compare(sac_native, sac_alien)
print(comp)
# Comparison: native vs alien
# Method: permutation (n=999)
# AUC difference: 1234.5 (p < 0.001)
# Native saturates faster (89 vs 134 plots to 90%)

plot(comp)  # both curves with difference highlighted
```

------------------------------------------------------------------------

### Rarefaction

#### `rarefy()`

Individual-based rarefaction (classic method).

``` r

rarefy(
  x,                          # species matrix (with abundances)
  n_individuals = NULL,       # target sample size (NULL = all levels)
  n_boot = 100L,              # bootstrap replicates for CI
  ...
)
```

**Returns:** Object of class `spacc_rare`

``` r

rare <- rarefy(species_abundance)
plot(rare)
```

------------------------------------------------------------------------

### Subsampling

#### `subsample()`

Spatially stratified subsampling.

``` r

subsample(
  coords,                     # coordinates
  n = NULL,                   # target number of points
  method = "grid",            # "grid", "random", "cluster"
  cell_size = NULL,           # for grid method
  min_dist = NULL,            # minimum distance between points
  ...
)
```

**Returns:** Integer vector of indices to keep

``` r

# Thin to reduce spatial autocorrelation
keep <- subsample(coords, method = "grid", cell_size = 10)
sac <- spacc(species[keep, ], coords[keep, ])
```

------------------------------------------------------------------------

### Utility Functions

#### `as_spacc()`

Convert from other formats.

``` r

as_spacc(x, ...)

# Methods:
as_spacc.specaccum()         # from vegan
as_spacc.mob_out()           # from mobr
as_spacc.data.frame()        # from tidy format
```

#### `c.spacc()`

Combine multiple spacc objects (for multi-group comparisons).

``` r

all_curves <- c(native = sac_native, alien = sac_alien, ...)
plot(all_curves)
```

------------------------------------------------------------------------

## S3 Methods

### For class `spacc`

| Method | Description |
|----|----|
| `print.spacc()` | One-line summary |
| `summary.spacc()` | Detailed statistics |
| [`plot.spacc()`](https://gillescolling.com/spacc/reference/plot.spacc.md) | ggplot2 visualization |
| `as.data.frame.spacc()` | Convert to tidy format |
| `[.spacc()` | Subset seeds |
| [`c.spacc()`](https://gillescolling.com/spacc/reference/c.spacc.md) | Combine objects |

### For class `spacc_dist`

| Method                   | Description                   |
|--------------------------|-------------------------------|
| `print.spacc_dist()`     | Dimensions + method           |
| `plot.spacc_dist()`      | Distance histogram or heatmap |
| `as.matrix.spacc_dist()` | Extract raw matrix            |

### For class `spacc_fit`

| Method                | Description                 |
|-----------------------|-----------------------------|
| `print.spacc_fit()`   | Asymptote estimate + CI     |
| `summary.spacc_fit()` | Model diagnostics           |
| `plot.spacc_fit()`    | Curve + extrapolation       |
| `predict.spacc_fit()` | Predict at new sample sizes |
| `coef.spacc_fit()`    | Model coefficients          |
| `confint.spacc_fit()` | Confidence intervals        |

### For class `spacc_comp`

| Method                 | Description              |
|------------------------|--------------------------|
| `print.spacc_comp()`   | Test result              |
| `summary.spacc_comp()` | Detailed comparison      |
| `plot.spacc_comp()`    | Both curves + difference |

### For class `spacc_rare`

| Method               | Description       |
|----------------------|-------------------|
| `print.spacc_rare()` | Summary           |
| `plot.spacc_rare()`  | Rarefaction curve |

------------------------------------------------------------------------

## Plot Customization

All [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods
return ggplot2 objects, so users can customize:

``` r

plot(sac) +
  labs(title = "My Custom Title") +
  theme_dark()
```

### Plot options

``` r

plot.spacc(
  x,
  ci = TRUE,                  # show confidence interval
  ci_level = 0.95,            # CI level
  ci_alpha = 0.3,             # ribbon transparency
  show_seeds = FALSE,         # show individual curves
  saturation = TRUE,          # mark saturation point
  saturation_level = 0.9,     # 90% of max
  ...
)
```

------------------------------------------------------------------------

## Progress & Messaging

Use cli package for pretty output:

``` r

sac <- spacc(species, coords, n_seeds = 100)
# ℹ Computing distances (5000 × 5000)
# ℹ Running kNN accumulation
# ■■■■■■■■■■■■■■■■■■■■ 100% | ETA: 0s
# ✔ Done in 2.3s
```

------------------------------------------------------------------------

## Internal (Not Exported)

``` r

# C++ functions (prefixed with cpp_)
cpp_dist_euclidean()
cpp_dist_haversine()
cpp_knn_single()
cpp_knn_parallel()
cpp_kncn_single()
cpp_kncn_parallel()
cpp_random_single()
cpp_random_parallel()

# R helpers
validate_species()
validate_coords()
detect_cores()
make_progress_bar()
```

------------------------------------------------------------------------

## Example Workflows

### Basic usage

``` r

library(spacc)

sac <- spacc(species, coords)
print(sac)
summary(sac)
plot(sac)
```

### Compare groups

``` r

sac_native <- spacc(native_sp, coords)
sac_alien <- spacc(alien_sp, coords)

comp <- compare(sac_native, sac_alien)
print(comp)
plot(comp)
```

### Full analysis

``` r

# Pre-compute distances (lat/lon data)
d <- distances(coords, method = "haversine")

# Accumulation curves
sac <- spacc(species, d, n_seeds = 100)

# Extrapolate total richness
fit <- extrapolate(sac, model = "lomolino")
predict(fit, n = c(100, 500, 1000))

# Export for paper
ggsave("figure1.pdf", plot(sac))
write.csv(as.data.frame(sac), "results.csv")
```

### Large dataset

``` r

# Subsample to reduce computation
keep <- subsample(coords, method = "grid", cell_size = 1)

# Run with progress
sac <- spacc(species[keep, ], coords[keep, ],
             n_seeds = 200, progress = TRUE)
```

------------------------------------------------------------------------

## Dependencies

**Imports:** - Rcpp, RcppParallel (C++ backend) - ggplot2 (plotting) -
cli (progress bars, messages)

**Suggests:** - vegan (for as_spacc.specaccum) - mobr (for
as_spacc.mob_out) - sf (for spatial operations in subsample)
