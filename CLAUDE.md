# spacc Development Guide

## Package Overview

Fast spatial species accumulation curves with C++ backend (Rcpp/RcppParallel).

## Build Commands

```r
# Load for development
devtools::load_all()

# Document
devtools::document()

# Check
devtools::check()

# Build pkgdown site
source("~/.R/build_pkgdown.R"); build_pkgdown_site()
```

## Project Structure

```
spacc/
├── R/                  # R source files
│   ├── spacc.R        # Main spacc(), wavefront(), distanceDecay()
│   ├── analytical.R   # coleman(), mao_tau(), spatialRarefaction()
│   ├── extrapolate.R  # extrapolate() asymptotic models
│   ├── compare.R      # compare() curve comparisons
│   ├── rarefy.R       # rarefy() rarefaction
│   ├── subsample.R    # subsample() spatial subsampling
│   ├── distances.R    # distances() pre-computation
│   ├── methods.R      # S3 methods (print, summary, plot)
│   ├── plot.R         # ggplot2-based plotting
│   └── zzz.R          # Package startup
├── src/               # C++ source files
│   ├── distance.cpp   # Distance matrix computation
│   ├── knn.cpp        # k-Nearest Neighbor accumulation
│   ├── kncn.cpp       # k-Nearest Centroid Neighbor
│   ├── random.cpp     # Random accumulation
│   └── methods.cpp    # Gaussian, cone, wavefront, distance-decay
└── tests/testthat/    # Unit tests
```

## Naming Conventions

- **Functions**: camelCase (`distanceDecay`, `spatialRarefaction`)
- **Arguments**: snake_case (`n_seeds`, `n_cores`)
- **Classes**: lowercase with underscore (`spacc`, `spacc_fit`, `spacc_decay`)
- **C++ functions**: snake_case with `cpp_` prefix (`cpp_knn_parallel`)

## Key Classes

- `spacc` - Main accumulation curve object
- `spacc_fit` - Fitted asymptotic model
- `spacc_comp` - Comparison result
- `spacc_rare` - Rarefaction result
- `spacc_dist` - Pre-computed distance matrix
- `spacc_wavefront` - Wavefront expansion result
- `spacc_decay` - Distance-decay result

## Testing

```r
devtools::test()
```

## C++ Development

After modifying C++ code:
```r
Rcpp::compileAttributes()
devtools::load_all()
```
