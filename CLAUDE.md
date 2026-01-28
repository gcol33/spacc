# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

Fast spatial species accumulation curves with C++ backend (Rcpp/RcppParallel). Supports Hill numbers, beta diversity, coverage-based rarefaction, phylogenetic/functional diversity, diversity-area relationships, and endemism curves.

## Build Commands

```r
devtools::load_all()                                    # Load for development
devtools::document()                                    # Regenerate man/ and NAMESPACE
devtools::test()                                        # Run tests
devtools::test(filter = "spacc")                        # Run single test file
devtools::check()                                       # Full R CMD check
source("~/.R/build_pkgdown.R"); build_pkgdown_site()    # Build pkgdown site
```

## C++ Development

After modifying any file in `src/`:
```r
Rcpp::compileAttributes()   # Regenerates R/RcppExports.R and src/RcppExports.cpp
devtools::document()        # Regenerates man pages
devtools::load_all()        # Recompile and load
```

- Only add `//` roxygen comments to C++ functions that are `@export`ed in NAMESPACE.
- All internal C++ wrapper functions must be listed in `utils::globalVariables()` in `R/zzz.R` — a partial list causes unlisted ones to fail R CMD check.
- In tests, call internal C++ wrappers with `spacc:::func()`, not bare `func()`.

## Architecture

### Pipeline Design

```
spacc(x, coords)  ─┬─> extrapolate()      -> spacc_fit      (asymptotic models)
                    ├─> compare()           -> spacc_comp     (curve comparison)
                    ├─> rarefy()            -> spacc_rare     (individual-based)
                    ├─> spaccHill()         -> spacc_hill     (Hill numbers q=0,1,2)
                    ├─> spaccBeta()         -> spacc_beta     (turnover/nestedness)
                    ├─> spaccCoverage()     -> spacc_coverage (Chao sample coverage)
                    ├─> spaccMetrics()      -> spacc_metrics  (per-site metrics)
                    ├─> spaccPhylo()        -> spacc_phylo    (Faith PD, MPD, MNTD)
                    ├─> spaccFunc()         -> spacc_func     (FDis, FRic)
                    ├─> wavefront()         -> spacc_wavefront
                    └─> distanceDecay()     -> spacc_decay
```

All downstream functions accept a `spacc` object and return typed S3 objects with `print()`, `summary()`, `plot()`, `as.data.frame()` methods.

### Two-Tier C++ Backend

- **Exact (brute-force)**: Uses precomputed distance matrix. Auto-selected for ≤500 sites.
- **KD-tree/Ball-tree**: Spatial indexing, no distance matrix needed. Auto-selected for >500 sites.
- Overridable via `backend` argument. Spatiotemporal distance always uses exact backend.

### Grouped Objects

Both `groups` (species split) and `support = "auto"` (country split) return grouped `spacc` objects where `curves` is a list-of-matrices instead of a single matrix. The function recursively calls itself per group, then combines results.

### Parallelism

Each seed runs independently via RcppParallel (no shared state). Multiple seeds produce a matrix of curves used for percentile-based confidence intervals.

## Naming Conventions

- **Functions**: camelCase (`distanceDecay`, `spatialRarefaction`)
- **Arguments**: snake_case (`n_seeds`, `n_cores`)
- **Classes**: lowercase with underscore (`spacc`, `spacc_fit`, `spacc_decay`)
- **C++ functions**: snake_case with `cpp_` prefix (`cpp_knn_parallel`)
