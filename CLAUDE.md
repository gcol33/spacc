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

### C++ Testing (Catch2)

Pure C++ algorithms are in `src/core/` headers and can be tested independently of R:

```bash
cd src/tests
./run_tests.bat          # Windows: compile and run tests
make all                 # Unix: compile and run tests
./run_coverage.bat       # Windows: run with gcov coverage
make coverage            # Unix: run with gcov coverage
```

Core headers (`src/core/`), pure C++ written against accessors (lambdas), used
by both the Rcpp layer and the Catch2 tests:
- `types.h` - Common type definitions
- `distance_core.h` - Euclidean/Haversine distance, `PointDistance`
- `accumulation_core.h` - Order generators (`nn_walk_order`, `kncn_order`,
  `focal_order`, `gaussian_order`, `cone_order`) and the accumulation kernels
  (`richness_along_order`, `accumulate_along_order`, `richness_within_radius`)
- `beta_core.h` - Baselga beta diversity partitioning
- `hill_core.h` - Hill numbers
- `coverage_core.h` - Chao-Jost and Chiu sample coverage
- `inext_core.h` - size- and coverage-based rarefaction/extrapolation of Hill
  numbers for one abundance sample (iNEXT.3D estimators)
- `phylo_core.h` - MPD, MNTD, Rao's Q, Faith's PD (`PhyloTree`)
- `func_core.h` - FDis, range-volume FRic, trait Rao's Q
- `arc_core.h` - accumulation rate curve richness estimator (Shestopaloff
  2025): rate curve, the three decay fits and their tails, the capped
  weights by closed-form projected gradient flow from zero, smoothing, ridge
  NNLS for the Poisson mixture

Rcpp layer (`src/*.cpp`):
- `orders.cpp` - one exported order generator per method, all run through
  the `OrderWorker` template (parallel over seeds; the tree variants build
  one spatial index per worker chunk)
- `accumulate.cpp` - `cpp_order_parallel()` (richness along orders) and
  `cpp_radius_parallel()` (wavefront / distance decay)
- `hill.cpp`, `coverage.cpp`, `beta.cpp`, `phylofunc.cpp` - a metric struct
  with `scratch()` and `at_step()` each, driven by `OrderMetricWorker`
  (`parallel_utils.h`). Adding a metric = one struct + one export.
- `distance.cpp` - distance matrix export
- `arc.cpp` - `arc()` point estimate, bootstrap refits and Poisson-mixture
  penalty fits (parallel over draws / penalties)

Adding an accumulation method: write the order generator in
`accumulation_core.h`, export it from `orders.cpp`, add its branch to
`.accumulation_orders()` in `R/spatial-orders.R`. Nothing else changes.

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
                    ├─> spaccWavefront()    -> spacc_wavefront
                    ├─> distanceDecay()     -> spacc_decay    (distance decay)
                    ├─> betaDecay()         -> spacc_beta_decay
                    ├─> zetaDiversity()     -> spacc_zeta     (zeta diversity)
                    └─> dar()               -> spacc_dar      (diversity-area)
```

All downstream functions accept a `spacc` object and return typed S3 objects with `print()`, `summary()`, `plot()`, `as.data.frame()` methods.

### Orders First

Every method is a rule for ordering sites. `.accumulation_orders()`
(`R/spatial-orders.R`) is the only place that turns a method name into an
`n_orders x n_sites` matrix of 1-based site indices; every front door then
hands those orders to a C++ kernel (`cpp_order_parallel()` for richness, the
`cpp_*_order_*_parallel()` workers for metrics). `all_sites = TRUE` gives the
per-site map variants (every site once as seed), `seed_pool` restricts seeds
to core sites under `support`.

Shared S3 helpers live in `R/curves-utils.R`: `.resolve_coords()`,
`.curve_stats()`, `.curve_frame()`, `.print_curves()`, `.plot_curves_by()`,
`.plot_curve()`, `.plot_site_map()`, `.as_sf_site_values()`, `.stack_seeds()`.
Use them instead of re-deriving quantiles or ribbons in a new class.

### Two-Tier C++ Backend

- **Exact (brute-force)**: Uses precomputed distance matrix. Auto-selected for ≤500 sites.
- **KD-tree/Ball-tree**: Spatial indexing, no distance matrix needed. Auto-selected for >500 sites.
- Overridable via `backend` argument. A supplied distance matrix (`spacc_dist`
  input, spatiotemporal composite) always takes the exact backend.

### Grouped Objects

Both `groups` (species split) and `support = "auto"` (country split) return grouped `spacc` objects where `curves` is a list-of-matrices instead of a single matrix. The function recursively calls itself per group, then combines results.

### Parallelism

Each seed runs independently via RcppParallel (no shared state); `run_workers()`
in `src/parallel_utils.h` runs a Worker on the thread pool for `n_cores > 1`
and serially otherwise, so there are no separate `_single` implementations.
Multiple seeds produce a matrix of curves used for percentile-based confidence
intervals. The Gaussian order generator is serial because it draws from R's
RNG.

## Naming Conventions

- **Functions**: camelCase (`spaccHill`, `spatialRarefaction`). Accumulation-curve front doors carry the `spacc` prefix (`spaccHill`, `spaccBeta`, `spaccWavefront`, ...); distance/area relationship functions keep their established literature names (`distanceDecay`, `betaDecay`, `zetaDiversity`, `dar`).
- **Arguments**: snake_case (`n_seeds`, `n_cores`)
- **Classes**: lowercase with underscore (`spacc`, `spacc_fit`, `spacc_decay`)
- **C++ functions**: snake_case with `cpp_` prefix (`cpp_knn_parallel`)
