# spacc v0.2.0 Implementation Plan

Four major features to make spacc paper-worthy for Methods in Ecology and Evolution.

## Overview

| Feature | Priority | Effort | Paper Impact |
|---------|----------|--------|--------------|
| Hill Numbers | 1 | Medium | High |
| Spatial Beta Diversity | 2 | Medium | High |
| Coverage-Based Rarefaction | 3 | Low | Medium |
| Phylogenetic/Functional | 4 | High | High |

---

## Feature 1: Hill Numbers for Spatial SACs

### Background

Hill numbers unify diversity measures:
- q = 0: Species richness (current spacc)
- q = 1: Exponential of Shannon entropy (effective number of common species)
- q = 2: Inverse Simpson (effective number of dominant species)

Formula: `qD = (Σ p_i^q)^(1/(1-q))` where p_i = relative abundance

### API Design

```r
# Extend existing spacc() with q parameter
spacc(x, coords, method = "knn", q = 0, ...)

# Or create dedicated function
hillspacc(x, coords, method = "knn", q = c(0, 1, 2), ...)

# Return object includes curves for each q value
```

**Recommendation**: Add `q` parameter to `spacc()` for simplicity. Default `q = 0` maintains backward compatibility.

### Implementation

#### R Changes (R/spacc.R)

```r
spacc <- function(x,
                  coords,
                  n_seeds = 50L,
                  method = c("knn", "kncn", "random", "radius", "gaussian", "cone", "collector"),
                  distance = c("euclidean", "haversine"),
                  q = 0,
                  ...) {

  # Validate q
  stopifnot("q must be >= 0" = all(q >= 0))

  # If abundance data and q > 0, need to track abundances not just presence
  use_abundance <- any(q > 0)

  # Call appropriate C++ function
  if (use_abundance) {
    curves <- cpp_knn_hill_parallel(x, dist_mat, n_seeds, q, n_cores, progress)
  } else {
    curves <- cpp_knn_parallel(species_pa, dist_mat, n_seeds, n_cores, progress)
  }
}
```

#### C++ Changes (src/hill.cpp)

```cpp
// Hill number calculation for a set of abundances
double calc_hill(const std::vector<int>& abundances, double q) {
  int N = 0;
  for (int a : abundances) N += a;
  if (N == 0) return 0.0;

  if (q == 0) {
    // Species richness
    int S = 0;
    for (int a : abundances) if (a > 0) S++;
    return (double)S;
  } else if (std::abs(q - 1.0) < 1e-10) {
    // Shannon (limit as q -> 1)
    double H = 0.0;
    for (int a : abundances) {
      if (a > 0) {
        double p = (double)a / N;
        H -= p * std::log(p);
      }
    }
    return std::exp(H);
  } else {
    // General Hill number
    double sum = 0.0;
    for (int a : abundances) {
      if (a > 0) {
        double p = (double)a / N;
        sum += std::pow(p, q);
      }
    }
    return std::pow(sum, 1.0 / (1.0 - q));
  }
}

// Single kNN accumulation with Hill numbers
NumericMatrix cpp_knn_hill_single(IntegerMatrix species_mat,
                                   NumericMatrix dist_mat,
                                   int seed,
                                   NumericVector q_values) {
  int n_sites = species_mat.nrow();
  int n_species = species_mat.ncol();
  int n_q = q_values.size();

  NumericMatrix curves(n_q, n_sites);
  std::vector<bool> visited(n_sites, false);
  std::vector<int> cumulative_abundance(n_species, 0);

  int current = seed;
  visited[current] = true;

  // Add first site abundances
  for (int sp = 0; sp < n_species; sp++) {
    cumulative_abundance[sp] += species_mat(current, sp);
  }

  // Calculate Hill numbers for each q
  for (int qi = 0; qi < n_q; qi++) {
    curves(qi, 0) = calc_hill(cumulative_abundance, q_values[qi]);
  }

  for (int step = 1; step < n_sites; step++) {
    // Find nearest unvisited
    double min_dist = R_PosInf;
    int next = -1;
    for (int j = 0; j < n_sites; j++) {
      if (!visited[j] && dist_mat(current, j) < min_dist) {
        min_dist = dist_mat(current, j);
        next = j;
      }
    }

    current = next;
    visited[current] = true;

    // Accumulate abundances
    for (int sp = 0; sp < n_species; sp++) {
      cumulative_abundance[sp] += species_mat(current, sp);
    }

    // Calculate Hill numbers
    for (int qi = 0; qi < n_q; qi++) {
      curves(qi, step) = calc_hill(cumulative_abundance, q_values[qi]);
    }
  }

  return curves;
}
```

#### New S3 Class: spacc_hill

```r
# Return structure
structure(
  list(
    curves = list(q0 = matrix, q1 = matrix, q2 = matrix),
    q = c(0, 1, 2),
    coords = coord_data,
    n_seeds = n_seeds,
    n_sites = n_sites,
    n_species = n_species_total,
    method = method,
    call = match.call()
  ),
  class = c("spacc_hill", "spacc")
)
```

### Plot Method

```r
plot.spacc_hill <- function(x, q = NULL, ci = TRUE, ...) {
  # If q not specified, plot all
  q_vals <- if (is.null(q)) x$q else q

  # Create faceted plot with one panel per q value
  # Or overlay with different colors
}
```

### Tests

```r
test_that("Hill numbers are correct", {
  # Known community: 10 species, abundances 1:10
  # q=0 should give 10
  # q=1 should give exp(Shannon)
  # q=2 should give 1/Simpson

  x <- matrix(1:10, nrow = 1)
  coords <- data.frame(x = 0, y = 0)

  result <- spacc(x, coords, method = "collector", q = c(0, 1, 2))

  expect_equal(result$curves$q0[1, 1], 10)
  # ... verify against vegan::diversity()
})
```

---

## Feature 2: Spatial Beta Diversity Accumulation

### Background

Beta diversity measures compositional change. As we accumulate sites spatially:
- **Turnover**: Species replacement (A has species X, B has species Y)
- **Nestedness**: Species loss (B is a subset of A)

Baselga (2010) framework: β_total = β_turnover + β_nestedness

### API Design

```r
# New function
spaccBeta(x, coords, method = "knn", n_seeds = 50,
          partition = TRUE, index = "sorensen", ...)

# Returns beta diversity vs cumulative sites/distance
```

### Implementation

#### Core Algorithm

At each step of spatial accumulation:
1. Track species composition of accumulated set
2. Compare to previous step (pairwise) OR to all accumulated (cumulative)
3. Partition into turnover and nestedness components

```r
# Sorensen-based dissimilarity
beta_sor <- function(a, b, c) {
  # a = shared species, b = unique to site 1, c = unique to site 2
  (b + c) / (2*a + b + c)
}

# Turnover component (Simpson)
beta_sim <- function(a, b, c) {
  min(b, c) / (a + min(b, c))
}

# Nestedness component
beta_nes <- function(a, b, c) {
  beta_sor(a, b, c) - beta_sim(a, b, c)
}
```

#### C++ Implementation (src/beta.cpp)

```cpp
struct BetaResult {
  double total;
  double turnover;
  double nestedness;
};

BetaResult calc_beta_sorensen(const std::set<int>& set1,
                               const std::set<int>& set2) {
  int a = 0;  // shared
  int b = 0;  // unique to set1
  int c = 0;  // unique to set2

  for (int sp : set1) {
    if (set2.count(sp)) a++;
    else b++;
  }
  for (int sp : set2) {
    if (!set1.count(sp)) c++;
  }

  BetaResult result;
  if (a + b + c == 0) {
    result.total = 0;
    result.turnover = 0;
    result.nestedness = 0;
  } else {
    result.total = (double)(b + c) / (2.0*a + b + c);
    result.turnover = (double)std::min(b, c) / (a + std::min(b, c));
    result.nestedness = result.total - result.turnover;
  }
  return result;
}

// Track cumulative beta as sites are added
List cpp_beta_accumulation(IntegerMatrix species_pa,
                           NumericMatrix dist_mat,
                           int seed,
                           std::string method) {
  int n_sites = species_pa.nrow();
  int n_species = species_pa.ncol();

  NumericVector beta_total(n_sites - 1);
  NumericVector beta_turn(n_sites - 1);
  NumericVector beta_nest(n_sites - 1);
  NumericVector distances(n_sites - 1);

  std::vector<bool> visited(n_sites, false);
  std::set<int> prev_species;
  std::set<int> curr_species;

  int current = seed;
  visited[current] = true;

  // Initialize with first site
  for (int sp = 0; sp < n_species; sp++) {
    if (species_pa(current, sp) > 0) {
      prev_species.insert(sp);
    }
  }

  double cum_dist = 0.0;

  for (int step = 0; step < n_sites - 1; step++) {
    // Find next site (kNN logic)
    double min_dist = R_PosInf;
    int next = -1;
    for (int j = 0; j < n_sites; j++) {
      if (!visited[j] && dist_mat(current, j) < min_dist) {
        min_dist = dist_mat(current, j);
        next = j;
      }
    }

    cum_dist += min_dist;
    distances[step] = cum_dist;

    // Get species at new site
    curr_species.clear();
    for (int sp = 0; sp < n_species; sp++) {
      if (species_pa(next, sp) > 0) {
        curr_species.insert(sp);
      }
    }

    // Calculate pairwise beta between accumulated and new site
    BetaResult beta = calc_beta_sorensen(prev_species, curr_species);
    beta_total[step] = beta.total;
    beta_turn[step] = beta.turnover;
    beta_nest[step] = beta.nestedness;

    // Update accumulated species set
    for (int sp : curr_species) {
      prev_species.insert(sp);
    }

    current = next;
    visited[current] = true;
  }

  return List::create(
    Named("beta_total") = beta_total,
    Named("beta_turnover") = beta_turn,
    Named("beta_nestedness") = beta_nest,
    Named("distance") = distances
  );
}
```

#### R Wrapper

```r
spaccBeta <- function(x, coords, method = "knn", n_seeds = 50,
                      partition = TRUE, index = c("sorensen", "jaccard"),
                      distance = "euclidean", progress = TRUE, seed = NULL) {

  index <- match.arg(index)

  # ... validation ...

  results <- cpp_beta_parallel(species_pa, dist_mat, n_seeds, index, n_cores)

  structure(
    list(
      beta_total = results$total,      # matrix: n_seeds x (n_sites-1)
      beta_turnover = results$turnover,
      beta_nestedness = results$nestedness,
      distance = results$distance,
      n_seeds = n_seeds,
      method = method,
      index = index,
      call = match.call()
    ),
    class = "spacc_beta"
  )
}
```

#### Plot Method

```r
plot.spacc_beta <- function(x, partition = TRUE, xaxis = c("sites", "distance"), ...) {
  xaxis <- match.arg(xaxis)

  # Summarize across seeds
  df <- data.frame(
    x = if (xaxis == "sites") 2:x$n_sites else colMeans(x$distance),
    total = colMeans(x$beta_total),
    turnover = colMeans(x$beta_turnover),
    nestedness = colMeans(x$beta_nestedness)
  )

  if (partition) {
    # Stacked area chart: turnover + nestedness = total
    # Or line plot with three lines
  } else {
    # Just total beta
  }
}
```

---

## Feature 3: Coverage-Based Spatial Rarefaction

### Background

Chao & Jost (2012): Sample completeness (coverage) is better than sample size for standardization.

Coverage = proportion of total community abundance represented by observed species.
`C = 1 - (f1/n)` where f1 = singletons, n = total individuals (simplified)

### API Design

```r
# Add standardize parameter to spacc()
spacc(x, coords, method = "knn", standardize = c("sites", "coverage"), ...)

# Or separate function
spaccCoverage(x, coords, method = "knn", target_coverage = 0.95, ...)
```

### Implementation

#### Coverage Calculation (Chao & Jost 2012)

```cpp
// Good-Turing coverage estimator
double calc_coverage(const std::vector<int>& abundances) {
  int n = 0;  // total individuals
  int f1 = 0; // singletons
  int f2 = 0; // doubletons

  for (int a : abundances) {
    n += a;
    if (a == 1) f1++;
    if (a == 2) f2++;
  }

  if (n == 0) return 0.0;

  // Chao & Jost estimator
  double C;
  if (f2 > 0) {
    C = 1.0 - (double)f1 / n * ((n - 1.0) * f1 / ((n - 1.0) * f1 + 2.0 * f2));
  } else {
    C = 1.0 - (double)f1 / n * ((n - 1.0) / n);
  }

  return std::max(0.0, std::min(1.0, C));
}
```

#### Spatial Accumulation with Coverage

```cpp
List cpp_knn_coverage(IntegerMatrix species_mat,
                      NumericMatrix dist_mat,
                      int seed) {
  int n_sites = species_mat.nrow();
  int n_species = species_mat.ncol();

  IntegerVector richness(n_sites);
  NumericVector coverage(n_sites);

  std::vector<bool> visited(n_sites, false);
  std::vector<int> cumulative(n_species, 0);

  int current = seed;
  visited[current] = true;

  for (int sp = 0; sp < n_species; sp++) {
    cumulative[sp] += species_mat(current, sp);
  }

  richness[0] = count_species(cumulative);
  coverage[0] = calc_coverage(cumulative);

  for (int step = 1; step < n_sites; step++) {
    // kNN logic to find next site
    // ...

    // Update cumulative abundances
    for (int sp = 0; sp < n_species; sp++) {
      cumulative[sp] += species_mat(current, sp);
    }

    richness[step] = count_species(cumulative);
    coverage[step] = calc_coverage(cumulative);
  }

  return List::create(
    Named("richness") = richness,
    Named("coverage") = coverage
  );
}
```

#### Interpolation to Target Coverage

```r
# Interpolate richness at target coverage levels
interpolate_coverage <- function(richness, coverage, target = c(0.9, 0.95, 0.99)) {
  sapply(target, function(t) {
    if (max(coverage) < t) return(NA)
    approx(coverage, richness, xout = t)$y
  })
}
```

---

## Feature 4: Phylogenetic & Functional Diversity

### Background

- **Faith's PD**: Sum of branch lengths connecting species in a sample
- **Functional Diversity (FD)**: Trait-based diversity (e.g., functional richness, dispersion)

### API Design

```r
# Phylogenetic
spaccPhylo(x, coords, tree, method = "knn", metric = c("PD", "MPD", "MNTD"), ...)

# Functional
spaccFunc(x, coords, traits, method = "knn", metric = c("FRic", "FDiv", "FDis"), ...)
```

### Dependencies

```r
Suggests:
    ape,        # phylogenetic trees
    picante,    # PD calculations (optional, can implement ourselves)
    FD,         # functional diversity
    geometry    # convex hulls for FRic
```

### Phylogenetic Implementation

#### Faith's PD Calculation

```r
# Using ape for tree manipulation
calc_pd <- function(tree, tips) {
  if (length(tips) == 0) return(0)
  if (length(tips) == 1) {
    # PD of single tip = edge length to root
    node <- which(tree$tip.label == tips)
    return(sum(ape::node.depth.edgelength(tree)[node]))
  }

  # Prune tree to tips
  subtree <- ape::keep.tip(tree, tips)
  sum(subtree$edge.length)
}
```

#### C++ Implementation for Performance

```cpp
// Pre-compute pairwise phylogenetic distances
// Then PD = sum of unique branch lengths

// Alternative: cophenetic distance matrix approach
// MPD = mean pairwise phylogenetic distance
double calc_mpd(const NumericMatrix& phylo_dist,
                const std::set<int>& species) {
  if (species.size() < 2) return 0.0;

  double sum = 0.0;
  int count = 0;

  std::vector<int> spp(species.begin(), species.end());
  for (size_t i = 0; i < spp.size(); i++) {
    for (size_t j = i + 1; j < spp.size(); j++) {
      sum += phylo_dist(spp[i], spp[j]);
      count++;
    }
  }

  return sum / count;
}
```

#### Spatial PD Accumulation

```cpp
NumericVector cpp_pd_accumulation(IntegerMatrix species_pa,
                                   NumericMatrix dist_mat,
                                   NumericMatrix phylo_dist,
                                   NumericVector edge_lengths,
                                   int seed) {
  // Track which edges are "activated" as species accumulate
  // PD = sum of activated edge lengths

  // Implementation depends on tree structure
  // May need to pass tree as edge list
}
```

### Functional Diversity Implementation

#### Functional Richness (Convex Hull Volume)

```r
calc_fric <- function(traits, species_present) {
  trait_subset <- traits[species_present, , drop = FALSE]

  if (nrow(trait_subset) <= ncol(trait_subset)) {
    return(NA)  # Need more species than traits for convex hull
  }

  # Convex hull volume
  hull <- geometry::convhulln(trait_subset, options = "FA")
  hull$vol
}
```

#### Functional Dispersion

```r
calc_fdis <- function(traits, species_present, abundances = NULL) {
  trait_subset <- traits[species_present, , drop = FALSE]

  if (is.null(abundances)) {
    abundances <- rep(1, nrow(trait_subset))
  }

  # Weighted centroid
  centroid <- colSums(trait_subset * abundances) / sum(abundances)

  # Mean distance to centroid
  dists <- sqrt(rowSums((trait_subset - centroid)^2))
  sum(dists * abundances) / sum(abundances)
}
```

---

## File Structure After Implementation

```
spacc/
├── R/
│   ├── spacc.R           # Updated: add q parameter
│   ├── spaccBeta.R       # NEW: beta diversity accumulation
│   ├── spaccCoverage.R   # NEW: coverage-based methods
│   ├── spaccPhylo.R      # NEW: phylogenetic diversity
│   ├── spaccFunc.R       # NEW: functional diversity
│   ├── hill.R            # NEW: Hill number utilities
│   └── ...
├── src/
│   ├── hill.cpp          # NEW: Hill number calculations
│   ├── beta.cpp          # NEW: beta diversity calculations
│   ├── coverage.cpp      # NEW: coverage estimation
│   ├── phylo.cpp         # NEW: PD calculations (optional)
│   └── ...
└── tests/testthat/
    ├── test-hill.R       # NEW
    ├── test-beta.R       # NEW
    ├── test-coverage.R   # NEW
    ├── test-phylo.R      # NEW
    └── ...
```

---

## Implementation Order

### Phase 1: Hill Numbers (Week 1-2)
1. Add `calc_hill()` C++ function
2. Modify `cpp_knn_single` to track abundances
3. Create `cpp_knn_hill_parallel`
4. Update R wrapper with `q` parameter
5. Add plot method for multiple q values
6. Write tests comparing to vegan::diversity()
7. Document

### Phase 2: Spatial Beta Diversity (Week 2-3)
1. Implement `calc_beta_sorensen()` in C++
2. Create `cpp_beta_accumulation()`
3. Add parallel version
4. Create `spaccBeta()` R function
5. Implement S3 methods (print, summary, plot)
6. Write tests comparing to betapart
7. Document

### Phase 3: Coverage-Based Rarefaction (Week 3)
1. Implement `calc_coverage()` in C++
2. Add coverage tracking to existing methods
3. Create interpolation function
4. Update spacc() or create spaccCoverage()
5. Write tests comparing to iNEXT
6. Document

### Phase 4: Phylogenetic/Functional (Week 4-5)
1. Design tree input format (ape::phylo)
2. Implement PD accumulation (start with R, optimize to C++ if needed)
3. Implement FD accumulation
4. Create spaccPhylo() and spaccFunc()
5. Write tests comparing to picante/FD
6. Document

---

## Testing Strategy

### Unit Tests
- Each C++ function has corresponding R test
- Compare against reference implementations (vegan, iNEXT, betapart, picante)
- Edge cases: empty data, single species, single site

### Integration Tests
- Full workflow tests with simulated data
- Verify plot methods don't error
- Check parallel vs sequential give same results

### Performance Tests
- Benchmark against mobr, vegan
- Document speedups in vignette

---

## Documentation Plan

### Function Documentation
- All exported functions have roxygen2 docs
- Examples for each function
- Cross-references with @seealso

### Vignette: "Comprehensive Spatial Diversity Analysis"

1. **Introduction**: Why spatial matters
2. **Basic Usage**: spacc() with richness
3. **Hill Numbers**: Comparing q=0,1,2
4. **Beta Diversity**: Turnover vs nestedness
5. **Coverage-Based**: When to use
6. **Phylogenetic/Functional**: Advanced applications
7. **Performance**: Benchmarks vs alternatives
8. **Case Study**: Real ecological dataset

---

## Paper Outline (MEE Application)

**Title**: "spacc: Fast spatial diversity accumulation with Hill numbers, beta partitioning, and phylogenetic extensions"

1. **Introduction** (500 words)
   - SACs and their limitations
   - Spatial structure matters
   - Gap: no fast, comprehensive tool

2. **Methods** (1500 words)
   - Hill number spatial accumulation
   - Beta diversity partitioning
   - Coverage standardization
   - Phylogenetic/functional extensions
   - C++ implementation

3. **Results** (1000 words)
   - Simulation study: bias of random vs spatial
   - Benchmark: performance vs mobr/vegan
   - Case study: real dataset

4. **Discussion** (500 words)
   - When to use which method
   - Limitations
   - Future directions

5. **Figures**
   - Fig 1: Conceptual diagram of spatial vs random SAC
   - Fig 2: Hill numbers reveal scale-dependent patterns
   - Fig 3: Beta partitioning example
   - Fig 4: Performance benchmarks

---

## Version Milestones

### v0.2.0 - Hill Numbers
- `q` parameter in spacc()
- Basic Hill number support

### v0.3.0 - Beta Diversity
- `spaccBeta()` function
- Turnover/nestedness partitioning

### v0.4.0 - Coverage
- Coverage-based standardization
- Interpolation to target coverage

### v0.5.0 - Phylogenetic/Functional
- `spaccPhylo()` and `spaccFunc()`
- Requires ape, geometry in Suggests

### v1.0.0 - Paper Release
- Full documentation
- Vignette
- CRAN submission
- MEE submission
