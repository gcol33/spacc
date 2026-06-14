# Diversity Accumulation

## Overview

Beyond simple species richness (q=0), spacc supports accumulation curves
for the full family of Hill numbers, beta diversity partitioning,
phylogenetic diversity, and functional diversity. This vignette
demonstrates each approach.

## Simulating Data

``` r

library(spacc)

set.seed(123)
n_sites <- 60
n_species <- 30

coords <- data.frame(
  x = runif(n_sites, 0, 100),
  y = runif(n_sites, 0, 100)
)

# Abundance matrix with spatial clustering
species <- matrix(0L, n_sites, n_species)
for (sp in seq_len(n_species)) {
  cx <- runif(1, 10, 90)
  cy <- runif(1, 10, 90)
  lambda <- 5 * exp(-0.001 * ((coords$x - cx)^2 + (coords$y - cy)^2))
  species[, sp] <- rpois(n_sites, lambda)
}
colnames(species) <- paste0("sp", seq_len(n_species))
```

## Alpha, Beta, and Gamma Diversity

spacc provides functions for the multiplicative decomposition of
diversity into alpha (local), beta (turnover), and gamma (regional)
components (Jost 2007):

``` r

# Alpha diversity: per-site Hill numbers
alpha <- alphaDiversity(species, q = c(0, 1, 2))
colMeans(alpha)
#>       q0       q1       q2 
#> 15.06667 12.39927 10.61170

# Gamma diversity: pooled regional diversity
gamma <- gammaDiversity(species, q = c(0, 1, 2))
gamma
#>       q0       q1       q2 
#> 30.00000 29.23513 28.56626

# Full partition: gamma = alpha * beta
partition <- diversityPartition(species, q = c(0, 1, 2))
partition
#> Alpha-Beta-Gamma Diversity Partitioning
#> 60 sites, 30 species
#> 
#>  q alpha beta gamma
#>  0 15.07 1.99 30.00
#>  1 11.42 2.56 29.24
#>  2  8.81 3.24 28.57
#> 
#> Interpretation:
#>   Alpha = mean effective species per site
#>   Beta  = effective number of communities (1 to n_sites)
#>   Gamma = regional effective species (gamma = alpha x beta)
```

## Hill Number Accumulation

Hill numbers unify richness (q=0), Shannon diversity (exponential, q=1),
and Simpson diversity (inverse, q=2) into a single framework (Jost 2007,
Chao et al. 2014):

``` r

hill <- spaccHill(species, coords, q = c(0, 1, 2), n_seeds = 20, progress = FALSE)
```

``` r

plot(hill)
```

![Hill number accumulation for q = 0, 1,
2.](diversity_files/figure-html/plot-hill-1.svg)

Hill number accumulation for q = 0, 1, 2.

Higher-order q values emphasise dominant species, so diversity
accumulates more slowly at q=2 compared to q=0.

## Beta Diversity

[`spaccBeta()`](https://gillescolling.com/spacc/reference/spaccBeta.md)
partitions spatial beta diversity into turnover and nestedness
components (Baselga 2010):

``` r

pa <- (species > 0) * 1L
beta <- spaccBeta(pa, coords, n_seeds = 20, progress = FALSE)
```

``` r

plot(beta)
```

![Spatial beta diversity accumulation with
turnover/nestedness.](diversity_files/figure-html/plot-beta-1.svg)

Spatial beta diversity accumulation with turnover/nestedness.

Turnover reflects species replacement along spatial gradients, while
nestedness captures diversity loss at species-poor sites.

### Functional Beta Diversity

[`spaccBetaFunc()`](https://gillescolling.com/spacc/reference/spaccBetaFunc.md)
weights beta diversity by trait dissimilarity (Baselga 2012):

``` r

# Simulate two continuous traits
traits <- data.frame(
  body_size = rnorm(n_species),
  wing_length = rnorm(n_species)
)
rownames(traits) <- colnames(species)

beta_func <- spaccBetaFunc(pa, coords, traits, n_seeds = 20, progress = FALSE)
```

``` r

plot(beta_func)
```

![Functional beta diversity
accumulation.](diversity_files/figure-html/plot-beta-func-1.svg)

Functional beta diversity accumulation.

### Phylogenetic Beta Diversity

[`spaccBetaPhylo()`](https://gillescolling.com/spacc/reference/spaccBetaPhylo.md)
uses phylogenetic distances to weight beta diversity (Chao et al. 2023):

``` r

library(ape)
#> 
#> Attaching package: 'ape'
#> The following object is masked from 'package:spacc':
#> 
#>     ace
tree <- rcoal(n_species, tip.label = colnames(species))

beta_phylo <- spaccBetaPhylo(pa, coords, tree, n_seeds = 20, progress = FALSE)
```

``` r

plot(beta_phylo)
```

![Phylogenetic beta diversity
accumulation.](diversity_files/figure-html/plot-beta-phylo-1.svg)

Phylogenetic beta diversity accumulation.

## Phylogenetic Diversity Accumulation

[`spaccPhylo()`](https://gillescolling.com/spacc/reference/spaccPhylo.md)
tracks phylogenetic diversity metrics (MPD, MNTD, Faith’s PD) as sites
accumulate spatially:

``` r

phylo_acc <- spaccPhylo(pa, coords, tree,
                        metric = c("mpd", "mntd", "pd"),
                        n_seeds = 20, progress = FALSE)
```

``` r

plot(phylo_acc)
```

![Phylogenetic diversity
accumulation.](diversity_files/figure-html/plot-phylo-1.svg)

Phylogenetic diversity accumulation.

## Functional Diversity Accumulation

[`spaccFunc()`](https://gillescolling.com/spacc/reference/spaccFunc.md)
tracks functional diversity metrics (FDis, FRic) as sites accumulate:

``` r

func_acc <- spaccFunc(species, coords, traits,
                      metric = c("fdis"),
                      n_seeds = 20, progress = FALSE)
```

``` r

plot(func_acc)
```

![Functional diversity
accumulation.](diversity_files/figure-html/plot-func-1.svg)

Functional diversity accumulation.

## Rao’s Quadratic Entropy

Rao’s quadratic entropy combines relative abundance with the pairwise
distance between species, giving an abundance-weighted measure of
phylogenetic or functional dispersion (Rao 1982). Both
[`spaccPhylo()`](https://gillescolling.com/spacc/reference/spaccPhylo.md)
and
[`spaccFunc()`](https://gillescolling.com/spacc/reference/spaccFunc.md)
expose it through `metric = "rao"`.

``` r

phylo_rao <- spaccPhylo(species, coords, tree,
                        metric = "rao",
                        n_seeds = 20, progress = FALSE)
```

``` r

plot(phylo_rao)
```

![Phylogenetic Rao's Q
accumulation.](diversity_files/figure-html/plot-rao-phylo-1.svg)

Phylogenetic Rao’s Q accumulation.

For functional Rao the species distance is the Euclidean distance in
trait space:

``` r

func_rao <- spaccFunc(species, coords, traits,
                      metric = "rao",
                      n_seeds = 20, progress = FALSE)
```

With presence/absence data Rao’s Q reduces to an equal-weight form; with
abundances it weights each species pair by the product of relative
abundances. The phylogenetic and functional backends carry abundances in
double precision, so continuous values such as percent cover or biomass
are used directly:

``` r

cover <- species / max(species)   # fractional values in [0, 1]
func_rao_cover <- spaccFunc(cover, coords, traits,
                            metric = "rao", n_seeds = 10, progress = FALSE)
tail(as.data.frame(func_rao_cover), 3)
#>    sites metric     mean    lower    upper           sd
#> 58    58    rao 1.762876 1.761032 1.764854 1.887203e-03
#> 59    59    rao 1.762211 1.759666 1.763262 1.292454e-03
#> 60    60    rao 1.761422 1.761422 1.761422 5.074204e-16
```

## Custom Diversity Metrics

[`spaccDiversity()`](https://gillescolling.com/spacc/reference/spaccDiversity.md)
accumulates any user-supplied index along a spatial ordering. At each
step the cumulative community is passed to a function that returns a
single number, so indices spacc does not implement directly can still be
tracked along the accumulation curve.

``` r

# Shannon entropy of the cumulative community
shannon <- function(comm) {
  p <- comm[comm > 0] / sum(comm)
  -sum(p * log(p))
}

div <- spaccDiversity(species, coords, shannon,
                      method = "knn", n_seeds = 20, progress = FALSE)
```

``` r

plot(div, ylab = "Cumulative Shannon entropy")
```

![Custom (Shannon entropy) accumulation
curve.](diversity_files/figure-html/plot-custom-1.svg)

Custom (Shannon entropy) accumulation curve.

The function receives a named vector of cumulative abundances (or 0/1
incidences when `incidence = TRUE`), plus any extra arguments passed
through `...`. The result inherits the `spacc` class, so
[`summary()`](https://rdrr.io/r/base/summary.html),
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), and
[`predict()`](https://rdrr.io/r/stats/predict.html) apply. The
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) method uses a
metric-neutral axis label by default; set `ylab` to name your index.
Available orderings are `"knn"`, `"kncn"`, `"random"`, `"radius"`, and
`"collector"`.

## Coverage-Based Rarefaction

[`spaccCoverage()`](https://gillescolling.com/spacc/reference/spaccCoverage.md)
computes accumulation curves with sample coverage tracking (Chao & Jost
2012), enabling standardization by completeness:

``` r

cov <- spaccCoverage(species, coords, n_seeds = 20, progress = FALSE)
```

``` r

plot(cov)
```

![Coverage-based spatial
rarefaction.](diversity_files/figure-html/plot-coverage-1.svg)

Coverage-based spatial rarefaction.

Interpolate richness at specific coverage targets:

``` r

interp <- interpolateCoverage(cov, target = c(0.90, 0.95))
summary(interp)
#>       C90             C95       
#>  Min.   :10.80   Min.   :15.78  
#>  1st Qu.:15.96   1st Qu.:19.75  
#>  Median :18.81   Median :21.41  
#>  Mean   :18.74   Mean   :21.53  
#>  3rd Qu.:22.12   3rd Qu.:24.14  
#>  Max.   :23.84   Max.   :25.36
```

## References

- Baselga, A. (2010). Partitioning the turnover and nestedness
  components of beta diversity. Global Ecology and Biogeography, 19,
  134-143.
- Baselga, A. (2012). The relationship between species replacement,
  dissimilarity derived from nestedness, and nestedness. Global Ecology
  and Biogeography, 21, 1223-1232.
- Chao, A. & Jost, L. (2012). Coverage-based rarefaction and
  extrapolation. Ecology, 93, 2533-2547.
- Chao, A., Gotelli, N.J., Hsieh, T.C., et al. (2014). Rarefaction and
  extrapolation with Hill numbers. Ecological Monographs, 84, 45-67.
- Faith, D.P. (1992). Conservation evaluation and phylogenetic
  diversity. Biological Conservation, 61, 1-10.
- Jost, L. (2007). Partitioning diversity into independent alpha and
  beta components. Ecology, 88, 2427-2439.
- Rao, C.R. (1982). Diversity and dissimilarity coefficients: a unified
  approach. Theoretical Population Biology, 21, 24-43.
