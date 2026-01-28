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
2.](diversity_files/figure-html/plot-hill-1.png)

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
turnover/nestedness.](diversity_files/figure-html/plot-beta-1.png)

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
accumulation.](diversity_files/figure-html/plot-beta-func-1.png)

Functional beta diversity accumulation.

### Phylogenetic Beta Diversity

[`spaccBetaPhylo()`](https://gillescolling.com/spacc/reference/spaccBetaPhylo.md)
uses phylogenetic distances to weight beta diversity (Chao et al. 2023):

``` r

library(ape)
tree <- rcoal(n_species, tip.label = colnames(species))

beta_phylo <- spaccBetaPhylo(pa, coords, tree, n_seeds = 20, progress = FALSE)
```

``` r

plot(beta_phylo)
```

![Phylogenetic beta diversity
accumulation.](diversity_files/figure-html/plot-beta-phylo-1.png)

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
accumulation.](diversity_files/figure-html/plot-phylo-1.png)

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
accumulation.](diversity_files/figure-html/plot-func-1.png)

Functional diversity accumulation.

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
rarefaction.](diversity_files/figure-html/plot-coverage-1.png)

Coverage-based spatial rarefaction.

Interpolate richness at specific coverage targets:

``` r

interp <- interpolateCoverage(cov, target = c(0.90, 0.95))
summary(interp)
#>       C90             C95       
#>  Min.   :12.57   Min.   :13.52  
#>  1st Qu.:17.65   1st Qu.:19.74  
#>  Median :19.40   Median :21.85  
#>  Mean   :19.57   Mean   :21.51  
#>  3rd Qu.:22.12   3rd Qu.:24.20  
#>  Max.   :25.00   Max.   :27.00
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
