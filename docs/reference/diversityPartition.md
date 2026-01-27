# Alpha-Beta-Gamma Diversity Partitioning

Decompose regional (gamma) diversity into local (alpha) and turnover
(beta) components using multiplicative partitioning of Hill numbers.

## Usage

``` r
diversityPartition(x, q = c(0, 1, 2), weights = "equal", coords = NULL)
```

## Arguments

- x:

  A site-by-species matrix (abundance data).

- q:

  Numeric vector. Orders of diversity. Default `c(0, 1, 2)`.

- weights:

  Character or numeric. Site weights for alpha calculation:

  - `"equal"`: Equal weights (default)

  - `"proportional"`: Weights proportional to site abundance

  - Numeric vector of custom weights

- coords:

  Optional data.frame with columns `x` and `y` for spatial mapping. When
  provided, enables
  [`as_sf()`](https://gillescolling.com/spacc/reference/as_sf.md) and
  `plot(type = "map")`.

## Value

An object of class `spacc_partition` containing:

- alpha:

  Mean alpha diversity (effective number of species per site)

- beta:

  Beta diversity (effective number of communities)

- gamma:

  Gamma diversity (regional species pool)

- q:

  Orders of diversity

- n_sites:

  Number of sites

- n_species:

  Total species count

## Details

This function implements Jost (2007) multiplicative partitioning:

\$\$\gamma = \alpha \times \beta\$\$

Where:

- **Alpha**: Mean effective number of species per site

- **Beta**: Effective number of distinct communities (1 = all identical,
  n_sites = all completely different)

- **Gamma**: Total effective number of species in the region

Beta diversity is interpreted as the effective number of communities:

- Beta = 1: All sites have identical composition

- Beta = n_sites: Sites share no species

## References

Jost, L. (2007). Partitioning diversity into independent alpha and beta
components. Ecology, 88, 2427-2439.

Chao, A., Chiu, C.H. & Jost, L. (2014). Unifying species diversity,
phylogenetic diversity, functional diversity, and related similarity and
differentiation measures through Hill numbers. Annual Review of Ecology,
Evolution, and Systematics, 45, 297-324.

## See also

[`alphaDiversity()`](https://gillescolling.com/spacc/reference/alphaDiversity.md),
[`gammaDiversity()`](https://gillescolling.com/spacc/reference/gammaDiversity.md),
[`spaccBeta()`](https://gillescolling.com/spacc/reference/spaccBeta.md)
for spatial beta diversity accumulation

## Examples

``` r
# Simulated community data
species <- matrix(rpois(50 * 30, 2), nrow = 50)

# Partition diversity
part <- diversityPartition(species, q = c(0, 1, 2))
print(part)

# Beta near 1 = homogeneous region
# Beta near n_sites = heterogeneous region
```
