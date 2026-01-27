# Alpha Diversity (Per-Site)

Compute Hill numbers for each site individually.

## Usage

``` r
alphaDiversity(x, q = c(0, 1, 2), coords = NULL)
```

## Arguments

- x:

  A site-by-species matrix (abundance data).

- q:

  Numeric vector. Orders of diversity. Default `c(0, 1, 2)`.

- coords:

  Optional data.frame with columns `x` and `y` for spatial mapping. When
  provided, returns a `spacc_alpha` object with
  [`as_sf()`](https://gillescolling.com/spacc/reference/as_sf.md) and
  `plot(type = "map")` support.

## Value

If `coords` is `NULL`, a data.frame with columns for each q value. If
`coords` is provided, a `spacc_alpha` object.

## Details

Alpha diversity represents local (within-site) diversity. For Hill
numbers:

- q=0: Species richness

- q=1: Exponential of Shannon entropy

- q=2: Inverse Simpson concentration

## References

Jost, L. (2007). Partitioning diversity into independent alpha and beta
components. Ecology, 88, 2427-2439.

## See also

[`gammaDiversity()`](https://gillescolling.com/spacc/reference/gammaDiversity.md)
for regional diversity,
[`diversityPartition()`](https://gillescolling.com/spacc/reference/diversityPartition.md)
for full alpha-beta-gamma decomposition

## Examples

``` r
species <- matrix(rpois(50 * 30, 2), nrow = 50)
alpha <- alphaDiversity(species, q = c(0, 1, 2))
head(alpha)

# Mean alpha diversity
colMeans(alpha)
```
