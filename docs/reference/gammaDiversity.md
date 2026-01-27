# Gamma Diversity (Regional)

Compute Hill numbers for the pooled community across all sites.

## Usage

``` r
gammaDiversity(x, q = c(0, 1, 2))
```

## Arguments

- x:

  A site-by-species matrix (abundance data).

- q:

  Numeric vector. Orders of diversity. Default `c(0, 1, 2)`.

## Value

A named numeric vector with gamma diversity for each q.

## Details

Gamma diversity represents regional (total) diversity across all sites.
It is computed by pooling abundances across all sites and calculating
Hill numbers on the combined community.

## References

Jost, L. (2007). Partitioning diversity into independent alpha and beta
components. Ecology, 88, 2427-2439.

## See also

[`alphaDiversity()`](https://gillescolling.com/spacc/reference/alphaDiversity.md)
for local diversity,
[`diversityPartition()`](https://gillescolling.com/spacc/reference/diversityPartition.md)
for full alpha-beta-gamma decomposition

## Examples

``` r
species <- matrix(rpois(50 * 30, 2), nrow = 50)
gammaDiversity(species, q = c(0, 1, 2))
```
