# Individual-Based Rarefaction

Compute individual-based rarefaction curves for Hill numbers at any
order q. This complements the sample-based accumulation in
[`spacc()`](https://gillescolling.com/spacc/reference/spacc.md).

## Usage

``` r
rarefy(x, n_individuals = NULL, q = 0, n_boot = 100L)
```

## Arguments

- x:

  A site-by-species matrix with abundance data (not presence/absence).

- n_individuals:

  Integer vector. Sample sizes to compute expected diversity for.
  Default `NULL` computes for all levels from 1 to total.

- q:

  Numeric. Order of Hill number; any value `>= 0`. Default 0 (species
  richness). q=1 gives rarefied Shannon diversity, q=2 gives rarefied
  Simpson diversity. q=0, 1, and 2 use exact estimators; other orders
  report the Hill number of order `q` of the sampled composition.

- n_boot:

  Integer. Number of bootstrap replicates for CI. Default 100.

## Value

An object of class `spacc_rare` containing:

- n:

  Sample sizes

- expected:

  Expected diversity (Hill number of order q)

- sd:

  Standard deviation

- lower, upper:

  95 percent confidence bounds

- q:

  Order of diversity used

## Details

For q=0 (species richness): uses the Hurlbert (1971) formula.

For q=1 (Shannon diversity): rarefied Shannon entropy is computed and
converted to effective number of species via exponentiation.

For q=2 (Simpson diversity): rarefied Simpson concentration is computed
and converted to effective number of species via inversion.

## References

Hurlbert, S.H. (1971). The nonconcept of species diversity: a critique
and alternative parameters. Ecology, 52, 577-586.

Chao, A., Gotelli, N.J., Hsieh, T.C., et al. (2014). Rarefaction and
extrapolation with Hill numbers: a framework for sampling and estimation
in species diversity studies. Ecological Monographs, 84, 45-67.

## Examples

``` r
# \donttest{
abundance_matrix <- matrix(rpois(50 * 30, 2), nrow = 50)
rare <- rarefy(abundance_matrix)
plot(rare)

# Shannon rarefaction
rare_q1 <- rarefy(abundance_matrix, q = 1)
plot(rare_q1)
# }
```
