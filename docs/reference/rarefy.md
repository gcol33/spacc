# Individual-Based Rarefaction

Compute classic individual-based rarefaction curves. This complements
the sample-based accumulation in
[`spacc()`](https://gillescolling.com/spacc/reference/spacc.md).

## Usage

``` r
rarefy(x, n_individuals = NULL, n_boot = 100L)
```

## Arguments

- x:

  A site-by-species matrix with abundance data (not presence/absence).

- n_individuals:

  Integer vector. Sample sizes to compute expected richness for. Default
  `NULL` computes for all levels from 1 to total.

- n_boot:

  Integer. Number of bootstrap replicates for CI. Default 100.

## Value

An object of class `spacc_rare` containing:

- n:

  Sample sizes

- expected:

  Expected species richness

- sd:

  Standard deviation

- lower, upper:

  95% confidence bounds

## Examples

``` r
if (FALSE) { # \dontrun{
# With abundance data
rare <- rarefy(abundance_matrix)
plot(rare)
} # }
```
