# Bootstrap Richness Estimator

Bootstrap estimator of total species richness (Smith & van Belle 1984).
Uses species detection probabilities to estimate undetected species.

## Usage

``` r
bootstrap_richness(x, n_boot = 200L)
```

## Arguments

- x:

  A site-by-species matrix (presence/absence or abundance).
  Automatically binarized.

- n_boot:

  Integer. Number of bootstrap replicates for SE. Default 200.

## Value

An object of class `spacc_estimate`.

## Details

The bootstrap estimator is: \$\$S\_{boot} = S\_{obs} +
\sum\_{i=1}^{S\_{obs}} (1 - p_i)^n\$\$ where \\p_i\\ is the proportion
of sites where species \\i\\ occurs.

## References

Smith, E.P. & van Belle, G. (1984). Nonparametric estimation of species
richness. Biometrics, 40, 119-129.

## See also

[`jackknife()`](https://gillescolling.com/spacc/reference/jackknife.md),
[`chao2()`](https://gillescolling.com/spacc/reference/chao2.md)

## Examples

``` r
species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)
bootstrap_richness(species, n_boot = 100)
```
