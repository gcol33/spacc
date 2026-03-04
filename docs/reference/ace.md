# ACE Richness Estimator

Abundance-based Coverage Estimator (ACE) of total species richness (Chao
& Lee 1992). Separates species into rare and abundant groups based on an
abundance threshold.

## Usage

``` r
ace(x, threshold = 10L)
```

## Arguments

- x:

  A site-by-species matrix (abundance data).

- threshold:

  Integer. Abundance threshold separating rare from abundant species.
  Default 10.

## Value

An object of class `spacc_estimate`.

## Details

ACE partitions species into rare (abundance \\\le\\ threshold) and
abundant (abundance \> threshold) groups. The estimate is: \$\$S\_{ACE}
= S\_{abund} + \frac{S\_{rare}}{C\_{ACE}} + \frac{f_1}{C\_{ACE}}
\gamma^2\$\$ where \\C\_{ACE}\\ is the sample coverage of rare species
and \\\gamma^2\\ is the estimated coefficient of variation.

## References

Chao, A. & Lee, S.M. (1992). Estimating the number of classes via sample
coverage. Journal of the American Statistical Association, 87, 210-217.

## See also

[`chao1()`](https://gillescolling.com/spacc/reference/chao1.md),
[`chao2()`](https://gillescolling.com/spacc/reference/chao2.md)

## Examples

``` r
species <- matrix(rpois(50 * 30, 2), nrow = 50)
ace(species)
```
