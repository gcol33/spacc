# Chao2 Richness Estimator

Estimate total species richness from incidence (presence/absence) data
using the Chao2 estimator (Chao 1987).

## Usage

``` r
chao2(x)
```

## Arguments

- x:

  A site-by-species matrix (presence/absence or abundance).
  Automatically binarized.

## Value

An object of class `spacc_estimate`.

## Details

The Chao2 estimator is the incidence-based analogue of Chao1:
\$\$S\_{Chao2} = S\_{obs} + \frac{Q_1^2}{2 Q_2}\$\$ where \\Q_1\\ is the
number of uniques (species found at exactly 1 site) and \\Q_2\\ is the
number of duplicates (species found at exactly 2 sites).

## References

Chao, A. (1987). Estimating the population size for capture-recapture
data with unequal catchability. Biometrics, 43, 783-791.

## See also

[`chao1()`](https://gillescolling.com/spacc/reference/chao1.md) for
abundance-based estimation

## Examples

``` r
species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)
chao2(species)
```
