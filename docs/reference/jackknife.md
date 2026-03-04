# Jackknife Richness Estimator

First- or second-order jackknife estimator of total species richness
(Burnham & Overton 1978, 1979).

## Usage

``` r
jackknife(x, order = 1L)
```

## Arguments

- x:

  A site-by-species matrix (presence/absence or abundance).
  Automatically binarized.

- order:

  Integer. Jackknife order: 1 (default) or 2.

## Value

An object of class `spacc_estimate`.

## Details

First-order jackknife: \$\$S\_{jack1} = S\_{obs} + Q_1 \frac{n-1}{n}\$\$

Second-order jackknife: \$\$S\_{jack2} = S\_{obs} +
\frac{Q_1(2n-3)}{n} - \frac{Q_2(n-2)^2}{n(n-1)}\$\$

## References

Burnham, K.P. & Overton, W.S. (1978). Estimation of the size of a closed
population when capture probabilities vary among animals. Biometrika,
65, 625-633.

Burnham, K.P. & Overton, W.S. (1979). Robust estimation of population
size when capture probabilities vary among animals. Ecology, 60,
927-936.

## See also

[`chao2()`](https://gillescolling.com/spacc/reference/chao2.md),
[`bootstrap_richness()`](https://gillescolling.com/spacc/reference/bootstrap_richness.md)

## Examples

``` r
species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)
jackknife(species, order = 1)
jackknife(species, order = 2)
```
