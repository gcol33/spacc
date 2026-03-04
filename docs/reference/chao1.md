# Chao1 Richness Estimator

Estimate total species richness from abundance data using the Chao1
estimator (Chao 1984). Uses the number of singletons and doubletons to
estimate undetected species.

## Usage

``` r
chao1(x)
```

## Arguments

- x:

  A site-by-species matrix (abundance data). Columns are pooled across
  sites.

## Value

An object of class `spacc_estimate` with components:

- estimator:

  Name of the estimator (`"chao1"`)

- estimate:

  Estimated total richness

- se:

  Standard error of the estimate

- lower:

  Lower 95 percent confidence bound

- upper:

  Upper 95 percent confidence bound

- S_obs:

  Observed species richness

- details:

  List with `f1` (singletons) and `f2` (doubletons)

## Details

The Chao1 estimator is: \$\$S\_{Chao1} = S\_{obs} + \frac{f_1^2}{2
f_2}\$\$ where \\f_1\\ is the number of singletons (species with total
abundance 1) and \\f_2\\ is the number of doubletons (abundance 2).

When \\f_2 = 0\\, the bias-corrected form is used: \$\$S\_{Chao1} =
S\_{obs} + \frac{f_1 (f_1 - 1)}{2}\$\$

## References

Chao, A. (1984). Nonparametric estimation of the number of classes in a
population. Scandinavian Journal of Statistics, 11, 265-270.

Chao, A. (1987). Estimating the population size for capture-recapture
data with unequal catchability. Biometrics, 43, 783-791.

## See also

[`chao2()`](https://gillescolling.com/spacc/reference/chao2.md) for
incidence-based estimation,
[`ace()`](https://gillescolling.com/spacc/reference/ace.md) for
abundance-based coverage estimator

## Examples

``` r
species <- matrix(rpois(50 * 30, 2), nrow = 50)
chao1(species)
```
