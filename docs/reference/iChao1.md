# Improved Chao1 (iChao1) Richness Estimator

Estimate total species richness using the improved Chao1 estimator (Chiu
et al. 2014). Uses singletons through quadrupletons (f1–f4) to reduce
bias for small samples.

## Usage

``` r
iChao1(x)
```

## Arguments

- x:

  A site-by-species matrix (abundance data). Columns are pooled across
  sites.

## Value

An object of class `spacc_estimate` with components:

- estimator:

  Name of the estimator (`"iChao1"`)

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

  List with `f1`, `f2`, `f3`, `f4`

## Details

The improved Chao1 estimator adds a correction term using f3 and f4:
\$\$S\_{iChao1} = S\_{Chao1} + \frac{f_3}{4 f_4} \max\left(f_1 -
\frac{f_2 f_3}{2 f_4}, 0\right)\$\$

When \\f_4 = 0\\, the estimator collapses to Chao1.

## References

Chiu, C.H., Wang, Y.T., Walther, B.A. & Chao, A. (2014). An improved
nonparametric lower bound of species richness via a modified Good-Turing
frequency formula. Biometrics, 70, 671-682.

## See also

[`chao1()`](https://gillescolling.com/spacc/reference/chao1.md) for the
standard estimator,
[`iChao2()`](https://gillescolling.com/spacc/reference/iChao2.md) for
incidence-based version

## Examples

``` r
species <- matrix(rpois(50 * 30, 2), nrow = 50)
iChao1(species)
```
