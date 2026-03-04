# Improved Chao2 (iChao2) Richness Estimator

Estimate total species richness from incidence data using the improved
Chao2 estimator (Chiu et al. 2014). Uses uniques through quadruplicates
(Q1–Q4) to reduce bias for small samples.

## Usage

``` r
iChao2(x)
```

## Arguments

- x:

  A site-by-species matrix (presence/absence or abundance).
  Automatically binarized.

## Value

An object of class `spacc_estimate` with components:

- estimator:

  Name of the estimator (`"iChao2"`)

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

  List with `Q1`, `Q2`, `Q3`, `Q4`, `n_sites`

## Details

The improved Chao2 estimator is the incidence-based analogue of iChao1:
\$\$S\_{iChao2} = S\_{Chao2} + \frac{Q_3}{4 Q_4} \max\left(Q_1 -
\frac{Q_2 Q_3}{2 Q_4}, 0\right)\$\$

When \\Q_4 = 0\\, the estimator collapses to Chao2.

## References

Chiu, C.H., Wang, Y.T., Walther, B.A. & Chao, A. (2014). An improved
nonparametric lower bound of species richness via a modified Good-Turing
frequency formula. Biometrics, 70, 671-682.

## See also

[`chao2()`](https://gillescolling.com/spacc/reference/chao2.md) for the
standard estimator,
[`iChao1()`](https://gillescolling.com/spacc/reference/iChao1.md) for
abundance-based version

## Examples

``` r
species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)
iChao2(species)
```
