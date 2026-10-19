# Accumulation Rate Curve (ARC) Richness Estimator

Estimate total species richness from the rate at which new species
appear along a random reordering of the observations (Shestopaloff
2025). The rate curve is extrapolated with three decay curves whose
tails give the number of species not yet seen, combined with bounded
weights.

## Usage

``` r
arc(
  x,
  data = c("abundance", "incidence"),
  K = 20L,
  adjust = FALSE,
  n_boot = 2000L,
  parallel = TRUE,
  n_cores = NULL
)
```

## Arguments

- x:

  A site-by-species matrix. For `data = "abundance"` the columns are
  pooled into species abundances; for `data = "incidence"` the matrix is
  binarized and each species counts the sites it occupies.

- data:

  Character. `"abundance"` (default) treats every individual as one
  observation; `"incidence"` treats every species-site occurrence as
  one.

- K:

  Integer. Species seen more than `K` times are left out of the rate
  curve. Default 20.

- adjust:

  Logical. Apply the smoothing adjustment of Eq. 11 with a log-linear
  fit to the frequency counts? Default `FALSE`. See Details.

- n_boot:

  Integer. Parametric bootstrap draws for the interval. Default 2000.
  `0` returns the point estimate only.

- parallel:

  Logical. Use parallel processing? Default `TRUE`.

- n_cores:

  Integer. Number of cores.

## Value

An object of class `spacc_estimate` with `estimator = "ARC"`. Its
`details` hold `undetected` (undetected species from each curve) and
`weights` (their weights), both averaged over the curve portions, `K`,
`data`, `penalty` (the selected Poisson-mixture penalty) and `boot` (the
bootstrap estimates).

## Details

Let \\y_k\\ be the number of species observed \\k\\ times and \\N\\ the
number of observations of the species with \\k \le K\\. Drawing the
observations in random order, the expected number of these species seen
after \\n\\ draws is \$\$S(n) = \sum\_{k=1}^{K} y_k \left\[1 - (1 -
n/N)^k\right\],\$\$ and the accumulation rate is \\s_n = S(n) -
S(n-1)\\, evaluated at up to 1000 points. Three curves are fitted to the
tail of \\s_n\\ by least squares with weights proportional to the
squared point index: the log-normal \\1 - \Phi((\log n - \mu)/\sigma)\\,
the log-logistic \\1/(1 + \alpha n^\gamma)\\ and the Weibull
\\\exp(-\lambda n^\gamma)\\, the last on the log scale. Each curve's
integral from \\N\\ to infinity is its number of undetected species; the
log-logistic integral is finite only for \\\gamma \> 1\\.

The estimate combines the three with weights fitted by the squared
misfit to \\s_n\\ plus a penalty on the product of the log-logistic and
Weibull weights (Eq. 8), each bounded by caps built from the three
estimates and the observed richness (Eq. 7). The undetected species are
averaged over fits to the last 100, 95, ..., 70 percent of the curve.

The penalty is there to make the fit prefer either the log-logistic or
the Weibull, and gives the loss one local minimum for each. The weights
are the end point of the projected gradient flow of the loss started at
zero weights, which `arc()` solves in closed form face by face of the
box of admissible weights. This is the minimum the author's reference
code reaches (L-BFGS within an augmented Lagrangian, started at zero).
In Poisson-Gamma samples of intermediate sparsity it is the log-logistic
one, with median bias between -1 and +4 percent; the global minimum of
the loss is the Weibull one there, with median bias between -9 and -17
percent
([`vignette("richness-estimation")`](https://gillescolling.com/spacc/articles/richness-estimation.md)).

With `adjust = TRUE` the counts \\y_1, \dots, y\_{k_1}\\ (the first run
of positive counts) are mixed half and half with a log-linear fit \\\log
y_k = b_0 + b_1 k\\, and fully replaced by it; the adjusted estimate is
\\2 S\_{0.5} - S_1\\. The adjustment assumes the ratios \\(k+1) y\_{k+1}
/ y_k\\ are linear in \\k\\.

The interval comes from a parametric bootstrap of the counts. Their
covariance is that of a mixture of 30 Poisson distributions whose
weights are fitted to \\(y_1, \dots, y_K, y\_{K+})\\ by penalised least
squares (Eq. 12-15). The penalty is chosen by resampling species: 1500
resamples are ordered by their distance to the left-out species and
split into deciles, and the first decile whose out-of-sample error is
smallest at a non-zero penalty sets it. `n_boot` count vectors are drawn
from a normal distribution with that covariance, rounded to non-negative
integers and re-estimated; the interval is their 2.5 and 97.5 percent
quantiles and the standard error their standard deviation. In
Poisson-Gamma simulations the 95 percent interval covered the truth in
93 to 96 percent of samples for sparse data and about 89 percent for
dense data, where the point estimate is biased upwards
([`vignette("richness-estimation")`](https://gillescolling.com/spacc/articles/richness-estimation.md)).

With `data = "incidence"` the observations are species-site occurrences
and \\y_k\\ becomes \\Q_k\\, the number of species found at \\k\\ sites,
the substitution that turns Chao1 into Chao2.

## References

Shestopaloff, K. (2025). An accumulation rate curve estimator for total
species. Environmental and Ecological Statistics, 32, 311-345.
[doi:10.1007/s10651-025-00644-y](https://doi.org/10.1007/s10651-025-00644-y)

## See also

[`chao1()`](https://gillescolling.com/spacc/reference/chao1.md),
[`chao2()`](https://gillescolling.com/spacc/reference/chao2.md),
[`ace()`](https://gillescolling.com/spacc/reference/ace.md)

## Examples

``` r
set.seed(1)
species <- matrix(rpois(40 * 300, rgamma(300, 0.4, 2)), nrow = 40, byrow = TRUE)
arc(species, n_boot = 0)
# \donttest{
arc(species, n_boot = 200)
# }
```
