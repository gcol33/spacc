# Chao2 Richness Estimator

Estimate total species richness from incidence (presence/absence) data
using the Chao2 estimator (Chao 1987).

## Usage

``` r
chao2(x, n_units = NULL)
```

## Arguments

- x:

  A site-by-species matrix (presence/absence or abundance).
  Automatically binarized.

- n_units:

  Total number of sampling units in the study region, when the rows of
  `x` are a sample drawn without replacement from a finite set of units
  (atlas cells, a plot network). `NULL` (default) or `Inf` assumes
  sampling with replacement and gives Chao2; a finite value gives the
  MoRE estimator of Chiu (2023).

## Value

An object of class `spacc_estimate`.

## Details

The Chao2 estimator is the incidence-based analogue of Chao1:
\$\$S\_{Chao2} = S\_{obs} + \frac{t-1}{t} \frac{Q_1^2}{2 Q_2}\$\$ where
\\t\\ is the number of sites, \\Q_1\\ is the number of uniques (species
found at exactly 1 site) and \\Q_2\\ is the number of duplicates
(species found at exactly 2 sites).

Chao2 assumes the sites are drawn with replacement, or from an
effectively infinite set of units. When they are a large fraction \\w =
t / T\\ of a finite set of \\T\\ units it overestimates richness and
does not reach the true richness as \\w \to 1\\. With `n_units = T`,
`chao2()` returns the moment-based estimator for sampling without
replacement (MoRE, Chiu 2023, Eq. 5): \$\$S\_{MoRE} = S\_{obs} +
\frac{t-1}{t} \frac{(1-w) Q_1^2}{Q_2 + w (Q_1 - Q_2)^+}\$\$ where
\\(A)^+\\ is \\A\\ for \\A \> 0\\ and 1 otherwise; for \\Q_2 = 0\\ the
undetected term is \\\frac{t-1}{t} \frac{1-w}{w} (Q_1 - 1)\\. MoRE
equals the observed richness at \\w = 1\\. Its standard error is the
delta-method variance of Chiu (2023) and the interval is the log-normal
one used for Chao2.

## References

Chao, A. (1987). Estimating the population size for capture-recapture
data with unequal catchability. Biometrics, 43, 783-791.

Chiu, C.-H. (2023). A species richness estimator for sample-based
incidence data sampled without replacement. Methods in Ecology and
Evolution, 14, 2482-2493.
[doi:10.1111/2041-210X.14146](https://doi.org/10.1111/2041-210X.14146)

## See also

[`chao1()`](https://gillescolling.com/spacc/reference/chao1.md) for
abundance-based estimation

## Examples

``` r
species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)
chao2(species)
# the 50 sites are 80% of a 62-cell grid
chao2(species, n_units = 62)
```
