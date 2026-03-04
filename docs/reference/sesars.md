# Sampling Effort Species-Area Relationship (SESARS)

Model the joint effect of sampling effort and area on species richness.
Corrects for unequal survey intensity across sites, common in atlas data
and citizen science datasets.

## Usage

``` r
sesars(object, effort, model = c("power", "additive"), ...)
```

## Arguments

- object:

  A `spacc` object.

- effort:

  Numeric vector. Sampling effort per site (e.g., hours, visits,
  trap-nights). Must have length equal to number of sites.

- model:

  Character. SESARS model:

  - `"power"` (default): S = c \* A^z \* E^w (multiplicative power law)

  - `"additive"`: S = c + z \* log(A) + w \* log(E)

- ...:

  Additional arguments passed to
  [`stats::nls()`](https://rdrr.io/r/stats/nls.html) or
  [`stats::lm()`](https://rdrr.io/r/stats/lm.html).

## Value

An object of class `spacc_sesars` containing:

- model:

  Model type

- fit:

  Fitted model object

- coef:

  Model coefficients

- data:

  Data frame used for fitting

## Details

Standard SARs assume complete sampling within each area unit. SESARS
incorporates sampling effort (E) alongside area (A) to provide unbiased
richness estimates across regions with unequal survey intensity.

## References

Dennstadt, F., Horak, J. & Martin, M.D. (2019). Predictive sampling
effort and species-area relationship models for estimating richness in
fragmented landscapes. Diversity and Distributions, 26, 1112-1123.

## See also

[`extrapolate()`](https://gillescolling.com/spacc/reference/extrapolate.md),
[`spacc()`](https://gillescolling.com/spacc/reference/spacc.md)

## Examples

``` r
# \donttest{
coords <- data.frame(x = runif(50), y = runif(50))
species <- matrix(rbinom(50 * 30, 1, 0.3), nrow = 50)
sac <- spacc(species, coords)
effort <- rpois(50, 10)
ses <- sesars(sac, effort, model = "power")
print(ses)
# }
```
