# Sample Completeness Profile

Compute the ratio of observed to estimated diversity across diversity
orders, measuring how complete a sample is at each level of the Hill
number spectrum.

## Usage

``` r
completenessProfile(x, q = seq(0, 2, by = 0.2), coords = NULL)
```

## Arguments

- x:

  A site-by-species matrix (abundance data).

- q:

  Numeric vector. Orders of diversity to evaluate. Default
  `seq(0, 2, by = 0.2)`.

- coords:

  Optional data.frame with columns `x` and `y` for spatial mapping. When
  provided, enables `plot(type = "map")` and
  [`as_sf()`](https://gillescolling.com/spacc/reference/as_sf.md).

## Value

An object of class `spacc_completeness` containing:

- completeness:

  Named numeric vector of completeness ratios per q

- observed:

  Named numeric vector of observed Hill numbers per q

- estimated:

  Named numeric vector of estimated asymptotic Hill numbers per q

- per_site:

  Matrix of per-site completeness (sites x q values), or `NULL`

- q:

  Vector of diversity orders

- coords:

  Coordinates if provided

- n_sites:

  Number of sites

- n_species:

  Number of species

## Details

Sample completeness at order q is: \$\$C_q =
\frac{D_q^{obs}}{D_q^{est}}\$\$ where \\D_q^{obs}\\ is the observed Hill
number and \\D_q^{est}\\ is the estimated asymptotic Hill number.

Completeness near 1 means the sample captures most of the true diversity
at that order. Completeness typically increases with q because dominant
species are detected early.

Asymptotic estimators used:

- q = 0: Chao1 estimator

- q = 1: Chao & Jost (2015) entropy estimator

- q = 2: Inverse Simpson estimator with bias correction

- Other q: Interpolated between adjacent integer estimates

When `coords` is provided, per-site completeness is computed by treating
each site's abundance vector as an independent sample.

## References

Chao, A. & Jost, L. (2012). Coverage-based rarefaction and
extrapolation: standardizing samples by completeness rather than size.
Ecology, 93, 2533-2547.

Chao, A. & Jost, L. (2015). Estimating diversity and entropy profiles
via discovery rates of new species. Methods in Ecology and Evolution, 6,
873-882.

## See also

[`diversityProfile()`](https://gillescolling.com/spacc/reference/diversityProfile.md)
for observed diversity profiles,
[`chao1()`](https://gillescolling.com/spacc/reference/chao1.md) for
richness estimation

## Examples

``` r
species <- matrix(rpois(50 * 30, 2), nrow = 50)
comp <- completenessProfile(species)
print(comp)

# \donttest{
plot(comp)
# }
```
