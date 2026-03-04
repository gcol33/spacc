# Diversity Profile

Compute Hill numbers across a continuous range of diversity orders (q),
producing a diversity profile for each site and the regional pool.

## Usage

``` r
diversityProfile(
  x,
  q = seq(0, 3, by = 0.1),
  type = c("both", "per_site", "regional"),
  coords = NULL
)
```

## Arguments

- x:

  A site-by-species matrix (abundance data).

- q:

  Numeric vector. Orders of diversity to evaluate. Default
  `seq(0, 3, by = 0.1)`.

- type:

  Character. What to compute: `"per_site"` (per-site profiles),
  `"regional"` (pooled gamma), or `"both"` (default).

- coords:

  Optional data.frame with columns `x` and `y` for spatial mapping. When
  provided, enables `plot(type = "map")`.

## Value

An object of class `spacc_profile` containing:

- per_site:

  Matrix of per-site diversity (sites x q values), or `NULL`

- regional:

  Named numeric vector of gamma diversity per q, or `NULL`

- q:

  Vector of diversity orders used

- coords:

  Coordinates if provided

- n_sites:

  Number of sites

- n_species:

  Number of species

## Details

A diversity profile plots effective number of species as a function of
the order q. The key property is that Hill numbers are non-increasing in
q: \\D_q \ge D\_{q'}\\ for \\q \le q'\\.

- q = 0: Species richness (insensitive to abundance)

- q = 1: Exponential of Shannon entropy (all species weighted equally by
  their proportional abundance)

- q = 2: Inverse Simpson concentration (emphasizes dominant species)

- q \> 2: Increasingly dominated by common species

## References

Leinster, T. & Cobbold, C.A. (2012). Measuring diversity: the importance
of species similarity. Ecology, 93, 477-489.

Chao, A., Chiu, C.H. & Jost, L. (2014). Unifying species diversity,
phylogenetic diversity, functional diversity, and related similarity and
differentiation measures through Hill numbers. Annual Review of Ecology,
Evolution, and Systematics, 45, 297-324.

## See also

[`alphaDiversity()`](https://gillescolling.com/spacc/reference/alphaDiversity.md)
for per-site values at specific q,
[`gammaDiversity()`](https://gillescolling.com/spacc/reference/gammaDiversity.md)
for regional diversity,
[`evenness()`](https://gillescolling.com/spacc/reference/evenness.md)
for evenness profiles

## Examples

``` r
species <- matrix(rpois(50 * 30, 2), nrow = 50)
prof <- diversityProfile(species)
print(prof)

# \donttest{
plot(prof)
# }
```
