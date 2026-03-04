# Functional Diversity Profile

Compute functional Hill numbers (Leinster & Cobbold 2012) across a
continuous range of diversity orders (q), producing a functional
diversity profile based on trait similarity.

## Usage

``` r
diversityProfileFunc(
  x,
  traits,
  q = seq(0, 3, by = 0.1),
  type = c("both", "per_site", "regional"),
  dist_method = c("euclidean", "gower"),
  normalize = TRUE,
  coords = NULL
)
```

## Arguments

- x:

  A site-by-species matrix (abundance data). Column names must match row
  names in `traits`.

- traits:

  A data.frame of species traits. Row names must match column names in
  `x`.

- q:

  Numeric vector. Orders of diversity. Default `seq(0, 3, by = 0.1)`.

- type:

  Character. What to compute: `"per_site"`, `"regional"`, or `"both"`
  (default).

- dist_method:

  Character. Distance method for trait matrix: `"euclidean"` (default)
  or `"gower"`.

- normalize:

  Logical. Normalize distances to \[0, 1\]? Default `TRUE`.

- coords:

  Optional data.frame with `x` and `y` for spatial mapping.

## Value

An object of class `spacc_profile` with `$profile_type = "functional"`.

## Details

Functional Hill numbers (Leinster & Cobbold 2012) incorporate trait
similarity via a similarity matrix Z = 1 - D. When all species are
maximally dissimilar (Z = identity), this reduces to standard Hill
numbers.

## References

Leinster, T. & Cobbold, C.A. (2012). Measuring diversity: the importance
of species similarity. Ecology, 93, 477-489.

## See also

[`diversityProfile()`](https://gillescolling.com/spacc/reference/diversityProfile.md)
for taxonomic profiles,
[`diversityProfilePhylo()`](https://gillescolling.com/spacc/reference/diversityProfilePhylo.md)
for phylogenetic profiles

## Examples

``` r
# \donttest{
species <- matrix(rpois(20 * 10, 2), nrow = 20,
                  dimnames = list(NULL, paste0("sp", 1:10)))
traits <- data.frame(
  body_size = rnorm(10), diet = rnorm(10),
  row.names = paste0("sp", 1:10)
)
prof <- diversityProfileFunc(species, traits)
print(prof)
# }
```
