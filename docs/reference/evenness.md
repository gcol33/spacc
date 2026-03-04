# Evenness Profiles

Compute species evenness across sites using Hill-based, Pielou, or
Simpson evenness measures.

## Usage

``` r
evenness(
  x,
  q = seq(0.1, 3, by = 0.1),
  type = c("hill", "pielou", "simpson"),
  coords = NULL
)
```

## Arguments

- x:

  A site-by-species matrix (abundance data).

- q:

  Numeric vector. Orders of diversity for Hill evenness. Default
  `seq(0.1, 3, by = 0.1)`. Note: q = 0 is excluded by default because
  Hill evenness at q = 0 is trivially S/S = 1.

- type:

  Character. Evenness type: `"hill"` (Hill evenness E_q = D_q / D_0,
  default), `"pielou"` (Pielou's J = log(D_1) / log(S)), or `"simpson"`
  (Simpson evenness = (1/D_2) / S).

- coords:

  Optional data.frame with columns `x` and `y` for spatial mapping. When
  provided, enables `plot(type = "map")`.

## Value

An object of class `spacc_evenness` containing:

- per_site:

  Matrix or vector of per-site evenness values

- regional:

  Regional (pooled) evenness

- q:

  Orders used (for Hill type)

- type:

  Evenness type

- coords:

  Coordinates if provided

- n_sites:

  Number of sites

- n_species:

  Number of species

## Details

All evenness measures are bounded in \\\[0, 1\]\\:

- 0 = maximally uneven (one dominant species)

- 1 = perfectly even (all species equally abundant)

**Hill evenness** (Jost 2010): \$\$E_q = D_q / D_0\$\$ This divides the
effective number of species at order q by species richness.

**Pielou's J** (Pielou 1966): \$\$J = \frac{\log(D_1)}{\log(S)} =
\frac{H'}{\log(S)}\$\$

**Simpson evenness**: \$\$E\_{1/D} = \frac{1}{D_2 \cdot S}\$\$

## References

Jost, L. (2010). The relation between evenness and diversity. Diversity,
2, 207-232.

Pielou, E.C. (1966). The measurement of diversity in different types of
biological collections. Journal of Theoretical Biology, 13, 131-144.

## See also

[`diversityProfile()`](https://gillescolling.com/spacc/reference/diversityProfile.md)
for Hill number profiles,
[`alphaDiversity()`](https://gillescolling.com/spacc/reference/alphaDiversity.md)
for raw diversity values

## Examples

``` r
species <- matrix(rpois(50 * 30, 2), nrow = 50)

# Hill evenness profile
ev <- evenness(species)
print(ev)

# Pielou's J
ev_j <- evenness(species, type = "pielou")
print(ev_j)
```
