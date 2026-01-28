# Species-Fragmented Area Relationship (SFAR)

Separate the effects of habitat loss (area reduction) from fragmentation
(splitting into patches) on species richness. Extends the classic
power-law SAR with an explicit fragmentation term.

## Usage

``` r
sfar(object, patches, model = c("power", "log"), ...)
```

## Arguments

- object:

  A `spacc` object.

- patches:

  Factor or integer vector assigning each site to a habitat fragment
  (patch). Must have length equal to the number of sites.

- model:

  Character. SFAR model:

  - `"power"` (default): S = c \* A^z \* n^(-f)

  - `"log"`: log(S) = log(c) + z*log(A) - f*log(n)

- ...:

  Additional arguments.

## Value

An object of class `spacc_sfar` containing:

- fit:

  Fitted model object

- coef:

  Coefficients: c (intercept), z (area exponent), f (fragmentation
  exponent)

- n_patches:

  Number of habitat fragments

## Details

The SFAR (Hanski et al. 2013) extends the power-law SAR to quantify the
additional effect of habitat fragmentation on species richness. The
model S = c \* A^z \* n^(-f) adds a penalty term for fragmentation (n =
number of fragments), where f \> 0 indicates that fragmentation reduces
richness beyond what area loss alone would predict.

## References

Hanski, I., Zurita, G.A., Bellocq, M.I. & Rybicki, J. (2013).
Species-fragmented area relationship. Proceedings of the National
Academy of Sciences, 110, 12715-12720.

Rybicki, J. & Hanski, I. (2013). Species-area relationships and
extinctions caused by habitat loss and fragmentation. Ecology Letters,
16, 27-38.

## See also

[`extrapolate()`](https://gillescolling.com/spacc/reference/extrapolate.md),
[`sesars()`](https://gillescolling.com/spacc/reference/sesars.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sac <- spacc(species, coords)
patches <- kmeans(coords, centers = 5)$cluster
sfar_result <- sfar(sac, patches)
print(sfar_result)
plot(sfar_result)
} # }
```
