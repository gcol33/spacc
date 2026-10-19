# Convert an Accumulation Curve from Another Package to a spacc Object

Bring an existing species accumulation curve into the spacc object
system so it can be summarised, plotted, extrapolated and compared with
[`compare()`](https://gillescolling.com/spacc/reference/compare.md)
alongside spatial curves.

## Usage

``` r
as_spacc(x, ...)

# Default S3 method
as_spacc(x, ...)

# S3 method for class 'specaccum'
as_spacc(x, ...)

# S3 method for class 'numeric'
as_spacc(x, ...)

# S3 method for class 'data.frame'
as_spacc(x, column = NULL, ...)
```

## Arguments

- x:

  An object to convert: a `specaccum` object from
  [`vegan::specaccum()`](https://vegandevs.github.io/vegan/reference/specaccum.html),
  a numeric vector holding one accumulation curve, or a data frame with
  such a curve in one column (for example the output of
  [`Rarefy::directionalSAC()`](https://rdrr.io/pkg/Rarefy/man/directionalSAC.html)
  or
  [`Rarefy::rare_alpha()`](https://rdrr.io/pkg/Rarefy/man/rare_alpha.html)).

- ...:

  Ignored, for method extensibility.

- column:

  For a data frame, the name of the column holding the curve. `NULL`
  (default) uses `N_SCR` (the spatially constrained curve of
  [`Rarefy::directionalSAC()`](https://rdrr.io/pkg/Rarefy/man/directionalSAC.html))
  or `Rarefaction`
  ([`Rarefy::rare_alpha()`](https://rdrr.io/pkg/Rarefy/man/rare_alpha.html)),
  whichever is present.

## Value

An object of class `spacc` with
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html),
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), and
[`compare()`](https://gillescolling.com/spacc/reference/compare.md)
support.

## Details

A `specaccum` object built with `method = "random"` stores one richness
curve per permutation in its `perm` matrix. Those curves become the rows
of the returned `curves` matrix and behave exactly like spacc seeds, so
the across-curve confidence band and the
[`compare()`](https://gillescolling.com/spacc/reference/compare.md)
permutation, bootstrap, and AUC tests carry over unchanged. A
`specaccum` object without a permutation matrix (for example
`method = "exact"` or `"collector"`), a numeric vector and a data frame
column each contribute a single mean curve, entry k being the value
after k sites.

`Rarefy::directionalSAC(community, gradient)` with a site distance
matrix as `gradient` orders the sites by distance from each plot in turn
and averages the collector curves; this is
`spacc(x, coords, method = "knn", focal_points = coords)`, so the two
curves agree up to the order of sites at tied distances.

Coordinates cannot be recovered from these objects, so the result
carries no `coords`: it is a curve container for plotting and
comparison, not a re-runnable spatial object.

## References

Thouverai, E., Pavoine, S., Tordoni, E., Chiarucci, A., Ricotta, C.,
Rocchini, D. & Bacaro, G. (2026). Rarefy: an R package for the
calculation of taxonomic, functional and phylogenetic sample-based
spatially explicit rarefaction curves. Ecological Informatics, 103803.
[doi:10.1016/j.ecoinf.2026.103803](https://doi.org/10.1016/j.ecoinf.2026.103803)

## See also

[`compare()`](https://gillescolling.com/spacc/reference/compare.md) to
test a converted curve against a spatial one.

## Examples

``` r
# \donttest{
if (requireNamespace("vegan", quietly = TRUE)) {
  data(BCI, package = "vegan")
  classic <- vegan::specaccum(BCI, method = "random")
  sp <- as_spacc(classic)
  summary(sp)
  plot(sp)
}
# }
as_spacc(c(5, 8, 10, 11, 12))
```
