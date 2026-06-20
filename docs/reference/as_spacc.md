# Convert an Accumulation Curve from Another Package to a spacc Object

Bring an existing species accumulation curve into the spacc object
system so it can be summarised, plotted, and compared with
[`compare()`](https://gillescolling.com/spacc/reference/compare.md)
alongside spatial curves. A vegan `specaccum` object is currently
supported.

## Usage

``` r
as_spacc(x, ...)

# Default S3 method
as_spacc(x, ...)

# S3 method for class 'specaccum'
as_spacc(x, ...)
```

## Arguments

- x:

  An object to convert. A `specaccum` object from
  [`vegan::specaccum()`](https://vegandevs.github.io/vegan/reference/specaccum.html)
  is supported.

- ...:

  Ignored, for method extensibility.

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
`method = "exact"` or `"collector"`) contributes its single mean curve.

Coordinates cannot be recovered from a `specaccum` object, so the result
carries no `coords`: it is a curve container for plotting and
comparison, not a re-runnable spatial object.

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
```
