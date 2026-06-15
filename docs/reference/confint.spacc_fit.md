# Confidence interval for an extrapolated asymptote

Returns the calibrated bootstrap interval stored on the fit by default.
With `method = "profile"` the (over-confident) `nls` profile interval is
computed instead. The interval is for the asymptote *parameter*; see the
bias caveat in
[`extrapolate()`](https://gillescolling.com/spacc/reference/extrapolate.md).

## Usage

``` r
# S3 method for class 'spacc_fit'
confint(object, parm, level = 0.95, method = c("stored", "profile"), ...)
```

## Arguments

- object:

  A `spacc_fit` object.

- parm:

  Ignored; the asymptote (`"a"`) is always returned.

- level:

  Confidence level. Recomputed from stored bootstrap draws.

- method:

  `"stored"` (default) or `"profile"`.

- ...:

  Passed to [`stats::confint()`](https://rdrr.io/r/stats/confint.html)
  when `method = "profile"`.

## Value

A 1-row matrix with the lower and upper bound for `a`.
