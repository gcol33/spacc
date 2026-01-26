# Coleman Expected Accumulation

Compute the expected species accumulation curve using the Coleman method
(Coleman et al. 1982). This is an analytical formula, no simulation
needed.

## Usage

``` r
coleman(x)
```

## Arguments

- x:

  A site-by-species matrix (presence/absence or abundance).

## Value

A data.frame with columns: sites, expected, sd

## References

Coleman, B.D., Mares, M.A., Willig, M.R. & Hsieh, Y.H. (1982).
Randomness, area, and species richness. Ecology, 63, 1121-1133.
