# Exact (Mao Tau) Expected Accumulation

Compute the expected species accumulation curve using sample-based
rarefaction (Mao Tau estimator). This is analytically identical to the
expected curve from random permutations.

## Usage

``` r
mao_tau(x)
```

## Arguments

- x:

  A site-by-species matrix (presence/absence or abundance).

## Value

A data.frame with columns: sites, expected, sd, lower, upper

## References

Colwell, R.K., Mao, C.X. & Chang, J. (2004). Interpolating,
extrapolating, and comparing incidence-based species accumulation
curves. Ecology, 85, 2717-2727.
