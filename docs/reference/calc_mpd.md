# Calculate Mean Pairwise Distance (MPD)

Mean phylogenetic or functional distance between species in a sample.

## Usage

``` r
calc_mpd(
  dist_mat,
  species_present,
  abundance_weighted = FALSE,
  abundances = numeric()
)
```

## Arguments

- dist_mat:

  Pairwise distance matrix (species x species)

- species_present:

  Logical vector indicating which species are present

- abundance_weighted:

  If true, weight by abundances

- abundances:

  Vector of abundances (used if abundance_weighted = true)

## Value

MPD value
