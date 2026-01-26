# Calculate Mean Nearest Taxon Distance (MNTD)

Mean distance to nearest relative for each species in sample.

## Usage

``` r
calc_mntd(
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

  Vector of abundances

## Value

MNTD value
