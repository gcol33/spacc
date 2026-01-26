# Single kNN Functional Diversity Accumulation

Single kNN Functional Diversity Accumulation

## Usage

``` r
cpp_func_knn_single(species_mat, site_dist_mat, traits, seed, metrics)
```

## Arguments

- species_mat:

  Integer matrix (sites x species) abundances

- site_dist_mat:

  Site distance matrix

- traits:

  Trait matrix (species x traits)

- seed:

  Starting site

- metrics:

  Metrics: "fdis", "fric"

## Value

List with curves for each metric
