# Calculate Faith's PD from branch lengths

Sum of branch lengths connecting species in sample to root. Requires
edge matrix and edge lengths from ape::phylo object.

## Usage

``` r
calc_faith_pd(edge, edge_length, n_tips, species_present)
```

## Arguments

- edge:

  Integer matrix (n_edges x 2) of tree edges

- edge_length:

  Numeric vector of edge lengths

- n_tips:

  Number of tips in tree

- species_present:

  Logical vector indicating which tips are present

## Value

Faith's PD
