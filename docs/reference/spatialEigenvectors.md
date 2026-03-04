# Spatial Eigenvectors (PCNM/dbMEM)

Compute spatial eigenvectors from site coordinates using Principal
Coordinates of Neighbour Matrices (PCNM) or distance-based Moran's
Eigenvector Maps (dbMEM).

## Usage

``` r
spatialEigenvectors(
  coords,
  threshold = NULL,
  method = c("pcnm", "dbmem"),
  distance = c("euclidean", "haversine")
)
```

## Arguments

- coords:

  A data.frame with columns `x` and `y`, or a `spacc_dist` object.

- threshold:

  Numeric. Distance threshold for truncation. Default `NULL` uses the
  largest edge in the minimum spanning tree.

- method:

  Character. `"pcnm"` (default) retains all positive eigenvectors;
  `"dbmem"` retains only those with positive Moran's I.

- distance:

  Character. `"euclidean"` (default) or `"haversine"`.

## Value

An object of class `spacc_mem` containing:

- vectors:

  Matrix of eigenvectors (sites x n_vectors)

- eigenvalues:

  Positive eigenvalues

- moran_i:

  Moran's I for each eigenvector

- threshold:

  Truncation distance used

- coords:

  Original coordinates

- n_sites:

  Number of sites

- method:

  Method used

## Details

PCNM (Borcard & Legendre 2002) decomposes spatial structure into
orthogonal eigenvectors representing patterns at different spatial
scales. Large eigenvalues correspond to broad-scale patterns (positive
spatial autocorrelation), while small eigenvalues represent fine-scale
patterns.

The algorithm:

1.  Compute pairwise distances

2.  Truncate distances beyond threshold (set to 4 \* threshold)

3.  Double-centre the squared truncated distance matrix

4.  Extract eigenvectors with positive eigenvalues

5.  For `"dbmem"`: additionally filter to Moran's I \> 0

## References

Borcard, D. & Legendre, P. (2002). All-scale spatial analysis of
ecological data by means of principal coordinates of neighbour matrices.
Ecological Modelling, 153, 51-68.

Dray, S., Legendre, P. & Peres-Neto, P.R. (2006). Spatial modelling: a
comprehensive framework for principal coordinate analysis of neighbour
matrices (PCNM). Ecological Modelling, 196, 483-493.

## See also

[`spatialPartition()`](https://gillescolling.com/spacc/reference/spatialPartition.md)
for variance partitioning with MEMs

## Examples

``` r
coords <- data.frame(x = runif(30), y = runif(30))
mem <- spatialEigenvectors(coords)
print(mem)
```
