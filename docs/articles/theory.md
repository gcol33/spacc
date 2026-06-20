# Methods and Algorithms

## Overview

This vignette describes what
[`spacc()`](https://gillescolling.com/spacc/reference/spacc.md) computes
and how. A classical species accumulation curve draws sampling units in
random order; `spacc` draws them in spatial order, expanding outward
from a focal site so the curve reflects how richness grows with surveyed
*area* rather than with effort alone. The expansion rule is a choice,
and the package offers seven of them. Each rule defines an ordering of
the sites; the curve at step \\k\\ is the number of distinct species
seen in the first \\k\\ sites of that ordering. Repeating the expansion
from many starting sites turns the spread of the curves into a
confidence band.

The sections below define the vocabulary, state the accumulation problem
formally, derive each expansion method from its ordering rule, give the
distance metrics those rules depend on, map the mathematics onto the
function arguments, and explain the two-tier nearest-neighbour backend
that keeps the computation fast on large datasets. The applied vignettes
build on this foundation: see
[`vignette("quickstart")`](https://gillescolling.com/spacc/articles/quickstart.md)
for the basic workflow,
[`vignette("diversity")`](https://gillescolling.com/spacc/articles/diversity.md)
for the Hill, beta, and coverage extensions, and
[`vignette("extrapolation")`](https://gillescolling.com/spacc/articles/extrapolation.md)
for fitting asymptotic models to the curves.

### Contents

1.  [Terminology](#terminology)

2.  [The spatial walk, concretely](#the-spatial-walk-concretely)

3.  [Problem formulation](#problem-formulation)

4.  [Expansion methods](#expansion-methods)

5.  [Distance metrics](#distance-metrics)

6.  [From theory to implementation](#from-theory-to-implementation)

7.  [The two-tier backend](#the-two-tier-backend)

8.  [Seeds and uncertainty](#seeds-and-uncertainty)

9.  [Relationship to random-order
    curves](#relationship-to-random-order-curves)

10. [From curves to other diversity
    measures](#from-curves-to-other-diversity-measures)

11. [Design notes](#design-notes)

12. [References](#references)

------------------------------------------------------------------------

## Terminology

These terms recur throughout the package documentation.

**Site.** One sampling unit: a plot, quadrat, grid cell, or locality.
Sites carry two coordinates and a row of species records.

**Occurrence matrix.** The input `x`: an \\n \times m\\ matrix of \\n\\
sites by \\m\\ species. Entries are presence-absence (\\0/1\\) or
abundance (counts). The accumulation algorithms reduce abundance to
presence-absence internally, since a curve counts whether a species is
*new*, not how many individuals it has.

**Coordinates.** The input `coords`: an \\n \times 2\\ table of
positions, supplied as a data frame with `x` and `y` columns, an `sf`
point layer, or a precomputed `spacc_dist` object.

**Ordering.** A permutation \\\pi = (\pi_1, \dots, \pi_n)\\ of the site
indices that fixes the sequence in which sites enter the curve. Every
expansion method is a rule for producing an ordering.

**Accumulation curve.** The vector \\S(1), \dots, S(n)\\ where \\S(k)\\
is the number of distinct species recorded across the first \\k\\ sites
of an ordering. The curve is non-decreasing and saturates at the total
species count.

**Seed.** The starting site \\\pi_1\\. Spatial methods grow an ordering
outward from the seed, so the seed determines the curve.

**Neighbourhood.** The set of unvisited sites considered for the next
step. For nearest-neighbour expansion it is all unvisited sites; for the
cone method it is restricted to a directional wedge.

**Seed band.** The set of curves obtained by repeating an expansion from
many seeds. Its pointwise quantiles form the confidence band.

**Backend.** The data structure used to answer nearest-neighbour
queries: a precomputed distance matrix (exact) or a spatial index (k-d
tree or ball tree).

------------------------------------------------------------------------

## The spatial walk, concretely

Start with twenty sites and watch one rule build an ordering. The
k-nearest-neighbour walk begins at a seed, steps to the closest
unvisited site, then to the closest site to *that*, and so on. The path
threads through the point cloud.

``` r

pts <- data.frame(
  x = c(1, 2, 1.5, 3, 3.5, 5, 5.5, 6, 2, 4,
        7, 7.5, 8, 1, 6.5, 4.5, 2.5, 8.5, 5, 3),
  y = c(1, 1.2, 2, 1, 2.5, 1, 2, 3, 3.5, 4,
        1.5, 3, 2, 4.5, 4, 3.5, 4.8, 4, 5, 5.5)
)
```

The walk is easy to trace by hand. From the seed, repeatedly pick the
nearest site that has not yet been visited.

``` r

knn_order <- function(coords, seed) {
  n <- nrow(coords); visited <- logical(n); ord <- integer(n)
  cur <- seed; visited[cur] <- TRUE; ord[1] <- cur
  for (k in 2:n) {
    d <- sqrt((coords$x - coords$x[cur])^2 + (coords$y - coords$y[cur])^2)
    d[visited] <- Inf
    cur <- which.min(d); visited[cur] <- TRUE; ord[k] <- cur
  }
  ord
}
ord <- knn_order(pts, seed = 1)
```

Drawing the path shows how the ordering hugs the local structure: the
first steps stay inside one cluster before the walk is forced to jump to
the next.

``` r

plot(pts$x, pts$y, pch = 19, col = "grey70", cex = 1.4,
     xlab = "x", ylab = "y", main = "kNN walk from seed (filled = seed)")
lines(pts$x[ord], pts$y[ord], col = "#2E7D32", lwd = 2)
points(pts$x[1], pts$y[1], pch = 19, col = "#C62828", cex = 2)
text(pts$x, pts$y, labels = match(seq_len(nrow(pts)), ord), pos = 3, cex = 0.7)
```

![Twenty points in the plane with a path connecting them in
nearest-neighbour visiting order, starting from the lower-left seed. The
path stays within local clusters before jumping to reach distant
points.](theory_files/figure-html/toy-knn-plot-1.svg)

A different seed gives a different path and a different curve. That
dependence on the starting point is not noise to be removed; it is the
spatial signal. When species are aggregated, a walk that starts inside a
rich cluster climbs steeply, a walk that starts in a sparse corner
climbs slowly, and the gap between those trajectories measures how much
composition turns over across the map. The package runs the walk from
many seeds precisely to capture that spread.

------------------------------------------------------------------------

## Problem formulation

**Input.** An occurrence matrix \\X \in \\0,1,2,\dots\\^{n \times m}\\
for \\n\\ sites and \\m\\ species, and coordinates \\c_1, \dots, c_n \in
\mathbb{R}^2\\. Let \\\mathrm{sp}(i) = \\\\ j : X\_{ij} \> 0 \\\\\\ be
the set of species present at site \\i\\.

**Ordering.** A method produces a permutation \\\pi\\ of \\\\1, \dots,
n\\\\. For the stochastic methods \\\pi\\ is a random variable whose
distribution depends on the seed and the rule.

**Accumulation function.** Given an ordering \\\pi\\, the curve is

\\ S\_\pi(k) \\=\\ \left\| \\ \bigcup\_{t=1}^{k} \mathrm{sp}(\pi_t) \\
\right\|, \qquad k = 1, \dots, n . \\

\\S\_\pi\\ is a step function: it rises by the count of species at
\\\pi_k\\ that were absent from the first \\k-1\\ sites, and never
falls. It reaches the total richness \\S\_\pi(n) = \|\bigcup_i
\mathrm{sp}(i)\|\\ regardless of the ordering, so methods differ only in
the *shape* of the approach, not the endpoint.

**Seed band.** Run the method from seeds \\s_1, \dots, s_B\\ (sampled
with replacement from the eligible sites), giving curves \\S^{(1)},
\dots, S^{(B)}\\ stacked as the \\B \times n\\ matrix returned in
`$curves`. The summary at step \\k\\ is the across-seed mean and the
empirical quantiles

\\ \bar S(k) = \frac{1}{B} \sum\_{b=1}^{B} S^{(b)}(k), \qquad \hat
q\_\alpha(k) = \text{quantile}\_\alpha\\\left( S^{(1)}(k), \dots,
S^{(B)}(k) \right). \\

The \\2.5\\\\ and \\97.5\\\\ quantiles give the default band. The
interval is not a parametric formula; it is the sampling distribution of
richness-at-effort induced by varying the starting point.

------------------------------------------------------------------------

## Expansion methods

Each method is a rule for choosing \\\pi\_{k+1}\\ given the sites
already visited. The seven rules fall into three families: walks that
chain through the cloud (`knn`, `kncn`, `gaussian`), orderings by
distance from a fixed seed (`radius`, `cone`), and geography-free
baselines (`random`, `collector`).

### k-nearest neighbour (`knn`)

From the current site \\\pi_k\\, move to the closest unvisited site:

\\ \pi\_{k+1} \\=\\ \arg\min\_{j \\\notin\\ \\\pi_1,\dots,\pi_k\\}
d(\pi_k, j). \\

The reference point is the *current* site, so the ordering is a
connected walk that follows local density. This is the default and the
curve the other methods are usually compared against.

### k-nearest centroid neighbour (`kncn`)

Track the centroid of the visited set and move to the unvisited site
closest to it:

\\ \bar c_k = \frac{1}{k} \sum\_{t=1}^{k} c\_{\pi_t}, \qquad \pi\_{k+1}
= \arg\min\_{j \\\notin\\ \\\pi_1,\dots,\pi_k\\} \\ c_j - \bar c_k \\ .
\\

Because the reference point is the centroid rather than the last site,
the visited set grows as a compact blob instead of a thread. `kncn`
resists the long jumps a `knn` walk makes when it exhausts a cluster, so
its early curve is smoother.

### Expanding radius (`radius`)

Sort every site by its distance from the seed and accumulate in that
order:

\\ \pi \\=\\ \text{argsort}\_j \\ d(s, j), \qquad s = \text{seed}. \\

The reference point is fixed at the seed for the whole curve, so the
ordering sweeps out a growing disc centred on the seed. This is the
cleanest “survey spreading outward from a point” interpretation. The
related
[`spaccWavefront()`](https://gillescolling.com/spacc/reference/spaccWavefront.md)
function parameterises the same idea by radius instead of by site count,
reporting richness as a function of the disc radius rather than the
number of sites included.

### Gaussian-weighted walk (`gaussian`)

A soft version of `knn`. From the current site, draw the next site at
random with probability proportional to a Gaussian kernel of distance:

\\ \Pr(\pi\_{k+1} = j) \\\propto\\ \exp\\\left( -\frac{d(\pi_k,
j)^2}{2\sigma^2} \right), \qquad j \notin \\\pi_1,\dots,\pi_k\\. \\

Nearby sites are favoured but not guaranteed, so the walk explores a
neighbourhood rather than always taking the single closest site. The
bandwidth \\\sigma\\ sets how sharp the preference is: small \\\sigma\\
approaches `knn`, large \\\sigma\\ approaches `random`. By default
\\\sigma\\ is the median of the non-zero pairwise distances.

### Directional cone (`cone`)

Pick a random direction \\\theta\\ for the seed. Sites whose bearing
from the seed falls within the half-width `cone_width` of \\\theta\\ are
accumulated first, in order of distance; the remaining sites follow
afterwards, also by distance:

\\ \text{in-cone}(j) \iff \bigl\| \angle(c_j - c_s) - \theta \bigr\| \le
\texttt{cone\\width}. \\

The cone models a survey that advances along a bearing, such as a
transect up a valley or along a coastline. The default half-width is
\\\pi/4\\ (a 90-degree wedge).

### Random order (`random`)

A uniform random permutation of the sites, independent of geography.
This is the classical accumulation curve and the null model the spatial
methods are measured against.

### Collector (`collector`)

The sites in the order they appear in the data, with no randomisation.
It produces a single curve and reproduces the “collector’s curve” of the
order in which a survey was actually conducted.

### Comparing the orderings

The same points, ordered by four rules, make the differences visible.
Colour encodes visiting order from first (dark) to last (light).

``` r

pal <- function(o) grDevices::hcl.colors(length(o), "Greens", rev = TRUE)[order(o)]
radius_order <- function(coords, seed)
  order(sqrt((coords$x - coords$x[seed])^2 + (coords$y - coords$y[seed])^2))

orders <- list(
  kNN       = knn_order(pts, 1),
  radius    = radius_order(pts, 1),
  random    = sample(nrow(pts)),
  collector = seq_len(nrow(pts))
)
op <- par(mfrow = c(2, 2), mar = c(2, 2, 2, 1))
for (nm in names(orders)) {
  rank_k <- match(seq_len(nrow(pts)), orders[[nm]])
  plot(pts$x, pts$y, pch = 19, cex = 1.6, col = pal(rank_k),
       xlab = "", ylab = "", main = nm, axes = FALSE); box()
}
```

![Four panels showing the same twenty points coloured by visiting order
under kNN, expanding radius, random, and collector rules. The kNN and
radius panels show smooth spatial gradients of colour; the random panel
shows no spatial pattern; the collector panel follows data
order.](theory_files/figure-html/method-schematic-1.svg)

``` r

par(op)
```

The `knn` and `radius` panels show a clear gradient: nearby sites are
visited at similar times. The `random` panel has no spatial pattern, and
`collector` follows whatever order the rows happen to be in.

------------------------------------------------------------------------

## Distance metrics

The walks depend on a distance \\d(i, j)\\ between sites, set by the
`distance` argument.

**Euclidean.** For projected or arbitrary planar coordinates,

\\ d(i,j) = \sqrt{(x_i - x_j)^2 + (y_i - y_j)^2}. \\

**Haversine.** For geographic coordinates given as longitude (`x`) and
latitude (`y`) in degrees, the great-circle distance in kilometres,

\\ d(i,j) = 2R \\\arcsin\\\sqrt{ \sin^2\\\tfrac{\Delta\varphi}{2} +
\cos\varphi_i \cos\varphi_j \\ \sin^2\\\tfrac{\Delta\lambda}{2} }, \\

with \\\varphi\\ latitude, \\\lambda\\ longitude, and \\R = 6371\\ km.
Haversine needs no projection step, so longitude-latitude data can be
used directly.

**Spatiotemporal.** When sites are sampled across time, supply a `time`
vector and the distance becomes a weighted sum of a spatial and a
temporal term:

\\ d(i,j) = w\_{\text{space}} \\ d\_{\text{space}}(i,j) +
w\_{\text{time}} \\ \|t_i - t_j\|. \\

The weights `w_space` and `w_time` trade geographic against temporal
proximity. A composite distance is not a metric a spatial tree can
index, so this mode always uses the exact backend and is available for
the methods that accept a distance matrix (`knn`, `radius`, `gaussian`).

------------------------------------------------------------------------

## From theory to implementation

The mathematical objects map onto arguments of
[`spacc()`](https://gillescolling.com/spacc/reference/spacc.md)
directly.

| Concept | Argument | Effect |
|----|----|----|
| Ordering rule | `method` | Selects one of the seven expansion rules |
| Distance \\d\\ | `distance` | `"euclidean"` or `"haversine"` |
| Number of seeds \\B\\ | `n_seeds` | Curves in the seed band |
| Bandwidth \\\sigma\\ | `sigma` | Gaussian kernel width (default: median distance) |
| Cone half-width | `cone_width` | Angular wedge for `cone` (default \\\pi/4\\) |
| Temporal axis | `time`, `w_space`, `w_time` | Switches to composite distance |
| Fixed ordering | `order` | Supplies \\\pi\\ directly, bypassing `method` |
| Species split | `groups` | One curve per group, same site ordering |
| Spatial support | `support`, `include_halo` | Seeds drawn from core sites only |
| Query backend | `backend` | `"auto"`, `"exact"`, or `"kdtree"` |

A typical call states the rule, the metric, and the number of seeds:

``` r

coords  <- data.frame(x = runif(120), y = runif(120))
species <- matrix(rbinom(120 * 40, 1, 0.25), nrow = 120)

sac <- spacc(species, coords, method = "knn", n_seeds = 100,
             progress = FALSE, seed = 42)
sac
#> spacc: 120 sites, 40 species, 100 seeds (knn)
```

Supplying `order` brings an externally computed ordering, for example
one from a sampling design or a different package, into the same
band-and-plot machinery. Each ordering is treated like a seed, so a
matrix of orderings yields a band.

``` r

ord_mat <- t(replicate(20, sample(nrow(species))))
sac_user <- spacc(species, coords, order = ord_mat, progress = FALSE)
sac_user$n_seeds
#> [1] 20
```

------------------------------------------------------------------------

## The two-tier backend

The cost of an accumulation walk is dominated by nearest-neighbour
queries. Two backends answer them, and `backend = "auto"` chooses
between them by site count.

**Exact.** Precompute the full \\n \times n\\ distance matrix once, then
answer each query by scanning a row. The matrix costs \\O(n^2)\\ memory
and the scan costs \\O(n)\\ per step, so a full walk is \\O(n^2)\\. For
up to a few hundred sites this is the fastest option and the matrix fits
comfortably in memory.

**Spatial tree.** Build a spatial index and query it in roughly \\O(\log
n)\\, avoiding the quadratic matrix entirely. For Euclidean distances
`spacc` uses a k-d tree (via nanoflann); for haversine distances it uses
a ball tree, whose spherical-cap geometry suits great-circle distance.
The build is \\O(n \log n)\\ and the walk is about \\O(n \log n)\\,
which wins decisively once \\n\\ grows past a few hundred.

| Backend   | Distance  | Memory     | Per query             | Selected when           |
|-----------|-----------|------------|-----------------------|-------------------------|
| Exact     | any       | \\O(n^2)\\ | \\O(n)\\              | \\n \le 500\\ (auto)    |
| k-d tree  | Euclidean | \\O(n)\\   | \\\approx O(\log n)\\ | \\n \> 500\\ (auto)     |
| Ball tree | Haversine | \\O(n)\\   | \\\approx O(\log n)\\ | \\n \> 500\\, haversine |

The `auto` rule switches to a tree above 500 sites. The two backends
compute the same ordering rule, so they agree on the curve up to ties;
only the speed differs.

``` r

e <- spacc(species, coords, method = "knn", backend = "exact",
           n_seeds = 30, progress = FALSE, seed = 7)
k <- spacc(species, coords, method = "knn", backend = "kdtree",
           n_seeds = 30, progress = FALSE, seed = 7)
c(exact_mean_end = mean(e$curves[, ncol(e$curves)]),
  kdtree_mean_end = mean(k$curves[, ncol(k$curves)]))
#>  exact_mean_end kdtree_mean_end 
#>              40              40
```

The `radius`, `gaussian`, and `cone` methods always use coordinates or a
distance matrix directly, and the spatiotemporal composite forces the
exact backend because a summed space-time distance has no tree to index.

------------------------------------------------------------------------

## Seeds and uncertainty

A single spatial curve is one realisation. The uncertainty comes from
the choice of starting point, so `spacc` repeats the expansion from
`n_seeds` seeds and reads the spread off the resulting band. Each seed
is an independent walk with no shared state, which makes the computation
embarrassingly parallel: the seeds are distributed across threads by
RcppParallel, and the per-step quantiles are taken after all walks
finish.

``` r

plot(sac)
```

![A spatial accumulation curve with a shaded confidence band. The mean
curve rises and saturates; the band is widest in the early, low-effort
region where the choice of starting site matters
most.](theory_files/figure-html/seed-band-1.svg)

The band is widest early and narrows as the curves converge on the
shared endpoint. Early on, the seed dominates: a walk that starts in a
rich patch and one that starts in a poor patch disagree most after a
handful of sites. By the time most of the map is covered, every walk has
seen nearly every species, so the curves meet. More seeds tighten the
*estimate* of the band but do not change its shape; the band reflects
real across-start variability, not Monte Carlo error that vanishes with
more replicates.

``` r

b50  <- spacc(species, coords, n_seeds = 50,  progress = FALSE, seed = 1)
b300 <- spacc(species, coords, n_seeds = 300, progress = FALSE, seed = 1)
mid <- round(ncol(b50$curves) / 4)
c(seeds_50  = diff(quantile(b50$curves[,  mid], c(.025, .975))),
  seeds_300 = diff(quantile(b300$curves[, mid], c(.025, .975))))
#>  seeds_50.97.5% seeds_300.97.5% 
#>               0               0
```

------------------------------------------------------------------------

## Relationship to random-order curves

The classical accumulation curve in
[`vegan::specaccum()`](https://vegandevs.github.io/vegan/reference/specaccum.html)
draws sites in random or collector order and averages over permutations.
`spacc` reproduces that curve with `method = "random"` or
`method = "collector"` and adds the spatial methods on top, so the two
can be read against each other. The
[`compare()`](https://gillescolling.com/spacc/reference/compare.md)
function tests the gap between any two curves by permutation, bootstrap,
or area-under-curve.

``` r

sp <- spacc(species, coords, method = "knn",    n_seeds = 100, progress = FALSE, seed = 3)
rd <- spacc(species, coords, method = "random", n_seeds = 100, progress = FALSE, seed = 3)
cmp <- compare(sp, rd)
plot(cmp)
```

![Two accumulation curves on the same axes: the random-order curve rises
faster and the spatial kNN curve lags behind it, with the shaded gap
between them representing spatial
turnover.](theory_files/figure-html/compare-random-1.svg)

When species are aggregated, the spatial curve lags the random curve:
neighbouring sites share species, so spreading outward discovers them
more slowly than mixing distant communities at every step. The size of
the gap is a measure of spatial turnover. A spatial curve that tracks
the random one signals a community with little spatial structure; a
curve that lags far behind signals strong distance decay in composition.
An existing
[`vegan::specaccum()`](https://vegandevs.github.io/vegan/reference/specaccum.html)
result can be imported with
[`as_spacc()`](https://gillescolling.com/spacc/reference/as_spacc.md) to
enter the same plotting and comparison machinery.

------------------------------------------------------------------------

## From curves to other diversity measures

The expansion rule fixes an ordering; what is counted at each step is
free to change. The downstream functions reuse the spatial walk and
replace the quantity accumulated, so the spatial logic in this vignette
carries through unchanged.

- [`spaccHill()`](https://gillescolling.com/spacc/reference/spaccHill.md)
  accumulates Hill numbers of order \\q = 0, 1, 2\\ (richness, the
  exponential of Shannon entropy, the inverse Simpson index), counting
  effective species instead of raw species.
- [`spaccBeta()`](https://gillescolling.com/spacc/reference/spaccBeta.md)
  accumulates beta diversity and partitions it into turnover and
  nestedness components.
- [`spaccCoverage()`](https://gillescolling.com/spacc/reference/spaccCoverage.md)
  tracks sample coverage, so curves can be standardised to a common
  completeness rather than a common site count.
- [`spaccPhylo()`](https://gillescolling.com/spacc/reference/spaccPhylo.md)
  and
  [`spaccFunc()`](https://gillescolling.com/spacc/reference/spaccFunc.md)
  accumulate phylogenetic and functional diversity along the same
  ordering.

``` r

hill <- spaccHill(species, coords, q = c(0, 1, 2), n_seeds = 50, progress = FALSE, seed = 9)
plot(hill)
```

![Hill-number accumulation curves for orders q=0, q=1, and q=2 along the
spatial ordering, with the q=0 richness curve highest and the q=2
inverse-Simpson curve
lowest.](theory_files/figure-html/downstream-taste-1.svg)

Each of these has its own vignette; the point here is that they share
the accumulation core. Choosing a `method` and `distance` configures the
walk once, and every diversity measure inherits it.

------------------------------------------------------------------------

## Design notes

**Why several expansion methods?** No single ordering rule is correct
for every question. A growing disc (`radius`) matches a survey expanding
from a point; a chained walk (`knn`) matches a surveyor moving to the
nearest accessible site; a directional cone matches a transect along a
bearing. Offering the rules as one argument lets the analysis state its
sampling model explicitly rather than defaulting to the random null.

**Why percentile bands instead of a parametric interval?** The
uncertainty in a spatial curve is dominated by where the survey starts,
which has no closed-form distribution. Resampling the seed and reading
off empirical quantiles measures that variability directly, without
assuming the curve follows a particular family.

**Why reduce abundance to presence-absence?** A species accumulation
curve counts first appearances. Whether a species has one individual or
a thousand at a site does not change whether it is new to the running
total, so the accumulation core works on presence-absence.
Abundance-aware questions are answered by the Hill and coverage
extensions, which keep the counts.

**Why a hard backend switch at 500 sites?** Below a few hundred sites
the quadratic distance matrix is both small and the fastest option;
above that the memory and the per-query scan grow faster than a tree’s
logarithmic lookup. The threshold is a default that `backend` overrides
when a particular dataset or benchmark calls for it.

------------------------------------------------------------------------

## References

Arrhenius, O. (1921). Species and area. *Journal of Ecology*, 9, 95-99.

Scheiner, S. M. (2003). Six types of species-area curves. *Global
Ecology and Biogeography*, 12, 441-447.

Chiarucci, A., Bacaro, G., & Scheiner, S. M. (2011). Old and new
challenges in using species diversity for assessing biodiversity.
*Philosophical Transactions of the Royal Society B*, 366, 2426-2437.

Gotelli, N. J., & Colwell, R. K. (2001). Quantifying biodiversity:
procedures and pitfalls in the measurement and comparison of species
richness. *Ecology Letters*, 4, 379-391.

Colwell, R. K., Chao, A., Gotelli, N. J., Lin, S.-Y., Mao, C. X.,
Chazdon, R. L., & Longino, J. T. (2012). Models and estimators linking
individual-based and sample-based rarefaction, extrapolation and
comparison of assemblages. *Journal of Plant Ecology*, 5, 3-21.

Chao, A., Gotelli, N. J., Hsieh, T. C., Sander, E. L., Ma, K. H.,
Colwell, R. K., & Ellison, A. M. (2014). Rarefaction and extrapolation
with Hill numbers: a framework for sampling and estimation in species
diversity studies. *Ecological Monographs*, 84, 45-67.
[doi:10.1890/13-0133.1](https://doi.org/10.1890/13-0133.1)

Ugland, K. I., Gray, J. S., & Ellingsen, K. E. (2003). The
species-accumulation curve and estimation of species richness. *Journal
of Animal Ecology*, 72, 888-897.
[doi:10.1046/j.1365-2656.2003.00748.x](https://doi.org/10.1046/j.1365-2656.2003.00748.x)

Nekola, J. C., & White, P. S. (1999). The distance decay of similarity
in biogeography and ecology. *Journal of Biogeography*, 26, 867-878.

Shigesada, N., & Kawasaki, K. (1997). *Biological Invasions: Theory and
Practice*. Oxford University Press.

Bentley, J. L. (1975). Multidimensional binary search trees used for
associative searching. *Communications of the ACM*, 18, 509-517.

Omohundro, S. M. (1989). *Five balltree construction algorithms*.
Technical Report 89-063, International Computer Science Institute,
Berkeley.

Blanco, J. L., & Rai, P. K. (2014). nanoflann: a C++ header-only library
for nearest neighbor (kNN) search with k-d trees.
<https://github.com/jlblancoc/nanoflann>

## See also

- [`vignette("quickstart")`](https://gillescolling.com/spacc/articles/quickstart.md)
  – the basic workflow end to end

- [`vignette("diversity")`](https://gillescolling.com/spacc/articles/diversity.md)
  – Hill numbers, beta diversity, coverage

- [`vignette("extrapolation")`](https://gillescolling.com/spacc/articles/extrapolation.md)
  – asymptotic richness models

- [`vignette("spatial-analysis")`](https://gillescolling.com/spacc/articles/spatial-analysis.md)
  – distance decay, endemism, zeta diversity

## Session info

``` r

sessionInfo()
#> R version 4.6.0 (2026-04-24 ucrt)
#> Platform: x86_64-w64-mingw32/x64
#> Running under: Windows 11 x64 (build 26200)
#> 
#> Matrix products: default
#>   LAPACK version 3.12.1
#> 
#> locale:
#> [1] LC_COLLATE=English_United States.utf8 
#> [2] LC_CTYPE=English_United States.utf8   
#> [3] LC_MONETARY=English_United States.utf8
#> [4] LC_NUMERIC=C                          
#> [5] LC_TIME=English_United States.utf8    
#> 
#> time zone: Europe/Luxembourg
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] spacc_0.9.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] gtable_0.3.6          jsonlite_2.0.0        dplyr_1.2.1          
#>  [4] compiler_4.6.0        tidyselect_1.2.1      Rcpp_1.1.1-1.1       
#>  [7] parallel_4.6.0        jquerylib_0.1.4       systemfonts_1.3.2    
#> [10] scales_1.4.0          textshaping_1.0.5     yaml_2.3.12          
#> [13] fastmap_1.2.0         ggplot2_4.0.3         R6_2.6.1             
#> [16] labeling_0.4.3        generics_0.1.4        knitr_1.51           
#> [19] htmlwidgets_1.6.4     tibble_3.3.1          desc_1.4.3           
#> [22] svglite_2.2.2         bslib_0.11.0          pillar_1.11.1        
#> [25] RColorBrewer_1.1-3    rlang_1.2.0           cachem_1.1.0         
#> [28] xfun_0.57             fs_2.1.0              sass_0.4.10          
#> [31] S7_0.2.2              RcppParallel_5.1.11-2 otel_0.2.0           
#> [34] cli_3.6.6             withr_3.0.2           pkgdown_2.2.0        
#> [37] magrittr_2.0.5        digest_0.6.39         grid_4.6.0           
#> [40] lifecycle_1.0.5       vctrs_0.7.3           evaluate_1.0.5       
#> [43] glue_1.8.1            farver_2.1.2          rmarkdown_2.31       
#> [46] pkgconfig_2.0.3       tools_4.6.0           htmltools_0.5.9
```
