# Story: spacc

## Claim
Spatial species accumulation curves: sites added in order of distance from many focal points, in C++, with the spread across focal points as uncertainty.

## Audience
Community ecologists and biogeographers who know vegan::specaccum() and iNEXT-style diversity curves.

## Relation
- to: vegan::specaccum() (classical random-order curve); iNEXT family (Hill numbers, coverage)
- phrase: the classical curve of `specaccum()` in vegan adds sites in random order
- same: method = "random" gives the random-order curve; as_spacc() imports specaccum objects
- differs: spatial ordering from continuous focal points; seven expansion orders
- home: vignettes/theory.Rmd (Relationship to random-order curves)

## Arcs

### A1 Where you start changes the curve
- question: does the spatial arrangement of my sites change the accumulation curve?
- answer: fixed-focus expansion from many focal points, compared with random order via compare()
- evidence: vignettes/quickstart.Rmd
- home: vignettes/quickstart.Rmd
- also: README Usage

### A2 Diversity beyond richness
- question: can I accumulate Hill numbers, beta, phylogenetic or functional diversity spatially?
- answer: spaccHill(), spaccBeta(), spaccCoverage(), spaccPhylo(), spaccFunc(), partitions
- home: vignettes/diversity.Rmd
- also: README

### A3 Extrapolation and total richness
- question: how many species are there beyond my sites?
- answer: extrapolate(), richness estimators
- home: vignettes/extrapolation.Rmd, vignettes/richness-estimation.Rmd
- also: README list

### A4 Scale
- question: does it handle many sites?
- answer: k-d tree (nanoflann) or ball tree above 500 sites, RcppParallel across focal points
- home: vignettes/theory.Rmd (The two-tier traversal backend)
- also: README list

## Non-goals
None stated in the docs.

## Documents
- DESCRIPTION: claim, methods with references
- README / landing: claim, relation, first curve, comparison, routing
- Get Started (quickstart): one path from matrix to curves, comparison, extrapolation
- articles: diversity, extrapolation, richness-estimation, rarefaction-standardization, spatial-analysis, community-assembly (how-to); theory (explanation)

## Open
- None.
