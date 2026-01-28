# spacc Feature Plan

Planned features for spacc v0.2.0, based on a review of the spatial
accumulation / SAR / diversity literature. Each feature extends the
existing C++ backend and S3 class system.

------------------------------------------------------------------------

## 1. Faith’s Phylogenetic Diversity Accumulation

**Status**: Stub exists in
[`spaccPhylo()`](https://gillescolling.com/spacc/reference/spaccPhylo.md)
(“not yet implemented”)

**What**: Accumulate Faith’s PD (sum of branch lengths connecting
species in the sample) as sites are added spatially. Currently
[`spaccPhylo()`](https://gillescolling.com/spacc/reference/spaccPhylo.md)
supports MPD and MNTD but not PD.

**Implementation**: - Accept a `phylo` object (ape), not just a distance
matrix - At each accumulation step, compute the total branch length of
the minimum spanning subtree connecting all observed species - C++
function `cpp_phylo_pd_knn_parallel` operating on an edge matrix
representation - Returns PD curve alongside existing MPD/MNTD in the
`spacc_phylo` object

**References**: - Faith, D.P. (1992). Conservation evaluation and
phylogenetic diversity. *Biological Conservation*, 61, 1–10. - Webb,
C.O. (2000). Exploring the phylogenetic structure of ecological
communities. *The American Naturalist*, 156, 145–155.

**Depends on**: ape (already in Suggests)

------------------------------------------------------------------------

## 2. Coverage-Based Extrapolation

**Status**:
[`spaccCoverage()`](https://gillescolling.com/spacc/reference/spaccCoverage.md)
and
[`interpolateCoverage()`](https://gillescolling.com/spacc/reference/interpolateCoverage.md)
exist for interpolation only

**What**: Extend coverage-based rarefaction to extrapolation beyond the
observed sample, following the Chao & Jost (2012) and Chao et al. (2014)
framework. Allow users to predict diversity at target coverage levels
that exceed the empirical maximum.

**Implementation**: - New function
`extrapolateCoverage(object, target_coverage, q = 0)` operating on
`spacc_coverage` objects - For q = 0: Chao1/Chao2 asymptotic estimator +
coverage-based extrapolation formula - For q = 1, 2: Closed-form Hill
number extrapolation (Chao et al. 2014, eq. 5–6) - Seamless
interpolation/extrapolation curve with a visual marker at the reference
sample - New class `spacc_coverage_ext` extending `spacc_coverage`

**References**: - Chao, A. & Jost, L. (2012). Coverage-based rarefaction
and extrapolation: standardizing samples by completeness rather than
size. *Ecology*, 93, 2533–2547. - Chao, A., Gotelli, N.J., Hsieh, T.C.,
et al. (2014). Rarefaction and extrapolation with Hill numbers: a
framework for sampling and estimation in species diversity studies.
*Ecological Monographs*, 84, 45–67. - Hsieh, T.C., Ma, K.H. & Chao, A.
(2016). iNEXT: an R package for rarefaction and extrapolation of species
diversity (Hill numbers). *Methods in Ecology and Evolution*, 7,
1451–1456.

**CRAN packages**: [iNEXT](https://cran.r-project.org/package=iNEXT),
[SpadeR](https://cran.r-project.org/package=SpadeR)

------------------------------------------------------------------------

## 3. Diversity-Area Relationship (DAR)

**Status**: Not implemented

**What**: Extend the classic SAR (species vs. area) to the DAR
framework, plotting Hill numbers of any order q against cumulative area
rather than cumulative site count. This captures how *effective*
diversity — not just richness — scales with area.

**Implementation**: - New function `dar(object, coords, q = c(0, 1, 2))`
that converts the site-based accumulation from
[`spaccHill()`](https://gillescolling.com/spacc/reference/spaccHill.md)
into an area-based relationship - Area estimation via Voronoi
tessellation (sf) or convex hull of accumulated sites at each step -
Optionally accept a cell_size argument for grid-based area estimation -
New class `spacc_dar` with plot method showing diversity-area curves for
each q - Log-log transformation option for power-law fitting

**References**: - Ma, Z.S. (2018). DAR (diversity-area relationship):
extending classic SAR (species-area relationship) for biodiversity and
biogeography analyses. *Ecology and Evolution*, 8, 10023–10038. -
Arrhenius, O. (1921). Species and area. *Journal of Ecology*, 9, 95–99.

**CRAN packages**: [vegan](https://cran.r-project.org/package=vegan)
(SAR fitting), [mobsim](https://cran.r-project.org/package=mobsim)
(divar function)

------------------------------------------------------------------------

## 4. Extreme Value Theory SAR Model

**Status**: Not implemented

**What**: Add an EVT-based model to
[`extrapolate()`](https://gillescolling.com/spacc/reference/extrapolate.md)
that captures the characteristic three-phase SAR pattern (rapid rise,
plateau, second rise) observed in log-log space. Standard asymptotic
models (Michaelis-Menten, etc.) cannot capture this triphasic shape.

**Implementation**: - Add `"evt"` as a new model option in
[`extrapolate()`](https://gillescolling.com/spacc/reference/extrapolate.md) -
The EVT model treats the SAR as a mixture of species-specific
minimum-distance distributions - Fit via maximum likelihood using a
Generalized Extreme Value (GEV) mixture - Parameters: location (mu),
scale (sigma), shape (xi) per species group, or a simplified
two-component mixture - Return `spacc_fit` object as usual, with
additional diagnostics for phase transitions - Log-log residual plot to
visualize the three phases

**References**: - Borda-de-Agua, L., Whittaker, R.J., Cardoso, P., et
al. (2025). Modelling the species-area relationship using extreme value
theory. *Nature Communications*, 16, 4045. - Scheather, S.J. &
Hettmansperger, T.P. (2009). Extreme value theory. In *The Species-Area
Relationship* (eds Storch, D., Marquet, P.A. & Brown, J.H.), Cambridge
University Press.

**CRAN packages**: [evd](https://cran.r-project.org/package=evd),
[extRemes](https://cran.r-project.org/package=extRemes)

------------------------------------------------------------------------

## 5. Functional Beta Diversity Accumulation

**Status**: Not implemented (taxonomic beta exists via
[`spaccBeta()`](https://gillescolling.com/spacc/reference/spaccBeta.md))

**What**: Extend beta diversity partitioning from taxonomic to
functional and phylogenetic dimensions. Decompose
functional/phylogenetic beta diversity into turnover and nestedness
components as sites accumulate spatially.

**Implementation**: - Extend
[`spaccBeta()`](https://gillescolling.com/spacc/reference/spaccBeta.md)
with a `type` parameter: `"taxonomic"` (default, current behavior),
`"functional"`, `"phylogenetic"` - Functional beta: Requires traits
matrix, computes functional Sorensen/Jaccard based on trait space
overlap (convex hull or kernel density) - Phylogenetic beta: Requires
phylo object, uses PhyloSor or UniFrac-based decomposition -
Turnover/nestedness partitioning following Baselga (2010, 2012) extended
to multiple facets - C++ backend `cpp_beta_func_knn_parallel` and
`cpp_beta_phylo_knn_parallel` - Returns `spacc_beta` object with type
attribute; plot method adjusts labels accordingly

**References**: - Baselga, A. (2010). Partitioning the turnover and
nestedness components of beta diversity. *Global Ecology and
Biogeography*, 19, 134–143. - Baselga, A. (2012). The relationship
between species replacement, dissimilarity derived from nestedness, and
nestedness. *Global Ecology and Biogeography*, 21, 1223–1232. - Chao,
A., Chiu, C.H., Vill'eger, S., et al. (2023). Rarefaction and
extrapolation with beta diversity under a framework of Hill numbers: the
iNEXT.beta3D standardization. *Ecological Monographs*, 93, e1588. -
Cardoso, P., Rigal, F. & Carvalho, J.C. (2015). BAT – Biodiversity
Assessment Tools, an R package for the measurement and estimation of
alpha and beta taxon, phylogenetic and functional diversity. *Methods in
Ecology and Evolution*, 6, 232–236.

**CRAN packages**:
[betapart](https://cran.r-project.org/package=betapart) (already in
Suggests),
[iNEXT.beta3D](https://cran.r-project.org/package=iNEXT.beta3D),
[BAT](https://cran.r-project.org/package=BAT)

------------------------------------------------------------------------

## 6. Sampling Effort SAR (SESARS)

**Status**: Not implemented

**What**: Model the joint effect of sampling effort and area on species
richness. Standard SARs assume complete sampling within each area unit;
SESARS corrects for unequal survey intensity across sites — common in
atlas data, citizen science datasets, and multi-year monitoring.

**Implementation**: - New function `sesars(object, effort)` where
`effort` is a numeric vector (e.g., sampling hours, number of visits,
trap-nights per site) - Fits a joint model: S ~ f(A, E) where A =
cumulative area and E = cumulative effort - Model options:
multiplicative power-law `S = c * A^z * E^w`, or additive with separate
area and effort terms - Predict richness at standardized effort levels -
New class `spacc_sesars` with print/summary/plot methods - Plot: 3D
surface or 2D slices at fixed effort levels

**References**: - Ribas, C.R., Sobrinho, T.G., Schoereder, J.H., et
al. (2012). How large is large enough for insects? Forest fragment size
and sampling effort for insects. *Environmental Entomology*, 41,
965–972. - Sousa-Baena, M.S., Garcia, L.C. & Peterson, A.T. (2014).
Completeness of digital accessible knowledge of the plants of Brazil.
*PLoS ONE*, 9, e84068. - Dennstadt, F., Horak, J. & Martin, M.D. (2019).
Predictive sampling effort and species-area relationship models for
estimating richness in fragmented landscapes. *Diversity and
Distributions*, 26, 1112–1123.

**CRAN packages**: [sars](https://cran.r-project.org/package=sars) (SAR
model fitting)

------------------------------------------------------------------------

## 7. Species-Fragmented Area Relationship (SFAR)

**Status**: Not implemented

**What**: Extend the power-law SAR to separately quantify the effects of
habitat loss (total area reduction) and fragmentation (splitting into
patches). The SFAR adds a fragmentation term to the classic SAR,
allowing users to disentangle these two drivers of species loss.

**Implementation**: - New function `sfar(object, patches)` where
`patches` is an sf polygon layer or a grouping vector assigning sites to
habitat fragments - Model: `S = c * A^z * n^(-f)` where A = total area,
n = number of fragments, f = fragmentation exponent - Alternatively:
Hanski et al. metapopulation-informed extension with patch isolation -
Fit via nls or maximum likelihood - New class `spacc_sfar` with
print/summary/plot methods - Plot: Observed vs predicted richness, with
separate area and fragmentation effect curves

**References**: - Hanski, I., Zurita, G.A., Bellocq, M.I. & Rybicki, J.
(2013). Species-fragmented area relationship. *Proceedings of the
National Academy of Sciences*, 110, 12715–12720. - Rybicki, J. & Hanski,
I. (2013). Species-area relationships and extinctions caused by habitat
loss and fragmentation. *Ecology Letters*, 16, 27–38.

**CRAN packages**: [sars](https://cran.r-project.org/package=sars)

------------------------------------------------------------------------

## 8. Endemism-Area Relationship

**Status**: Not implemented

**What**: Compute the number of endemic species (species found *only*
within the accumulated area) as a function of cumulative area or site
count. Complements the standard SAR by tracking how many species are
unique to each spatial extent — critical for conservation
prioritization.

**Implementation**: - New function `spaccEndemism(x, coords, ...)` or
add `metric = "endemism"` option to
[`spaccHill()`](https://gillescolling.com/spacc/reference/spaccHill.md) -
At each accumulation step, count species present in the accumulated set
but absent from all remaining (unvisited) sites - This is the
complement: endemic_k = species in sites 1..k that are absent from sites
k+1..n - C++ function `cpp_knn_endemism_parallel` tracking both
cumulative and exclusive species counts - New class `spacc_endemism` or
extend `spacc` with an endemism column - Plot: Endemism curve overlaid
on or alongside the standard SAC

**References**: - Hobohm, C. (ed.) (2014). *Endemism in Vascular
Plants*. Springer. - Kier, G., Kreft, H., Lee, T.M., et al. (2009). A
global assessment of endemism and species richness across island and
mainland regions. *Proceedings of the National Academy of Sciences*,
106, 9322–9327. - May, F., Gerstner, K., McGlinn, D.J., et al. (2018).
mobsim: An R package for the simulation and measurement of biodiversity
across spatial scales. *Methods in Ecology and Evolution*, 9, 1401–1408.

**CRAN packages**: [mobsim](https://cran.r-project.org/package=mobsim)
(endemic species in divar)

------------------------------------------------------------------------

## Implementation Priority

| \# | Feature | Complexity | Builds on |
|----|----|----|----|
| 1 | Faith’s PD accumulation | Low | [`spaccPhylo()`](https://gillescolling.com/spacc/reference/spaccPhylo.md), ape |
| 2 | Coverage-based extrapolation | Medium | [`spaccCoverage()`](https://gillescolling.com/spacc/reference/spaccCoverage.md), [`interpolateCoverage()`](https://gillescolling.com/spacc/reference/interpolateCoverage.md) |
| 3 | Diversity-Area Relationship | Medium | [`spaccHill()`](https://gillescolling.com/spacc/reference/spaccHill.md), sf |
| 4 | EVT SAR model | Medium | [`extrapolate()`](https://gillescolling.com/spacc/reference/extrapolate.md) |
| 5 | Functional/phylogenetic beta | High | [`spaccBeta()`](https://gillescolling.com/spacc/reference/spaccBeta.md), [`spaccFunc()`](https://gillescolling.com/spacc/reference/spaccFunc.md), [`spaccPhylo()`](https://gillescolling.com/spacc/reference/spaccPhylo.md) |
| 6 | SESARS | Medium | [`extrapolate()`](https://gillescolling.com/spacc/reference/extrapolate.md), new effort data |
| 7 | SFAR | Medium | New, sf patches |
| 8 | Endemism-area | Low | [`spacc()`](https://gillescolling.com/spacc/reference/spacc.md) C++ backend |
