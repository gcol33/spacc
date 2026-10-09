// core/func_core.h
// Functional diversity of a species set from a species-by-traits matrix:
// functional dispersion (FDis), range-product functional richness and Rao's
// quadratic entropy on Euclidean trait distances.
// Pure C++: traits are accessed through traits(species, trait).

#ifndef SPACC_CORE_FUNC_H
#define SPACC_CORE_FUNC_H

#include <cmath>
#include <limits>
#include <vector>

namespace spacc {

// Euclidean distance between two species in trait space.
template<class TraitFn>
double trait_distance(const TraitFn& traits, int n_traits, int i, int j) {
  double d2 = 0.0;
  for (int t = 0; t < n_traits; t++) {
    double diff = traits(i, t) - traits(j, t);
    d2 += diff * diff;
  }
  return std::sqrt(d2);
}

// Functional dispersion (Laliberte & Legendre 2010): abundance-weighted mean
// distance of the present species to their abundance-weighted centroid.
template<class TraitFn, class AbundFn>
double functional_dispersion(const TraitFn& traits, int n_traits,
                             const std::vector<int>& present,
                             const AbundFn& abundance) {
  double total = 0.0;
  for (int sp : present) total += abundance(sp);
  if (present.size() < 2 || total == 0.0) return 0.0;

  std::vector<double> centroid(n_traits, 0.0);
  for (int sp : present) {
    double w = abundance(sp) / total;
    for (int t = 0; t < n_traits; t++) centroid[t] += traits(sp, t) * w;
  }

  double sum_dist = 0.0;
  for (int sp : present) {
    double d2 = 0.0;
    for (int t = 0; t < n_traits; t++) {
      double diff = traits(sp, t) - centroid[t];
      d2 += diff * diff;
    }
    sum_dist += std::sqrt(d2) * (abundance(sp) / total);
  }
  return sum_dist;
}

// Range-product approximation of functional richness: the product of the
// per-trait ranges spanned by the present species. Zero unless there are more
// species than traits and every trait has positive range.
template<class TraitFn>
double functional_range_volume(const TraitFn& traits, int n_traits,
                               const std::vector<int>& present) {
  if (present.size() <= static_cast<size_t>(n_traits)) return 0.0;
  const double inf = std::numeric_limits<double>::infinity();
  double volume = 1.0;
  for (int t = 0; t < n_traits; t++) {
    double lo = inf, hi = -inf;
    for (int sp : present) {
      double v = traits(sp, t);
      if (v < lo) lo = v;
      if (v > hi) hi = v;
    }
    double range = hi - lo;
    if (range <= 0.0) return 0.0;
    volume *= range;
  }
  return volume;
}

// Rao's quadratic entropy on Euclidean trait distances over all n species.
template<class TraitFn, class AbundFn>
double rao_trait_entropy(const TraitFn& traits, int n_traits, int n_species,
                         const AbundFn& abundance) {
  double total = 0.0;
  for (int i = 0; i < n_species; i++) total += abundance(i);
  if (total <= 0.0) return 0.0;
  double q = 0.0;
  for (int i = 0; i < n_species; i++) {
    double ai = abundance(i);
    if (ai <= 0.0) continue;
    double pi = ai / total;
    for (int j = i + 1; j < n_species; j++) {
      double aj = abundance(j);
      if (aj <= 0.0) continue;
      q += pi * (aj / total) * trait_distance(traits, n_traits, i, j);
    }
  }
  return 2.0 * q;
}

} // namespace spacc

#endif // SPACC_CORE_FUNC_H
