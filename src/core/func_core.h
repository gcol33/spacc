// core/func_core.h
// Functional diversity of a species set from a species-by-traits matrix:
// functional dispersion (FDis), range-product functional richness and Rao's
// quadratic entropy on Euclidean trait distances, and the aggregated
// functional diversity K with its four facets.
// Pure C++: traits are accessed through traits(species, trait).

#ifndef SPACC_CORE_FUNC_H
#define SPACC_CORE_FUNC_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
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

// Facets of the aggregated functional diversity K (Wojcik et al. 2024).
struct KFacets {
  double richness, biomass_evenness, trait_evenness, dispersion, K;
};

// K for the present species, one trait axis at a time and combined as the
// geometric mean over axes; facets are reported as their geometric means over
// axes, so K is the geometric mean of the reported facets. Per axis, with the
// n species sorted by trait value T and relative abundances b:
//   FRic = (max T - min T) / L,          L = axis_length[t]
//   BE   = 1 - sqrt(sum (b_i - 1/n)^2 / ((n-1)/n))
//   TE   = 1 - sqrt(sum (gap_i / range - 1/(n-1))^2 / ((n-2)/(n-1)))
//   Dis  = min(V_obs / V_uni, 1),  V_obs = sum b_i (T_i - mean_b T)^2,
//          V_uni = range^2 / 12
// TE needs three species; with two, K is the cube root of the other three
// facets and TE is NaN. With fewer than two species, or a zero trait range on
// any axis, K and its facets are 0 (TE NaN).
template<class TraitFn, class AbundFn>
KFacets functional_facets_K(const TraitFn& traits, int n_traits,
                            const std::vector<int>& present,
                            const AbundFn& abundance,
                            const std::vector<double>& axis_length) {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const int n = static_cast<int>(present.size());
  double total = 0.0;
  for (int sp : present) total += abundance(sp);
  if (n < 2 || total <= 0.0) return KFacets{0.0, 0.0, nan, 0.0, 0.0};

  double be_sum = 0.0;
  for (int sp : present) {
    const double d = abundance(sp) / total - 1.0 / n;
    be_sum += d * d;
  }
  const double be = std::max(0.0, 1.0 - std::sqrt(be_sum / ((n - 1.0) / n)));

  double log_fric = 0.0, log_te = 0.0, log_dis = 0.0;
  std::vector<std::pair<double, double>> axis(n);   // (trait value, relative abundance)
  for (int t = 0; t < n_traits; t++) {
    for (int i = 0; i < n; i++) {
      axis[i] = {traits(present[i], t), abundance(present[i]) / total};
    }
    std::sort(axis.begin(), axis.end());
    const double range = axis[n - 1].first - axis[0].first;
    if (range <= 0.0) return KFacets{0.0, be, nan, 0.0, 0.0};

    double mean = 0.0;
    for (const auto& a : axis) mean += a.second * a.first;
    double v_obs = 0.0;
    for (const auto& a : axis) v_obs += a.second * (a.first - mean) * (a.first - mean);
    const double dis = std::min(v_obs / (range * range / 12.0), 1.0);

    log_fric += std::log(range / axis_length[t]);
    log_dis += std::log(dis);
    if (n >= 3) {
      double te_sum = 0.0;
      for (int i = 0; i + 1 < n; i++) {
        const double d = (axis[i + 1].first - axis[i].first) / range - 1.0 / (n - 1);
        te_sum += d * d;
      }
      log_te += std::log(std::max(0.0, 1.0 - std::sqrt(te_sum / ((n - 2.0) / (n - 1.0)))));
    }
  }

  KFacets f;
  f.richness = std::exp(log_fric / n_traits);
  f.dispersion = std::exp(log_dis / n_traits);
  f.biomass_evenness = be;
  f.trait_evenness = n >= 3 ? std::exp(log_te / n_traits) : nan;
  f.K = n >= 3
    ? std::pow(f.richness * be * f.trait_evenness * f.dispersion, 0.25)
    : std::cbrt(f.richness * be * f.dispersion);
  return f;
}

} // namespace spacc

#endif // SPACC_CORE_FUNC_H
