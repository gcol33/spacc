// core/accumulation_core.h
// Site-visitation orders and accumulation along an order.
// Pure C++ (no Rcpp): every routine takes accessors (lambdas or functors) so
// the same code runs on Rcpp views inside the package and on std::vector
// containers in the Catch2 tests.
//
// Accessor conventions
//   species(site, sp)   abundance (> 0 means present)
//   d(i, j)             distance between sites i and j
//   order(step)         0-based site visited at `step`
//   next(current, vis)  nearest unvisited site to site `current`
//   nearest(qx, qy, vis) nearest unvisited site to the point (qx, qy)

#ifndef SPACC_CORE_ACCUMULATION_H
#define SPACC_CORE_ACCUMULATION_H

#include "types.h"
#include "distance_core.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

namespace spacc {

// ---------------------------------------------------------------------------
// Nearest-unvisited searches (brute force)
// ---------------------------------------------------------------------------

template<class DistFn>
int nearest_unvisited(const DistFn& d, int current, int n_sites,
                      const VisitedFlags& visited) {
  double best = std::numeric_limits<double>::infinity();
  int nearest = -1;
  for (int j = 0; j < n_sites; j++) {
    if (!visited[j] && d(current, j) < best) {
      best = d(current, j);
      nearest = j;
    }
  }
  return nearest;
}

template<class PointDistFn>
int nearest_unvisited_point(double qx, double qy,
                            const std::vector<double>& x,
                            const std::vector<double>& y,
                            const VisitedFlags& visited,
                            const PointDistFn& dist) {
  double best = std::numeric_limits<double>::infinity();
  int nearest = -1;
  const int n = static_cast<int>(x.size());
  for (int j = 0; j < n; j++) {
    if (visited[j]) continue;
    double dj = dist(qx, qy, x[j], y[j]);
    if (dj < best) {
      best = dj;
      nearest = j;
    }
  }
  return nearest;
}

// ---------------------------------------------------------------------------
// Order generators. Each fills `order` (length n_sites) with 0-based sites.
// ---------------------------------------------------------------------------

// Nearest-neighbour walk: the reference point moves to the last visited site.
template<class NextFn>
void nn_walk_order(int n_sites, int seed, const NextFn& next,
                   std::vector<int>& order) {
  order.assign(n_sites, -1);
  VisitedFlags visited(n_sites, false);
  int current = seed;
  visited[current] = true;
  order[0] = current;
  for (int step = 1; step < n_sites; step++) {
    current = next(current, visited);
    visited[current] = true;
    order[step] = current;
  }
}

// Nearest-centroid walk: the reference point is the running centroid of the
// visited sites.
template<class NearestFn>
void kncn_order(const std::vector<double>& x, const std::vector<double>& y,
                int seed, const NearestFn& nearest, std::vector<int>& order) {
  const int n_sites = static_cast<int>(x.size());
  order.assign(n_sites, -1);
  VisitedFlags visited(n_sites, false);
  double sum_x = x[seed], sum_y = y[seed];
  int current = seed;
  visited[current] = true;
  order[0] = current;
  for (int step = 1; step < n_sites; step++) {
    double cx = sum_x / step, cy = sum_y / step;
    current = nearest(cx, cy, visited);
    visited[current] = true;
    sum_x += x[current];
    sum_y += y[current];
    order[step] = current;
  }
}

// Fixed-focus order: every site ranked by its distance from one focal point,
// ties broken by site index.
template<class DistToFocusFn>
void focal_order(int n_sites, const DistToFocusFn& dist_to_focus,
                 std::vector<int>& order) {
  std::vector<std::pair<double, int>> ranked(n_sites);
  for (int site = 0; site < n_sites; site++) {
    ranked[site] = std::make_pair(dist_to_focus(site), site);
  }
  std::stable_sort(ranked.begin(), ranked.end());
  order.resize(n_sites);
  for (int step = 0; step < n_sites; step++) order[step] = ranked[step].second;
}

// Gaussian-weighted walk: the next site is drawn among the unvisited ones with
// probability proportional to exp(-d^2 / (2 sigma^2)) from the current site.
// `unif()` returns a uniform deviate in [0, 1).
template<class DistFn, class UnifFn>
void gaussian_order(const DistFn& d, int n_sites, int seed, double sigma,
                    const UnifFn& unif, std::vector<int>& order) {
  order.assign(n_sites, -1);
  VisitedFlags visited(n_sites, false);
  const double two_sigma2 = 2.0 * sigma * sigma;
  std::vector<double> weight(n_sites);
  int current = seed;
  visited[current] = true;
  order[0] = current;
  for (int step = 1; step < n_sites; step++) {
    double total = 0.0;
    for (int j = 0; j < n_sites; j++) {
      double w = 0.0;
      if (!visited[j]) {
        double dj = d(current, j);
        w = std::exp(-(dj * dj) / two_sigma2);
      }
      weight[j] = w;
      total += w;
    }
    int chosen = -1;
    if (total > 0.0) {
      double u = unif() * total;
      double acc = 0.0;
      for (int j = 0; j < n_sites; j++) {
        if (visited[j]) continue;
        acc += weight[j];
        if (u < acc) { chosen = j; break; }
      }
    }
    if (chosen < 0) {
      // All weights underflowed to zero (or rounding at the top end): take the
      // nearest unvisited site, which is also the argmax of the weights.
      chosen = nearest_unvisited(d, current, n_sites, visited);
    }
    current = chosen;
    visited[current] = true;
    order[step] = current;
  }
}

// Directional cone: sites inside the angular window around `angle` (seen from
// the seed) are visited first by increasing distance, then the remaining sites
// by increasing distance.
inline void cone_order(const std::vector<double>& x, const std::vector<double>& y,
                       int seed, double angle, double width,
                       std::vector<int>& order) {
  const int n_sites = static_cast<int>(x.size());
  std::vector<std::pair<double, int>> in_cone, out_cone;
  in_cone.reserve(n_sites);
  out_cone.reserve(n_sites);
  for (int i = 0; i < n_sites; i++) {
    double dx = x[i] - x[seed], dy = y[i] - y[seed];
    double dist = std::sqrt(dx * dx + dy * dy);
    double diff = std::fabs(std::atan2(dy, dx) - angle);
    if (diff > PI) diff = 2.0 * PI - diff;
    if (diff <= width || i == seed) in_cone.emplace_back(dist, i);
    else out_cone.emplace_back(dist, i);
  }
  std::sort(in_cone.begin(), in_cone.end());
  std::sort(out_cone.begin(), out_cone.end());
  order.clear();
  order.reserve(n_sites);
  for (const auto& p : in_cone) order.push_back(p.second);
  for (const auto& p : out_cone) order.push_back(p.second);
}

// ---------------------------------------------------------------------------
// Accumulation along an order
// ---------------------------------------------------------------------------

// Add each visited site's row into `cumulative`, then call at_step(step, site).
template<class SpeciesFn, class OrderFn, class T, class StepFn>
void accumulate_along_order(const SpeciesFn& species, int n_sites, int n_species,
                            const OrderFn& order, std::vector<T>& cumulative,
                            const StepFn& at_step) {
  for (int step = 0; step < n_sites; step++) {
    const int site = order(step);
    for (int sp = 0; sp < n_species; sp++) {
      cumulative[sp] += static_cast<T>(species(site, sp));
    }
    at_step(step, site);
  }
}

// Cumulative species richness along an order; curve(step) receives the count.
template<class SpeciesFn, class OrderFn, class CurveFn>
void richness_along_order(const SpeciesFn& species, int n_sites, int n_species,
                          const OrderFn& order, const CurveFn& curve) {
  std::vector<bool> seen(n_species, false);
  int richness = 0;
  for (int step = 0; step < n_sites; step++) {
    const int site = order(step);
    for (int sp = 0; sp < n_species; sp++) {
      if (!seen[sp] && species(site, sp) > 0) {
        seen[sp] = true;
        richness++;
      }
    }
    curve(step, richness);
  }
}

// Species richness and site count within `radius` of `seed`.
template<class SpeciesFn, class DistFn>
std::pair<int, int> richness_within_radius(const SpeciesFn& species,
                                           const DistFn& d,
                                           int n_sites, int n_species,
                                           int seed, double radius) {
  std::vector<bool> seen(n_species, false);
  int richness = 0, included = 0;
  for (int i = 0; i < n_sites; i++) {
    if (d(seed, i) > radius) continue;
    included++;
    for (int sp = 0; sp < n_species; sp++) {
      if (!seen[sp] && species(i, sp) > 0) {
        seen[sp] = true;
        richness++;
      }
    }
  }
  return std::make_pair(richness, included);
}

} // namespace spacc

#endif // SPACC_CORE_ACCUMULATION_H
