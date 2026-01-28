// core/accumulation_core.h
// Pure C++ species accumulation algorithms for spacc
// No Rcpp dependencies - testable with Catch2

#ifndef SPACC_CORE_ACCUMULATION_H
#define SPACC_CORE_ACCUMULATION_H

#include "types.h"
#include "distance_core.h"
#include <algorithm>
#include <limits>

namespace spacc {

// Accumulate species from a site into the seen set
// species_pa: presence-absence row for site
// Returns number of new species added
inline int accumulate_species(
    const std::vector<int>& species_pa,
    SpeciesSet& seen) {

  int new_species = 0;
  for (size_t sp = 0; sp < species_pa.size(); sp++) {
    if (species_pa[sp] > 0) {
      auto result = seen.insert(static_cast<int>(sp));
      if (result.second) new_species++;
    }
  }
  return new_species;
}

// Accumulate abundances from a site into cumulative vector
inline void accumulate_abundances(
    const std::vector<int>& site_abundances,
    std::vector<int>& cumulative) {

  for (size_t sp = 0; sp < site_abundances.size(); sp++) {
    cumulative[sp] += site_abundances[sp];
  }
}

// Run kNN accumulation from a single seed
// Returns curve of cumulative species richness at each step
inline Curve knn_accumulate_single(
    const SiteSpeciesMatrix& species_pa,
    const DistanceMatrix& dist_mat,
    int seed) {

  size_t n_sites = species_pa.size();
  size_t n_species = species_pa[0].size();

  Curve curve(n_sites);
  VisitedFlags visited(n_sites, false);
  SpeciesSet species_seen;

  int current = seed;
  visited[current] = true;

  // Record species at starting site
  for (size_t sp = 0; sp < n_species; sp++) {
    if (species_pa[current][sp] > 0) {
      species_seen.insert(static_cast<int>(sp));
    }
  }
  curve[0] = static_cast<int>(species_seen.size());

  // Visit remaining sites
  for (size_t step = 1; step < n_sites; step++) {
    // Find nearest unvisited
    int nearest = find_nearest_unvisited(dist_mat, current, visited);

    // Move to nearest
    current = nearest;
    visited[current] = true;

    // Add new species
    for (size_t sp = 0; sp < n_species; sp++) {
      if (species_pa[current][sp] > 0) {
        species_seen.insert(static_cast<int>(sp));
      }
    }
    curve[step] = static_cast<int>(species_seen.size());
  }

  return curve;
}

// Run kNCN (k-Nearest Centroid Neighbor) accumulation from a single seed
// After each step, find nearest to centroid of visited sites
inline Curve kncn_accumulate_single(
    const SiteSpeciesMatrix& species_pa,
    const std::vector<double>& x,
    const std::vector<double>& y,
    int seed) {

  size_t n_sites = species_pa.size();
  size_t n_species = species_pa[0].size();

  Curve curve(n_sites);
  VisitedFlags visited(n_sites, false);
  SpeciesSet species_seen;

  // Centroid tracking
  double cx = x[seed];
  double cy = y[seed];
  int n_visited = 1;

  int current = seed;
  visited[current] = true;

  // Record species at starting site
  for (size_t sp = 0; sp < n_species; sp++) {
    if (species_pa[current][sp] > 0) {
      species_seen.insert(static_cast<int>(sp));
    }
  }
  curve[0] = static_cast<int>(species_seen.size());

  // Visit remaining sites
  for (size_t step = 1; step < n_sites; step++) {
    // Find nearest unvisited to centroid
    double min_dist = std::numeric_limits<double>::infinity();
    int nearest = -1;

    for (size_t j = 0; j < n_sites; j++) {
      if (!visited[j]) {
        double d = euclidean_distance(cx, cy, x[j], y[j]);
        if (d < min_dist) {
          min_dist = d;
          nearest = static_cast<int>(j);
        }
      }
    }

    // Move to nearest
    current = nearest;
    visited[current] = true;

    // Update centroid
    cx = (cx * n_visited + x[current]) / (n_visited + 1);
    cy = (cy * n_visited + y[current]) / (n_visited + 1);
    n_visited++;

    // Add new species
    for (size_t sp = 0; sp < n_species; sp++) {
      if (species_pa[current][sp] > 0) {
        species_seen.insert(static_cast<int>(sp));
      }
    }
    curve[step] = static_cast<int>(species_seen.size());
  }

  return curve;
}

// Run multiple kNN curves and return matrix
inline CurveMatrix knn_accumulate_multi(
    const SiteSpeciesMatrix& species_pa,
    const DistanceMatrix& dist_mat,
    const std::vector<int>& seeds) {

  size_t n_seeds = seeds.size();
  CurveMatrix curves(n_seeds);

  for (size_t s = 0; s < n_seeds; s++) {
    curves[s] = knn_accumulate_single(species_pa, dist_mat, seeds[s]);
  }

  return curves;
}

// Calculate mean curve from multiple curves
inline CurveDouble mean_curve(const CurveMatrix& curves) {
  if (curves.empty()) return CurveDouble();

  size_t n_steps = curves[0].size();
  size_t n_curves = curves.size();
  CurveDouble mean(n_steps, 0.0);

  for (size_t step = 0; step < n_steps; step++) {
    double sum = 0.0;
    for (size_t c = 0; c < n_curves; c++) {
      sum += curves[c][step];
    }
    mean[step] = sum / n_curves;
  }

  return mean;
}

// Calculate quantile of curves at each step
inline CurveDouble curve_quantile(const CurveMatrix& curves, double prob) {
  if (curves.empty()) return CurveDouble();

  size_t n_steps = curves[0].size();
  size_t n_curves = curves.size();
  CurveDouble result(n_steps, 0.0);

  for (size_t step = 0; step < n_steps; step++) {
    std::vector<double> values(n_curves);
    for (size_t c = 0; c < n_curves; c++) {
      values[c] = static_cast<double>(curves[c][step]);
    }
    std::sort(values.begin(), values.end());

    // Simple quantile calculation
    double idx = prob * (n_curves - 1);
    size_t lo = static_cast<size_t>(idx);
    size_t hi = std::min(lo + 1, n_curves - 1);
    double frac = idx - lo;

    result[step] = values[lo] * (1 - frac) + values[hi] * frac;
  }

  return result;
}

} // namespace spacc

#endif // SPACC_CORE_ACCUMULATION_H
