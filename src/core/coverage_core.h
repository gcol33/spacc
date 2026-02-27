// core/coverage_core.h
// Pure C++ coverage calculations for spacc
// Based on Chao & Jost (2012) sample coverage estimator
// No Rcpp dependencies - testable with Catch2

#ifndef SPACC_CORE_COVERAGE_H
#define SPACC_CORE_COVERAGE_H

#include "types.h"
#include <algorithm>

namespace spacc {

// Count singletons (species with exactly 1 individual)
template<typename T>
int count_singletons(const std::vector<T>& abundances) {
  int f1 = 0;
  for (size_t i = 0; i < abundances.size(); i++) {
    if (static_cast<int>(abundances[i]) == 1) f1++;
  }
  return f1;
}

// Count doubletons (species with exactly 2 individuals)
template<typename T>
int count_doubletons(const std::vector<T>& abundances) {
  int f2 = 0;
  for (size_t i = 0; i < abundances.size(); i++) {
    if (static_cast<int>(abundances[i]) == 2) f2++;
  }
  return f2;
}

// Count total individuals
template<typename T>
int count_total_individuals(const std::vector<T>& abundances) {
  int n = 0;
  for (size_t i = 0; i < abundances.size(); i++) {
    n += static_cast<int>(abundances[i]);
  }
  return n;
}

// Calculate Chao-Jost sample coverage estimator
// Returns value in [0, 1] where 1 = complete coverage
template<typename T>
double calc_chao_coverage(const std::vector<T>& abundances) {
  int n = 0;   // total individuals
  int f1 = 0;  // singletons
  int f2 = 0;  // doubletons

  for (size_t i = 0; i < abundances.size(); i++) {
    int a = static_cast<int>(abundances[i]);
    n += a;
    if (a == 1) f1++;
    if (a == 2) f2++;
  }

  if (n == 0) return 0.0;
  if (n == 1) return 0.0;

  double C;
  if (f2 > 0) {
    // Chao & Jost estimator with doubletons
    C = 1.0 - (static_cast<double>(f1) / n) *
              ((static_cast<double>(n - 1) * f1) /
               (static_cast<double>(n - 1) * f1 + 2.0 * f2));
  } else if (f1 > 0) {
    // Without doubletons
    C = 1.0 - (static_cast<double>(f1) / n) *
              (static_cast<double>(n - 1) / n);
  } else {
    C = 1.0;
  }

  // Clamp to [0, 1]
  return std::max(0.0, std::min(1.0, C));
}

// Calculate Chiu (2023) sample-based coverage estimator
// Uses incidence frequency counts (Q1, Q2) and abundance of uniques (G1)
// instead of individual-based singletons/doubletons.
// Near-unbiased for plot-based (spatially aggregated) data.
//
// Parameters:
//   abundances  - cumulative abundance per species (length = n_species)
//   incidences  - number of sampling units each species occurs in (length = n_species)
//   T           - number of accumulated sampling units (sites)
//
// Returns value in [0, 1] where 1 = complete coverage
template<typename T>
double calc_chiu_coverage(const std::vector<T>& abundances,
                          const std::vector<int>& incidences,
                          int T_sites) {
  if (T_sites <= 1) return 0.0;

  int n_plus = 0;  // total individuals across all accumulated sites
  int Q1 = 0;      // species detected in exactly 1 unit
  int Q2 = 0;      // species detected in exactly 2 units
  int G1 = 0;      // total abundance of species detected in exactly 1 unit

  for (size_t i = 0; i < abundances.size(); i++) {
    int a = static_cast<int>(abundances[i]);
    if (a == 0) continue;
    n_plus += a;
    int q = incidences[i];
    if (q == 1) {
      Q1++;
      G1 += a;
    }
    if (q == 2) Q2++;
  }

  if (n_plus == 0) return 0.0;

  // Cap Q1 and Q2 at minimum 1 to avoid division by zero
  int Q1_safe = std::max(1, Q1);
  int Q2_safe = std::max(1, Q2);

  // Chiu (2023) formula:
  // a0_hat = ((T-1)/T) * G1 * Q1 / (2 * Q2)
  // C = n+ / (n+ + a0_hat)
  // When Q1 == 0 (no uniques), coverage is 1.0
  if (Q1 == 0) return 1.0;

  double a0_hat = (static_cast<double>(T_sites - 1) / T_sites) *
                  static_cast<double>(G1) * Q1_safe /
                  (2.0 * Q2_safe);

  double C = static_cast<double>(n_plus) / (n_plus + a0_hat);

  return std::max(0.0, std::min(1.0, C));
}

// Calculate coverage deficit (1 - coverage)
// Represents the proportion of species not yet observed
template<typename T>
double calc_coverage_deficit(const std::vector<T>& abundances) {
  return 1.0 - calc_chao_coverage(abundances);
}

// Interpolate richness at target coverage level
// Given vectors of richness and coverage at each step,
// find the richness at a target coverage
inline double interpolate_richness_at_coverage(
    const std::vector<double>& richness,
    const std::vector<double>& coverage,
    double target) {

  size_t n = richness.size();
  if (n == 0) return 0.0;

  if (target > coverage[n-1]) {
    // Target above max coverage - cannot interpolate
    return -1.0;  // Indicates NA
  }

  if (target <= coverage[0]) {
    return richness[0];
  }

  // Find bracketing indices
  for (size_t i = 1; i < n; i++) {
    if (coverage[i] >= target) {
      // Linear interpolation
      double c0 = coverage[i-1];
      double c1 = coverage[i];
      double r0 = richness[i-1];
      double r1 = richness[i];

      if (c1 == c0) {
        return r1;
      } else {
        return r0 + (target - c0) * (r1 - r0) / (c1 - c0);
      }
    }
  }

  return richness[n-1];
}

} // namespace spacc

#endif // SPACC_CORE_COVERAGE_H
