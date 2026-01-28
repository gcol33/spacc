// core/beta_core.h
// Pure C++ beta diversity calculations for spacc
// Based on Baselga (2010) partitioning
// No Rcpp dependencies - testable with Catch2

#ifndef SPACC_CORE_BETA_H
#define SPACC_CORE_BETA_H

#include "types.h"
#include <algorithm>
#include <tuple>

namespace spacc {

// Structure to hold beta diversity components
struct BetaComponents {
  double total;
  double turnover;
  double nestedness;

  BetaComponents() : total(0.0), turnover(0.0), nestedness(0.0) {}
  BetaComponents(double t, double tu, double n) : total(t), turnover(tu), nestedness(n) {}
};

// Calculate Sorensen-based beta diversity and its components
// a = shared species, b = unique to set1, c = unique to set2
inline BetaComponents calc_beta_sorensen(int a, int b, int c) {
  BetaComponents result;

  if (a + b + c == 0) return result;

  // Total beta (Sorensen dissimilarity)
  result.total = static_cast<double>(b + c) / (2.0 * a + b + c);

  // Turnover component (Simpson)
  int min_bc = std::min(b, c);
  if (a + min_bc > 0) {
    result.turnover = static_cast<double>(min_bc) / (a + min_bc);
  }

  // Nestedness component (total - turnover)
  result.nestedness = result.total - result.turnover;

  return result;
}

// Calculate Jaccard-based beta diversity and its components
inline BetaComponents calc_beta_jaccard(int a, int b, int c) {
  BetaComponents result;

  if (a + b + c == 0) return result;

  // Total beta (Jaccard dissimilarity)
  result.total = static_cast<double>(b + c) / (a + b + c);

  // Turnover component
  int min_bc = std::min(b, c);
  if (a + min_bc > 0) {
    result.turnover = 2.0 * static_cast<double>(min_bc) / (a + 2.0 * min_bc);
  }

  // Nestedness component
  result.nestedness = result.total - result.turnover;

  return result;
}

// Count shared and unique species between two sets
// Returns tuple<a, b, c> where:
//   a = shared species
//   b = unique to set1
//   c = unique to set2
inline std::tuple<int, int, int> count_abc(
    const SpeciesSet& set1,
    const SpeciesSet& set2) {

  int a = 0, b = 0, c = 0;

  for (int sp : set1) {
    if (set2.count(sp) > 0) {
      a++;
    } else {
      b++;
    }
  }

  for (int sp : set2) {
    if (set1.count(sp) == 0) {
      c++;
    }
  }

  return std::make_tuple(a, b, c);
}

// Convenience function that takes vectors instead of sets
inline std::tuple<int, int, int> count_abc_vectors(
    const std::vector<int>& vec1,
    const std::vector<int>& vec2) {

  SpeciesSet set1(vec1.begin(), vec1.end());
  SpeciesSet set2(vec2.begin(), vec2.end());
  return count_abc(set1, set2);
}

} // namespace spacc

#endif // SPACC_CORE_BETA_H
