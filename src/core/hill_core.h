// core/hill_core.h
// Pure C++ Hill number calculations for spacc
// No Rcpp dependencies - testable with Catch2

#ifndef SPACC_CORE_HILL_H
#define SPACC_CORE_HILL_H

#include "types.h"
#include <cmath>
#include <numeric>

namespace spacc {

// Calculate Hill number of order q
// q = 0: Species richness
// q = 1: Exponential of Shannon entropy (limit as q -> 1)
// q = 2: Inverse Simpson concentration
template<typename T>
double calc_hill_number(const std::vector<T>& abundances, double q) {
  double N = 0.0;
  int S = 0;

  for (size_t i = 0; i < abundances.size(); i++) {
    if (abundances[i] > 0) {
      N += static_cast<double>(abundances[i]);
      S++;
    }
  }

  if (N == 0 || S == 0) return 0.0;

  if (q == 0) {
    // Species richness
    return static_cast<double>(S);
  } else if (std::abs(q - 1.0) < 1e-10) {
    // Shannon diversity (limit as q -> 1): exp(H)
    double H = 0.0;
    for (size_t i = 0; i < abundances.size(); i++) {
      if (abundances[i] > 0) {
        double p = static_cast<double>(abundances[i]) / N;
        H -= p * std::log(p);
      }
    }
    return std::exp(H);
  } else {
    // General Hill number: (sum(p_i^q))^(1/(1-q))
    double sum_pq = 0.0;
    for (size_t i = 0; i < abundances.size(); i++) {
      if (abundances[i] > 0) {
        double p = static_cast<double>(abundances[i]) / N;
        sum_pq += std::pow(p, q);
      }
    }
    return std::pow(sum_pq, 1.0 / (1.0 - q));
  }
}

// Calculate Shannon entropy (H)
template<typename T>
double calc_shannon_entropy(const std::vector<T>& abundances) {
  double N = 0.0;
  for (size_t i = 0; i < abundances.size(); i++) {
    if (abundances[i] > 0) {
      N += static_cast<double>(abundances[i]);
    }
  }

  if (N == 0) return 0.0;

  double H = 0.0;
  for (size_t i = 0; i < abundances.size(); i++) {
    if (abundances[i] > 0) {
      double p = static_cast<double>(abundances[i]) / N;
      H -= p * std::log(p);
    }
  }

  return H;
}

// Calculate Simpson index (D = sum(p_i^2))
template<typename T>
double calc_simpson_index(const std::vector<T>& abundances) {
  double N = 0.0;
  for (size_t i = 0; i < abundances.size(); i++) {
    if (abundances[i] > 0) {
      N += static_cast<double>(abundances[i]);
    }
  }

  if (N == 0) return 0.0;

  double D = 0.0;
  for (size_t i = 0; i < abundances.size(); i++) {
    if (abundances[i] > 0) {
      double p = static_cast<double>(abundances[i]) / N;
      D += p * p;
    }
  }

  return D;
}

// Calculate inverse Simpson (1/D) = Hill number q=2
template<typename T>
double calc_inverse_simpson(const std::vector<T>& abundances) {
  double D = calc_simpson_index(abundances);
  return (D > 0) ? 1.0 / D : 0.0;
}

// Count species richness (number of species with abundance > 0)
template<typename T>
int count_species_richness(const std::vector<T>& abundances) {
  int S = 0;
  for (size_t i = 0; i < abundances.size(); i++) {
    if (abundances[i] > 0) S++;
  }
  return S;
}

} // namespace spacc

#endif // SPACC_CORE_HILL_H
