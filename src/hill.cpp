// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace Rcpp;
using namespace RcppParallel;


// ============================================================================
// HILL NUMBER CALCULATIONS
// ============================================================================

//' Calculate Hill Number for a vector of abundances
//'
//' @param abundances Vector of species abundances
//' @param q Order of diversity (0 = richness, 1 = Shannon, 2 = Simpson)
//' @return Hill number (effective number of species)
// [[Rcpp::export]]
double calc_hill_number(NumericVector abundances, double q) {
  double N = 0.0;
  int S = 0;

  for (int i = 0; i < abundances.size(); i++) {
    if (abundances[i] > 0) {
      N += abundances[i];
      S++;
    }
  }

  if (N == 0 || S == 0) return 0.0;

  if (q == 0) {
    // Species richness
    return (double)S;
  } else if (std::abs(q - 1.0) < 1e-10) {
    // Shannon diversity (limit as q -> 1): exp(H)
    double H = 0.0;
    for (int i = 0; i < abundances.size(); i++) {
      if (abundances[i] > 0) {
        double p = abundances[i] / N;
        H -= p * std::log(p);
      }
    }
    return std::exp(H);
  } else {
    // General Hill number: (sum(p_i^q))^(1/(1-q))
    double sum_pq = 0.0;
    for (int i = 0; i < abundances.size(); i++) {
      if (abundances[i] > 0) {
        double p = abundances[i] / N;
        sum_pq += std::pow(p, q);
      }
    }
    return std::pow(sum_pq, 1.0 / (1.0 - q));
  }
}


// Internal version for std::vector
double calc_hill_internal(const std::vector<double>& abundances, double q) {
  double N = 0.0;
  int S = 0;

  for (size_t i = 0; i < abundances.size(); i++) {
    if (abundances[i] > 0) {
      N += abundances[i];
      S++;
    }
  }

  if (N == 0 || S == 0) return 0.0;

  if (q == 0) {
    return (double)S;
  } else if (std::abs(q - 1.0) < 1e-10) {
    double H = 0.0;
    for (size_t i = 0; i < abundances.size(); i++) {
      if (abundances[i] > 0) {
        double p = abundances[i] / N;
        H -= p * std::log(p);
      }
    }
    return std::exp(H);
  } else {
    double sum_pq = 0.0;
    for (size_t i = 0; i < abundances.size(); i++) {
      if (abundances[i] > 0) {
        double p = abundances[i] / N;
        sum_pq += std::pow(p, q);
      }
    }
    return std::pow(sum_pq, 1.0 / (1.0 - q));
  }
}


//' Single kNN Accumulation with Hill Numbers
//'
//' @param species_mat Integer matrix (sites x species) with abundances
//' @param dist_mat Numeric matrix of pairwise distances
//' @param seed Starting site index (0-based)
//' @param q_values Vector of Hill number orders to compute
//' @return Numeric matrix (length(q) x n_sites) of Hill numbers
// [[Rcpp::export]]
NumericMatrix cpp_knn_hill_single(IntegerMatrix species_mat,
                                   NumericMatrix dist_mat,
                                   int seed,
                                   NumericVector q_values) {
  int n_sites = species_mat.nrow();
  int n_species = species_mat.ncol();
  int n_q = q_values.size();

  NumericMatrix curves(n_q, n_sites);
  std::vector<bool> visited(n_sites, false);
  std::vector<double> cumulative(n_species, 0.0);

  int current = seed;
  visited[current] = true;

  // Add first site abundances
  for (int sp = 0; sp < n_species; sp++) {
    cumulative[sp] += species_mat(current, sp);
  }

  // Calculate Hill numbers for each q
  for (int qi = 0; qi < n_q; qi++) {
    curves(qi, 0) = calc_hill_internal(cumulative, q_values[qi]);
  }

  for (int step = 1; step < n_sites; step++) {
    // Find nearest unvisited
    double min_dist = R_PosInf;
    int next = -1;
    for (int j = 0; j < n_sites; j++) {
      if (!visited[j] && dist_mat(current, j) < min_dist) {
        min_dist = dist_mat(current, j);
        next = j;
      }
    }

    current = next;
    visited[current] = true;

    // Accumulate abundances
    for (int sp = 0; sp < n_species; sp++) {
      cumulative[sp] += species_mat(current, sp);
    }

    // Calculate Hill numbers
    for (int qi = 0; qi < n_q; qi++) {
      curves(qi, step) = calc_hill_internal(cumulative, q_values[qi]);
    }
  }

  return curves;
}


// Worker for parallel Hill number accumulation
struct HillKnnWorker : public Worker {
  const RMatrix<int> species_mat;
  const RMatrix<double> dist_mat;
  const RVector<int> seeds;
  const RVector<double> q_values;
  RMatrix<double> curves;  // 3D flattened: (n_seeds * n_q) x n_sites
  const int n_q;
  const int n_sites;

  HillKnnWorker(const IntegerMatrix& sp, const NumericMatrix& dm,
                const IntegerVector& s, const NumericVector& q,
                NumericMatrix& c, int nq, int ns)
    : species_mat(sp), dist_mat(dm), seeds(s), q_values(q),
      curves(c), n_q(nq), n_sites(ns) {}

  void operator()(std::size_t begin, std::size_t end) {
    int n_species = species_mat.ncol();

    for (std::size_t s = begin; s < end; s++) {
      std::vector<bool> visited(n_sites, false);
      std::vector<double> cumulative(n_species, 0.0);

      int current = seeds[s];
      visited[current] = true;

      for (int sp = 0; sp < n_species; sp++) {
        cumulative[sp] += species_mat(current, sp);
      }

      for (int qi = 0; qi < n_q; qi++) {
        curves(s * n_q + qi, 0) = calc_hill_internal(cumulative, q_values[qi]);
      }

      for (int step = 1; step < n_sites; step++) {
        double min_dist = R_PosInf;
        int next = -1;
        for (int j = 0; j < n_sites; j++) {
          if (!visited[j] && dist_mat(current, j) < min_dist) {
            min_dist = dist_mat(current, j);
            next = j;
          }
        }

        current = next;
        visited[current] = true;

        for (int sp = 0; sp < n_species; sp++) {
          cumulative[sp] += species_mat(current, sp);
        }

        for (int qi = 0; qi < n_q; qi++) {
          curves(s * n_q + qi, step) = calc_hill_internal(cumulative, q_values[qi]);
        }
      }
    }
  }
};


//' Parallel kNN Accumulation with Hill Numbers
//'
//' @param species_mat Integer matrix (sites x species) with abundances
//' @param dist_mat Numeric matrix of pairwise distances
//' @param n_seeds Number of random starting points
//' @param q_values Vector of Hill number orders
//' @param n_cores Number of cores
//' @param progress Show progress
//' @return List with curves for each q value
// [[Rcpp::export]]
List cpp_knn_hill_parallel(IntegerMatrix species_mat,
                           NumericMatrix dist_mat,
                           int n_seeds,
                           NumericVector q_values,
                           int n_cores = 1,
                           bool progress = false) {
  int n_sites = species_mat.nrow();
  int n_q = q_values.size();

  IntegerVector seeds = Rcpp::sample(n_sites, n_seeds, true) - 1;

  // Flattened output: (n_seeds * n_q) x n_sites
  NumericMatrix curves_flat(n_seeds * n_q, n_sites);

  if (n_cores > 1) {
    HillKnnWorker worker(species_mat, dist_mat, seeds, q_values,
                         curves_flat, n_q, n_sites);
    parallelFor(0, n_seeds, worker);
  } else {
    for (int s = 0; s < n_seeds; s++) {
      NumericMatrix single = cpp_knn_hill_single(species_mat, dist_mat,
                                                  seeds[s], q_values);
      for (int qi = 0; qi < n_q; qi++) {
        for (int st = 0; st < n_sites; st++) {
          curves_flat(s * n_q + qi, st) = single(qi, st);
        }
      }
    }
  }

  // Reshape to list of matrices, one per q
  List result(n_q);
  CharacterVector names(n_q);

  for (int qi = 0; qi < n_q; qi++) {
    NumericMatrix q_curves(n_seeds, n_sites);
    for (int s = 0; s < n_seeds; s++) {
      for (int st = 0; st < n_sites; st++) {
        q_curves(s, st) = curves_flat(s * n_q + qi, st);
      }
    }
    result[qi] = q_curves;

    // Name the list element
    char name[10];
    snprintf(name, sizeof(name), "q%.1f", q_values[qi]);
    names[qi] = name;
  }

  result.attr("names") = names;
  return result;
}
