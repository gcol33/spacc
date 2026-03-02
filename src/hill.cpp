// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include "core/hill_core.h"
using namespace Rcpp;
using namespace RcppParallel;


// ============================================================================
// HILL NUMBER CALCULATIONS
// ============================================================================

// [[Rcpp::export]]
double calc_hill_number(NumericVector abundances, double q) {
  std::vector<double> v(abundances.begin(), abundances.end());
  return spacc::calc_hill_number(v, q);
}


// Internal version for std::vector — delegates to core
inline double calc_hill_internal(const std::vector<double>& abundances, double q) {
  return spacc::calc_hill_number(abundances, q);
}


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


// ============================================================================
// HILL BETA DIVERSITY: gamma/alpha decomposition along accumulation
// ============================================================================

// Worker for parallel Hill beta accumulation
struct HillBetaKnnWorker : public Worker {
  const RMatrix<int> species_mat;
  const RMatrix<double> dist_mat;
  const RVector<int> seeds;
  const RVector<double> q_values;

  // Output: flattened (n_seeds * 3 * n_q) x n_sites
  // Layout: [seed][component=0:gamma,1:alpha,2:beta][q] x sites
  RMatrix<double> curves;
  const int n_q;
  const int n_sites;

  HillBetaKnnWorker(const IntegerMatrix& sp, const NumericMatrix& dm,
                    const IntegerVector& s, const NumericVector& q,
                    NumericMatrix& c, int nq, int ns)
    : species_mat(sp), dist_mat(dm), seeds(s), q_values(q),
      curves(c), n_q(nq), n_sites(ns) {}

  void operator()(std::size_t begin, std::size_t end) {
    int n_species = species_mat.ncol();

    for (std::size_t s = begin; s < end; s++) {
      std::vector<bool> visited(n_sites, false);
      std::vector<double> pooled(n_species, 0.0);

      // Track per-site Hill numbers for alpha computation
      // We store running power-mean accumulators per q
      // For q=1: sum of log(hill_i), for q!=1: sum of hill_i^(1-q)
      std::vector<double> alpha_accum(n_q, 0.0);

      int current = seeds[s];
      visited[current] = true;

      // Add first site
      std::vector<double> site_abund(n_species, 0.0);
      for (int sp = 0; sp < n_species; sp++) {
        double val = (double)species_mat(current, sp);
        pooled[sp] += val;
        site_abund[sp] = val;
      }

      // Gamma = Hill of pooled after 1 site
      // Alpha = Hill of first site (= gamma for k=1)
      // Beta = 1
      for (int qi = 0; qi < n_q; qi++) {
        double gamma_q = calc_hill_internal(pooled, q_values[qi]);
        double hill_site = calc_hill_internal(site_abund, q_values[qi]);

        // Initialize alpha accumulator
        if (std::abs(q_values[qi] - 1.0) < 1e-10) {
          alpha_accum[qi] = (hill_site > 0) ? std::log(hill_site) : 0.0;
        } else {
          alpha_accum[qi] = (hill_site > 0) ? std::pow(hill_site, 1.0 - q_values[qi]) : 0.0;
        }

        double alpha_q = hill_site;
        double beta_q = (alpha_q > 0) ? gamma_q / alpha_q : 1.0;

        int base = s * 3 * n_q;
        curves(base + 0 * n_q + qi, 0) = gamma_q;
        curves(base + 1 * n_q + qi, 0) = alpha_q;
        curves(base + 2 * n_q + qi, 0) = beta_q;
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

        // Compute per-site abundances and accumulate pooled
        std::fill(site_abund.begin(), site_abund.end(), 0.0);
        for (int sp = 0; sp < n_species; sp++) {
          double val = (double)species_mat(current, sp);
          pooled[sp] += val;
          site_abund[sp] = val;
        }

        int k = step + 1;  // number of sites so far

        for (int qi = 0; qi < n_q; qi++) {
          double gamma_q = calc_hill_internal(pooled, q_values[qi]);
          double hill_site = calc_hill_internal(site_abund, q_values[qi]);

          // Update alpha accumulator
          double alpha_q;
          if (std::abs(q_values[qi] - 1.0) < 1e-10) {
            // Geometric mean: exp(sum(log(hill_i)) / k)
            alpha_accum[qi] += (hill_site > 0) ? std::log(hill_site) : 0.0;
            alpha_q = std::exp(alpha_accum[qi] / k);
          } else {
            // Power mean: (sum(hill_i^(1-q)) / k)^(1/(1-q))
            alpha_accum[qi] += (hill_site > 0) ? std::pow(hill_site, 1.0 - q_values[qi]) : 0.0;
            alpha_q = std::pow(alpha_accum[qi] / k, 1.0 / (1.0 - q_values[qi]));
          }

          double beta_q = (alpha_q > 0) ? gamma_q / alpha_q : 1.0;

          int base = s * 3 * n_q;
          curves(base + 0 * n_q + qi, step) = gamma_q;
          curves(base + 1 * n_q + qi, step) = alpha_q;
          curves(base + 2 * n_q + qi, step) = beta_q;
        }
      }
    }
  }
};


// [[Rcpp::export]]
List cpp_knn_hill_beta_parallel(IntegerMatrix species_mat,
                                 NumericMatrix dist_mat,
                                 int n_seeds,
                                 NumericVector q_values,
                                 int n_cores = 1,
                                 bool progress = false) {
  int n_sites = species_mat.nrow();
  int n_q = q_values.size();

  IntegerVector seeds = Rcpp::sample(n_sites, n_seeds, true) - 1;

  // Flattened: (n_seeds * 3 * n_q) x n_sites
  NumericMatrix curves_flat(n_seeds * 3 * n_q, n_sites);

  if (n_cores > 1) {
    HillBetaKnnWorker worker(species_mat, dist_mat, seeds, q_values,
                              curves_flat, n_q, n_sites);
    parallelFor(0, n_seeds, worker);
  } else {
    HillBetaKnnWorker worker(species_mat, dist_mat, seeds, q_values,
                              curves_flat, n_q, n_sites);
    worker(0, n_seeds);
  }

  // Reshape to named lists: $gamma, $alpha, $beta
  // Each is a list of matrices (one per q)
  List gamma_list(n_q), alpha_list(n_q), beta_list(n_q);
  CharacterVector q_names(n_q);

  for (int qi = 0; qi < n_q; qi++) {
    NumericMatrix gamma_mat(n_seeds, n_sites);
    NumericMatrix alpha_mat(n_seeds, n_sites);
    NumericMatrix beta_mat(n_seeds, n_sites);

    for (int s = 0; s < n_seeds; s++) {
      int base = s * 3 * n_q;
      for (int st = 0; st < n_sites; st++) {
        gamma_mat(s, st) = curves_flat(base + 0 * n_q + qi, st);
        alpha_mat(s, st) = curves_flat(base + 1 * n_q + qi, st);
        beta_mat(s, st)  = curves_flat(base + 2 * n_q + qi, st);
      }
    }

    gamma_list[qi] = gamma_mat;
    alpha_list[qi] = alpha_mat;
    beta_list[qi]  = beta_mat;

    char name[10];
    snprintf(name, sizeof(name), "q%.1f", q_values[qi]);
    q_names[qi] = name;
  }

  gamma_list.attr("names") = q_names;
  alpha_list.attr("names") = q_names;
  beta_list.attr("names")  = q_names;

  return List::create(
    Named("gamma") = gamma_list,
    Named("alpha") = alpha_list,
    Named("beta")  = beta_list
  );
}
