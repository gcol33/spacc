// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include "core/hill_core.h"
#include "core/coverage_core.h"
using namespace Rcpp;
using namespace RcppParallel;


// Thin wrappers around core implementations (single source of truth)
static inline double hc_calc_hill(const std::vector<double>& abundances, double q) {
  return spacc::calc_hill_number(abundances, q);
}

static inline double hc_calc_coverage(const std::vector<int>& abundances) {
  return spacc::calc_chao_coverage(abundances);
}


// ============================================================================
// COMBINED HILL + COVERAGE WORKER
// ============================================================================

struct HillCoverageOrderWorker : public Worker {
  const RMatrix<int> species_mat;
  const RMatrix<int> orders;
  const RVector<double> q_values;

  // Output matrices
  RMatrix<double> hills_flat;  // (n_seeds * n_q) x n_sites
  RMatrix<double> coverage;    // n_seeds x n_sites

  const int n_q;
  const int n_sites;

  HillCoverageOrderWorker(const IntegerMatrix& sp, const IntegerMatrix& ord,
                        const NumericVector& q,
                        NumericMatrix& h, NumericMatrix& c,
                        int nq, int ns)
    : species_mat(sp), orders(ord), q_values(q),
      hills_flat(h), coverage(c), n_q(nq), n_sites(ns) {}

  void operator()(std::size_t begin, std::size_t end) {
    int n_species = species_mat.ncol();

    for (std::size_t s = begin; s < end; s++) {
      std::vector<int> cumul_int(n_species, 0);
      std::vector<double> cumul_dbl(n_species, 0.0);

      int current = orders(s, 0);

      // Add first site
      for (int sp = 0; sp < n_species; sp++) {
        cumul_int[sp] += species_mat(current, sp);
        cumul_dbl[sp] += (double)species_mat(current, sp);
      }

      // Compute coverage and Hill for first step
      coverage(s, 0) = hc_calc_coverage(cumul_int);
      for (int qi = 0; qi < n_q; qi++) {
        hills_flat(s * n_q + qi, 0) = hc_calc_hill(cumul_dbl, q_values[qi]);
      }

      for (int step = 1; step < n_sites; step++) {
        current = orders(s, step);

        // Accumulate
        for (int sp = 0; sp < n_species; sp++) {
          cumul_int[sp] += species_mat(current, sp);
          cumul_dbl[sp] += (double)species_mat(current, sp);
        }

        coverage(s, step) = hc_calc_coverage(cumul_int);
        for (int qi = 0; qi < n_q; qi++) {
          hills_flat(s * n_q + qi, step) = hc_calc_hill(cumul_dbl, q_values[qi]);
        }
      }
    }
  }
};


// [[Rcpp::export]]
List cpp_order_hill_coverage_parallel(IntegerMatrix species_mat,
                                      IntegerMatrix orders,
                                      NumericVector q_values,
                                      int n_cores = 1,
                                      bool progress = false) {
  int n_sites = species_mat.nrow();
  int n_q = q_values.size();
  int n_seeds = orders.nrow();

  NumericMatrix hills_flat(n_seeds * n_q, n_sites);
  NumericMatrix coverage(n_seeds, n_sites);

  if (n_cores > 1) {
    HillCoverageOrderWorker worker(species_mat, orders, q_values,
                                  hills_flat, coverage, n_q, n_sites);
    parallelFor(0, n_seeds, worker);
  } else {
    // Single-threaded fallback using the worker directly
    HillCoverageOrderWorker worker(species_mat, orders, q_values,
                                  hills_flat, coverage, n_q, n_sites);
    worker(0, n_seeds);
  }

  // Reshape hills to list of matrices, one per q
  List hills(n_q);
  CharacterVector names(n_q);

  for (int qi = 0; qi < n_q; qi++) {
    NumericMatrix q_curves(n_seeds, n_sites);
    for (int s = 0; s < n_seeds; s++) {
      for (int st = 0; st < n_sites; st++) {
        q_curves(s, st) = hills_flat(s * n_q + qi, st);
      }
    }
    hills[qi] = q_curves;

    char name[10];
    snprintf(name, sizeof(name), "q%.1f", q_values[qi]);
    names[qi] = name;
  }

  hills.attr("names") = names;

  return List::create(
    Named("hills") = hills,
    Named("coverage") = coverage
  );
}
