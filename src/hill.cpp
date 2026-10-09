// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <cmath>
#include <vector>
#include "core/hill_core.h"
#include "core/coverage_core.h"
#include "parallel_utils.h"
using namespace Rcpp;
using namespace RcppParallel;


// [[Rcpp::export]]
double calc_hill_number(NumericVector abundances, double q) {
  std::vector<double> v(abundances.begin(), abundances.end());
  return spacc::calc_hill_number(v, q);
}


// ============================================================================
// Hill numbers of the pooled community along each order
// ============================================================================

struct HillMetric {
  RMatrix<double> out;       // (n_seeds * n_q) x n_sites
  RVector<double> q;

  HillMetric(NumericMatrix& out_, const NumericVector& q_) : out(out_), q(q_) {}

  struct Scratch {};
  Scratch scratch(int) const { return Scratch(); }

  void at_step(std::size_t s, int step, int, const std::vector<int>& cumulative,
               Scratch&) {
    const std::size_t n_q = q.size();
    for (std::size_t qi = 0; qi < n_q; qi++) {
      out(s * n_q + qi, step) = spacc::calc_hill_number(cumulative, q[qi]);
    }
  }
};

// [[Rcpp::export]]
List cpp_order_hill_parallel(IntegerMatrix species_mat, IntegerMatrix orders,
                             NumericVector q_values, int n_cores = 1) {
  const int n_seeds = orders.nrow(), n_q = q_values.size();
  NumericMatrix flat(n_seeds * n_q, species_mat.nrow());
  HillMetric metric(flat, q_values);
  run_order_metric<int>(species_mat, orders, metric, n_cores);
  return split_components(flat, n_seeds, n_q, q_names(q_values));
}


// ============================================================================
// Hill gamma / alpha / beta decomposition (Jost 2007) along each order
// ============================================================================

struct HillBetaMetric {
  RMatrix<int> species;
  RMatrix<double> out;       // (n_seeds * 3 * n_q) x n_sites: [gamma | alpha | beta] per seed
  RVector<double> q;

  HillBetaMetric(const IntegerMatrix& species_, NumericMatrix& out_,
                 const NumericVector& q_)
    : species(species_), out(out_), q(q_) {}

  struct Scratch {
    std::vector<double> alpha_accum;   // per q: sum of log(hill_i) (q = 1) or hill_i^(1-q)
    std::vector<int> site_abund;
  };
  Scratch scratch(int n_species) const {
    Scratch sc;
    sc.alpha_accum.assign(q.size(), 0.0);
    sc.site_abund.assign(n_species, 0);
    return sc;
  }

  void at_step(std::size_t s, int step, int site, const std::vector<int>& pooled,
               Scratch& sc) {
    const int n_species = species.ncol();
    for (int sp = 0; sp < n_species; sp++) sc.site_abund[sp] = species(site, sp);

    const std::size_t n_q = q.size();
    const int k = step + 1;
    const std::size_t base = s * 3 * n_q;
    for (std::size_t qi = 0; qi < n_q; qi++) {
      const double gamma_q = spacc::calc_hill_number(pooled, q[qi]);
      const double hill_site = spacc::calc_hill_number(sc.site_abund, q[qi]);
      const bool shannon = std::abs(q[qi] - 1.0) < 1e-10;

      if (hill_site > 0) {
        sc.alpha_accum[qi] += shannon ? std::log(hill_site)
                                      : std::pow(hill_site, 1.0 - q[qi]);
      }
      double alpha_q;
      if (step == 0) {
        alpha_q = hill_site;
      } else if (shannon) {
        alpha_q = std::exp(sc.alpha_accum[qi] / k);
      } else {
        alpha_q = std::pow(sc.alpha_accum[qi] / k, 1.0 / (1.0 - q[qi]));
      }
      const double beta_q = (alpha_q > 0) ? gamma_q / alpha_q : 1.0;

      out(base + 0 * n_q + qi, step) = gamma_q;
      out(base + 1 * n_q + qi, step) = alpha_q;
      out(base + 2 * n_q + qi, step) = beta_q;
    }
  }
};

// [[Rcpp::export]]
List cpp_order_hill_beta_parallel(IntegerMatrix species_mat, IntegerMatrix orders,
                                  NumericVector q_values, int n_cores = 1) {
  const int n_seeds = orders.nrow(), n_q = q_values.size(), n_sites = species_mat.nrow();
  NumericMatrix flat(n_seeds * 3 * n_q, n_sites);
  HillBetaMetric metric(species_mat, flat, q_values);
  run_order_metric<int>(species_mat, orders, metric, n_cores);

  // Regroup the per-seed [gamma | alpha | beta] blocks into one flat matrix per
  // component, then split each by q.
  CharacterVector names = q_names(q_values);
  List out(3);
  const char* comp_names[3] = {"gamma", "alpha", "beta"};
  for (int comp = 0; comp < 3; comp++) {
    NumericMatrix comp_flat(n_seeds * n_q, n_sites);
    for (int s = 0; s < n_seeds; s++) {
      for (int qi = 0; qi < n_q; qi++) {
        for (int st = 0; st < n_sites; st++) {
          comp_flat(s * n_q + qi, st) = flat(s * 3 * n_q + comp * n_q + qi, st);
        }
      }
    }
    out[comp] = split_components(comp_flat, n_seeds, n_q, names);
  }
  out.attr("names") = CharacterVector::create(comp_names[0], comp_names[1], comp_names[2]);
  return out;
}


// ============================================================================
// Hill numbers with Chao-Jost coverage along each order
// ============================================================================

struct HillCoverageMetric {
  RMatrix<double> hills;     // (n_seeds * n_q) x n_sites
  RMatrix<double> coverage;  // n_seeds x n_sites
  RVector<double> q;

  HillCoverageMetric(NumericMatrix& hills_, NumericMatrix& coverage_,
                     const NumericVector& q_)
    : hills(hills_), coverage(coverage_), q(q_) {}

  struct Scratch {};
  Scratch scratch(int) const { return Scratch(); }

  void at_step(std::size_t s, int step, int, const std::vector<int>& cumulative,
               Scratch&) {
    coverage(s, step) = spacc::calc_chao_coverage(cumulative);
    const std::size_t n_q = q.size();
    for (std::size_t qi = 0; qi < n_q; qi++) {
      hills(s * n_q + qi, step) = spacc::calc_hill_number(cumulative, q[qi]);
    }
  }
};

// [[Rcpp::export]]
List cpp_order_hill_coverage_parallel(IntegerMatrix species_mat, IntegerMatrix orders,
                                      NumericVector q_values, int n_cores = 1) {
  const int n_seeds = orders.nrow(), n_q = q_values.size();
  NumericMatrix hills_flat(n_seeds * n_q, species_mat.nrow());
  NumericMatrix coverage(n_seeds, species_mat.nrow());
  HillCoverageMetric metric(hills_flat, coverage, q_values);
  run_order_metric<int>(species_mat, orders, metric, n_cores);
  return List::create(
    Named("hills") = split_components(hills_flat, n_seeds, n_q, q_names(q_values)),
    Named("coverage") = coverage
  );
}
