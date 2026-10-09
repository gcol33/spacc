// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <cmath>
#include <vector>
#include "core/hill_core.h"
#include "core/coverage_core.h"
#include "core/inext_core.h"
#include "parallel_utils.h"
using namespace Rcpp;
using namespace RcppParallel;


// [[Rcpp::export]]
double calc_hill_number(NumericVector abundances, double q) {
  std::vector<double> v(abundances.begin(), abundances.end());
  return spacc::calc_hill_number(v, q);
}


// Rarefied / extrapolated Hill numbers of one abundance sample at sample
// sizes m (one row per q) and their expected coverage.
// [[Rcpp::export]]
List calc_hill_at_size(NumericVector abundances, NumericVector m, NumericVector q) {
  spacc::AbundanceSample s(std::vector<double>(abundances.begin(), abundances.end()));
  NumericMatrix hill(q.size(), m.size());
  NumericVector coverage(m.size());
  for (int j = 0; j < m.size(); j++) {
    coverage[j] = spacc::coverage_at_size(s, m[j]);
    for (int i = 0; i < q.size(); i++) hill(i, j) = spacc::hill_at_size(s, m[j], q[i]);
  }
  return List::create(Named("hill") = hill, Named("coverage") = coverage);
}


// Hill numbers of one abundance sample at target coverages C (one row per q)
// and the sample sizes those coverages correspond to.
// [[Rcpp::export]]
List calc_hill_at_coverage(NumericVector abundances, NumericVector C, NumericVector q) {
  spacc::AbundanceSample s(std::vector<double>(abundances.begin(), abundances.end()));
  NumericMatrix hill(q.size(), C.size());
  NumericVector size(C.size());
  for (int j = 0; j < C.size(); j++) {
    size[j] = spacc::size_at_coverage(s, C[j]);
    for (int i = 0; i < q.size(); i++) hill(i, j) = spacc::hill_at_size(s, size[j], q[i]);
  }
  return List::create(Named("hill") = hill, Named("size") = size);
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
// Hill gamma / alpha / beta decomposition along each order: Jost (2007) on the
// observed samples, or the iNEXT.beta3D coverage standardisation (Chao et al.
// 2023)
// ============================================================================

// Dissimilarities of a beta value over N assemblages (Chao et al. 2023):
// 1 - C_qN (Sorensen-type), 1 - U_qN (Jaccard-type), 1 - V_qN, 1 - S_qN.
static inline void beta_dissimilarities(double beta, double N, double q, double* d) {
  if (N <= 1) {
    for (int i = 0; i < 4; i++) d[i] = NA_REAL;
    return;
  }
  if (q == 1.0) {
    d[0] = d[1] = std::log(beta) / std::log(N);
  } else {
    d[0] = (std::pow(beta, 1.0 - q) - 1.0) / (std::pow(N, 1.0 - q) - 1.0);
    d[1] = (std::pow(beta, q - 1.0) - 1.0) / (std::pow(N, q - 1.0) - 1.0);
  }
  d[2] = (beta - 1.0) / (N - 1.0);
  d[3] = (1.0 / beta - 1.0) / (1.0 / N - 1.0);
}

struct HillBetaMetric {
  static constexpr int kObserved = 3;       // gamma, alpha, beta
  static constexpr int kStandardized = 7;   // + four dissimilarities

  RMatrix<int> species;
  RMatrix<double> out;       // (n_seeds * n_comp * n_q) x n_sites, per seed [comp][q]
  RMatrix<double> coverage;  // n_seeds x n_sites, standardisation coverage per step
  RVector<double> q;
  bool standardize;
  double target;             // fixed coverage, or NaN for the per-step reference
  int n_comp;

  HillBetaMetric(const IntegerMatrix& species_, NumericMatrix& out_,
                 NumericMatrix& coverage_, const NumericVector& q_,
                 bool standardize_, double target_)
    : species(species_), out(out_), coverage(coverage_), q(q_),
      standardize(standardize_), target(target_),
      n_comp(standardize_ ? kStandardized : kObserved) {}

  struct Scratch {
    std::vector<double> alpha_accum;   // per q: sum of log(hill_i) (q = 1) or hill_i^(1-q)
    std::vector<int> site_abund;
    spacc::AbundanceSample joint;      // species-by-site cells of the accumulated sites
  };
  Scratch scratch(int n_species) const {
    Scratch sc;
    sc.alpha_accum.assign(q.size(), 0.0);
    sc.site_abund.assign(n_species, 0);
    return sc;
  }

  void at_step(std::size_t s, int step, int site, const std::vector<int>& pooled,
               Scratch& sc) {
    if (standardize) standardized_step(s, step, site, pooled, sc);
    else observed_step(s, step, site, pooled, sc);
  }

  // Gamma of the pooled sample and alpha of the joint sample (Chiu et al.
  // 2014: joint Hill number / N), both at the same coverage C. Without a
  // fixed target, C is the lower of the two observed coverages, so both are
  // read off the observed sample or its rarefaction.
  void standardized_step(std::size_t s, int step, int site, const std::vector<int>& pooled,
                         Scratch& sc) {
    const int n_species = species.ncol();
    for (int sp = 0; sp < n_species; sp++) sc.joint.add(species(site, sp));
    spacc::AbundanceSample gamma_sample(pooled);

    const double N = step + 1;
    const double C = std::isnan(target)
      ? std::min(spacc::coverage_at_size(gamma_sample, gamma_sample.n),
                 spacc::coverage_at_size(sc.joint, sc.joint.n))
      : target;
    coverage(s, step) = C;
    const double m_gamma = spacc::size_at_coverage(gamma_sample, C);
    const double m_alpha = spacc::size_at_coverage(sc.joint, C);

    const std::size_t n_q = q.size();
    const std::size_t base = s * n_comp * n_q;
    for (std::size_t qi = 0; qi < n_q; qi++) {
      const double gamma_q = spacc::hill_at_size(gamma_sample, m_gamma, q[qi]);
      const double alpha_q = spacc::hill_at_size(sc.joint, m_alpha, q[qi]) / N;
      const double beta_q = gamma_q / alpha_q;
      double d[4];
      beta_dissimilarities(beta_q, N, q[qi], d);
      out(base + 0 * n_q + qi, step) = gamma_q;
      out(base + 1 * n_q + qi, step) = alpha_q;
      out(base + 2 * n_q + qi, step) = beta_q;
      for (int c = 0; c < 4; c++) out(base + (3 + c) * n_q + qi, step) = d[c];
    }
  }

  void observed_step(std::size_t s, int step, int site, const std::vector<int>& pooled,
                     Scratch& sc) {
    const int n_species = species.ncol();
    for (int sp = 0; sp < n_species; sp++) sc.site_abund[sp] = species(site, sp);

    const std::size_t n_q = q.size();
    const int k = step + 1;
    const std::size_t base = s * kObserved * n_q;
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

// coverage_target: NA for the per-step reference coverage, used only when
// standardize is true.
// [[Rcpp::export]]
List cpp_order_hill_beta_parallel(IntegerMatrix species_mat, IntegerMatrix orders,
                                  NumericVector q_values, int n_cores = 1,
                                  bool standardize = false,
                                  double coverage_target = NA_REAL) {
  const int n_seeds = orders.nrow(), n_q = q_values.size(), n_sites = species_mat.nrow();
  const int n_comp = standardize ? HillBetaMetric::kStandardized : HillBetaMetric::kObserved;
  NumericMatrix flat(n_seeds * n_comp * n_q, n_sites);
  NumericMatrix coverage(standardize ? n_seeds : 0, standardize ? n_sites : 0);
  HillBetaMetric metric(species_mat, flat, coverage, q_values, standardize,
                        ISNAN(coverage_target) ? std::nan("") : coverage_target);
  run_order_metric<int>(species_mat, orders, metric, n_cores);

  const char* all_names[HillBetaMetric::kStandardized] =
    {"gamma", "alpha", "beta", "1-CqN", "1-UqN", "1-VqN", "1-SqN"};
  CharacterVector comp_names(n_comp);
  for (int c = 0; c < n_comp; c++) comp_names[c] = all_names[c];
  List out = split_component_blocks(flat, n_seeds, comp_names, q_names(q_values));
  if (standardize) {
    out.push_back(coverage, "coverage");
  }
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
