// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include "core/coverage_core.h"
#include "parallel_utils.h"
using namespace Rcpp;
using namespace RcppParallel;


// Chao & Jost (2012) individual-based sample coverage.
// [[Rcpp::export]]
double calc_coverage(NumericVector abundances) {
  std::vector<int> v(abundances.size());
  for (int i = 0; i < abundances.size(); i++) v[i] = static_cast<int>(abundances[i]);
  return spacc::calc_chao_coverage(v);
}


// Chiu (2023) sample-based coverage from incidence frequencies.
// [[Rcpp::export]]
double calc_chiu_coverage(NumericVector abundances, IntegerVector incidences, int T_sites) {
  std::vector<int> abund(abundances.size());
  for (int i = 0; i < abundances.size(); i++) abund[i] = static_cast<int>(abundances[i]);
  std::vector<int> incid(incidences.begin(), incidences.end());
  return spacc::calc_chiu_coverage(abund, incid, T_sites);
}


// Richness, individuals and coverage of the pooled community along each order.
// coverage_type: 0 = Chao-Jost (individual-based), 1 = Chiu (sample-based)
struct CoverageMetric {
  RMatrix<int> species;
  RMatrix<int> richness;
  RMatrix<int> individuals;
  RMatrix<double> coverage;
  int coverage_type;

  CoverageMetric(const IntegerMatrix& species_, IntegerMatrix& richness_,
                 IntegerMatrix& individuals_, NumericMatrix& coverage_, int type)
    : species(species_), richness(richness_), individuals(individuals_),
      coverage(coverage_), coverage_type(type) {}

  struct Scratch {
    std::vector<int> incidence;   // sites per species so far
    int total_individuals;
  };
  Scratch scratch(int n_species) const {
    Scratch sc;
    sc.incidence.assign(n_species, 0);
    sc.total_individuals = 0;
    return sc;
  }

  void at_step(std::size_t s, int step, int site, const std::vector<int>& cumulative,
               Scratch& sc) {
    const int n_species = species.ncol();
    int n_present = 0;
    for (int sp = 0; sp < n_species; sp++) {
      const int a = species(site, sp);
      sc.total_individuals += a;
      if (a > 0) sc.incidence[sp]++;
      if (cumulative[sp] > 0) n_present++;
    }
    richness(s, step) = n_present;
    individuals(s, step) = sc.total_individuals;
    coverage(s, step) = coverage_type == 1
      ? spacc::calc_chiu_coverage(cumulative, sc.incidence, step + 1)
      : spacc::calc_chao_coverage(cumulative);
  }
};

// [[Rcpp::export]]
List cpp_order_coverage_parallel(IntegerMatrix species_mat, IntegerMatrix orders,
                                 int coverage_type = 0, int n_cores = 1) {
  const int n_seeds = orders.nrow(), n_sites = species_mat.nrow();
  IntegerMatrix richness(n_seeds, n_sites);
  IntegerMatrix individuals(n_seeds, n_sites);
  NumericMatrix coverage(n_seeds, n_sites);
  CoverageMetric metric(species_mat, richness, individuals, coverage, coverage_type);
  run_order_metric<int>(species_mat, orders, metric, n_cores);
  return List::create(
    Named("richness") = richness,
    Named("individuals") = individuals,
    Named("coverage") = coverage
  );
}


// Linear interpolation of a richness curve at target coverage values; NA
// where the target exceeds the final coverage.
// [[Rcpp::export]]
NumericVector interpolate_at_coverage(NumericVector richness, NumericVector coverage,
                                      NumericVector targets) {
  std::vector<double> r(richness.begin(), richness.end());
  std::vector<double> c(coverage.begin(), coverage.end());
  NumericVector result(targets.size());
  for (int t = 0; t < targets.size(); t++) {
    double value = spacc::interpolate_richness_at_coverage(r, c, targets[t]);
    result[t] = value < 0 ? NA_REAL : value;
  }
  return result;
}
