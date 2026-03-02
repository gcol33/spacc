// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <algorithm>
#include "core/coverage_core.h"
using namespace Rcpp;
using namespace RcppParallel;


// ============================================================================
// COVERAGE CALCULATIONS (Chao & Jost 2012)
// ============================================================================

// [[Rcpp::export]]
double calc_coverage(NumericVector abundances) {
  std::vector<int> v(abundances.size());
  for (int i = 0; i < abundances.size(); i++) {
    v[i] = (int)abundances[i];
  }
  return spacc::calc_chao_coverage(v);
}


// Internal version for std::vector — delegates to core
inline double calc_coverage_internal(const std::vector<int>& abundances) {
  return spacc::calc_chao_coverage(abundances);
}


// Internal Chiu (2023) sample-based coverage — delegates to core
inline double calc_chiu_coverage_internal(const std::vector<int>& abundances,
                                          const std::vector<int>& incidences,
                                          int T_sites) {
  return spacc::calc_chiu_coverage(abundances, incidences, T_sites);
}


// [[Rcpp::export]]
double calc_chiu_coverage(NumericVector abundances,
                          IntegerVector incidences,
                          int T_sites) {
  std::vector<int> abund(abundances.size());
  std::vector<int> incid(incidences.size());
  for (int i = 0; i < abundances.size(); i++) {
    abund[i] = (int)abundances[i];
  }
  for (int i = 0; i < incidences.size(); i++) {
    incid[i] = incidences[i];
  }
  return calc_chiu_coverage_internal(abund, incid, T_sites);
}


// coverage_type: 0 = Chao-Jost (individual-based), 1 = Chiu (sample-based)
// [[Rcpp::export]]
List cpp_knn_coverage_single(IntegerMatrix species_mat,
                              NumericMatrix dist_mat,
                              int seed,
                              int coverage_type = 0) {
  int n_sites = species_mat.nrow();
  int n_species = species_mat.ncol();

  IntegerVector richness(n_sites);
  IntegerVector individuals(n_sites);
  NumericVector coverage(n_sites);

  std::vector<bool> visited(n_sites, false);
  std::vector<int> cumulative(n_species, 0);
  std::vector<int> incidence(n_species, 0);  // per-species site count

  int current = seed;
  visited[current] = true;
  int T_accum = 1;  // number of accumulated sites

  // Add first site
  int total_ind = 0;
  int total_sp = 0;
  for (int sp = 0; sp < n_species; sp++) {
    cumulative[sp] += species_mat(current, sp);
    total_ind += species_mat(current, sp);
    if (species_mat(current, sp) > 0) incidence[sp]++;
    if (cumulative[sp] > 0) total_sp++;
  }

  richness[0] = total_sp;
  individuals[0] = total_ind;
  if (coverage_type == 1) {
    coverage[0] = calc_chiu_coverage_internal(cumulative, incidence, T_accum);
  } else {
    coverage[0] = calc_coverage_internal(cumulative);
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
    T_accum++;

    // Accumulate
    total_sp = 0;
    for (int sp = 0; sp < n_species; sp++) {
      cumulative[sp] += species_mat(current, sp);
      total_ind += species_mat(current, sp);
      if (species_mat(current, sp) > 0) incidence[sp]++;
      if (cumulative[sp] > 0) total_sp++;
    }

    richness[step] = total_sp;
    individuals[step] = total_ind;
    if (coverage_type == 1) {
      coverage[step] = calc_chiu_coverage_internal(cumulative, incidence, T_accum);
    } else {
      coverage[step] = calc_coverage_internal(cumulative);
    }
  }

  return List::create(
    Named("richness") = richness,
    Named("individuals") = individuals,
    Named("coverage") = coverage
  );
}


// Worker for parallel coverage accumulation
struct CoverageKnnWorker : public Worker {
  const RMatrix<int> species_mat;
  const RMatrix<double> dist_mat;
  const RVector<int> seeds;

  RMatrix<int> richness;
  RMatrix<int> individuals;
  RMatrix<double> coverage;
  int coverage_type;  // 0 = Chao-Jost, 1 = Chiu

  CoverageKnnWorker(const IntegerMatrix& sp, const NumericMatrix& dm,
                    const IntegerVector& s,
                    IntegerMatrix& r, IntegerMatrix& ind, NumericMatrix& c,
                    int cov_type = 0)
    : species_mat(sp), dist_mat(dm), seeds(s),
      richness(r), individuals(ind), coverage(c),
      coverage_type(cov_type) {}

  void operator()(std::size_t begin, std::size_t end) {
    int n_sites = species_mat.nrow();
    int n_species = species_mat.ncol();

    for (std::size_t s = begin; s < end; s++) {
      std::vector<bool> visited(n_sites, false);
      std::vector<int> cumulative(n_species, 0);
      std::vector<int> incidence(n_species, 0);

      int current = seeds[s];
      visited[current] = true;
      int T_accum = 1;

      int total_ind = 0;
      int total_sp = 0;
      for (int sp = 0; sp < n_species; sp++) {
        cumulative[sp] += species_mat(current, sp);
        total_ind += species_mat(current, sp);
        if (species_mat(current, sp) > 0) incidence[sp]++;
        if (cumulative[sp] > 0) total_sp++;
      }

      richness(s, 0) = total_sp;
      individuals(s, 0) = total_ind;
      if (coverage_type == 1) {
        coverage(s, 0) = calc_chiu_coverage_internal(cumulative, incidence, T_accum);
      } else {
        coverage(s, 0) = calc_coverage_internal(cumulative);
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
        T_accum++;

        total_sp = 0;
        for (int sp = 0; sp < n_species; sp++) {
          cumulative[sp] += species_mat(current, sp);
          total_ind += species_mat(current, sp);
          if (species_mat(current, sp) > 0) incidence[sp]++;
          if (cumulative[sp] > 0) total_sp++;
        }

        richness(s, step) = total_sp;
        individuals(s, step) = total_ind;
        if (coverage_type == 1) {
          coverage(s, step) = calc_chiu_coverage_internal(cumulative, incidence, T_accum);
        } else {
          coverage(s, step) = calc_coverage_internal(cumulative);
        }
      }
    }
  }
};


// coverage_type: 0 = Chao-Jost (individual-based), 1 = Chiu (sample-based)
// [[Rcpp::export]]
List cpp_knn_coverage_parallel(IntegerMatrix species_mat,
                                NumericMatrix dist_mat,
                                int n_seeds,
                                int n_cores = 1,
                                bool progress = false,
                                int coverage_type = 0) {
  int n_sites = species_mat.nrow();

  IntegerVector seeds = Rcpp::sample(n_sites, n_seeds, true) - 1;

  IntegerMatrix richness(n_seeds, n_sites);
  IntegerMatrix individuals(n_seeds, n_sites);
  NumericMatrix coverage(n_seeds, n_sites);

  if (n_cores > 1) {
    CoverageKnnWorker worker(species_mat, dist_mat, seeds,
                             richness, individuals, coverage,
                             coverage_type);
    parallelFor(0, n_seeds, worker);
  } else {
    for (int s = 0; s < n_seeds; s++) {
      List single = cpp_knn_coverage_single(species_mat, dist_mat, seeds[s],
                                            coverage_type);
      IntegerVector r = single["richness"];
      IntegerVector ind = single["individuals"];
      NumericVector c = single["coverage"];

      for (int step = 0; step < n_sites; step++) {
        richness(s, step) = r[step];
        individuals(s, step) = ind[step];
        coverage(s, step) = c[step];
      }
    }
  }

  return List::create(
    Named("richness") = richness,
    Named("individuals") = individuals,
    Named("coverage") = coverage
  );
}


// [[Rcpp::export]]
NumericVector interpolate_at_coverage(NumericVector richness,
                                       NumericVector coverage,
                                       NumericVector targets) {
  int n = richness.size();
  int n_targets = targets.size();
  NumericVector result(n_targets);

  for (int t = 0; t < n_targets; t++) {
    double target = targets[t];

    if (target > coverage[n-1]) {
      // Target above max coverage
      result[t] = NA_REAL;
      continue;
    }

    if (target <= coverage[0]) {
      result[t] = richness[0];
      continue;
    }

    // Find bracketing indices
    for (int i = 1; i < n; i++) {
      if (coverage[i] >= target) {
        // Linear interpolation
        double c0 = coverage[i-1];
        double c1 = coverage[i];
        double r0 = richness[i-1];
        double r1 = richness[i];

        if (c1 == c0) {
          result[t] = r1;
        } else {
          result[t] = r0 + (target - c0) * (r1 - r0) / (c1 - c0);
        }
        break;
      }
    }
  }

  return result;
}
