#include <Rcpp.h>
#include <vector>
#include <set>
using namespace Rcpp;

//' Single kNN Accumulation Curve
//'
//' Traverse sites in nearest-neighbor order from a starting seed,
//' accumulating species counts.
//'
//' @param species_pa Integer matrix (sites × species), presence/absence
//' @param dist_mat Numeric matrix of pairwise distances
//' @param seed Starting site index (0-based)
//' @return Integer vector of cumulative species counts
//'
// [[Rcpp::export]]
IntegerVector cpp_knn_single(IntegerMatrix species_pa, NumericMatrix dist_mat, int seed) {
  int n_sites = species_pa.nrow();
  int n_species = species_pa.ncol();

  IntegerVector curve(n_sites);
  std::vector<bool> visited(n_sites, false);
  std::set<int> species_seen;

  int current = seed;
  visited[current] = true;

  // Record species at starting site
  for (int sp = 0; sp < n_species; sp++) {
    if (species_pa(current, sp) > 0) {
      species_seen.insert(sp);
    }
  }
  curve[0] = species_seen.size();

  // Visit remaining sites
 for (int step = 1; step < n_sites; step++) {
    // Find nearest unvisited
    double min_dist = R_PosInf;
    int nearest = -1;

    for (int j = 0; j < n_sites; j++) {
      if (!visited[j] && dist_mat(current, j) < min_dist) {
        min_dist = dist_mat(current, j);
        nearest = j;
      }
    }

    // Move to nearest
    current = nearest;
    visited[current] = true;

    // Add new species
    for (int sp = 0; sp < n_species; sp++) {
      if (species_pa(current, sp) > 0) {
        species_seen.insert(sp);
      }
    }
    curve[step] = species_seen.size();
  }

  return curve;
}


//' Parallel kNN Accumulation (multiple seeds)
//'
//' Run kNN accumulation from multiple random starting points.
//' TODO: Implement with RcppParallel for true parallelism.
//'
//' @param species_pa Integer matrix (sites × species)
//' @param dist_mat Numeric matrix of pairwise distances
//' @param n_seeds Number of random starting points
//' @return Integer matrix (n_seeds × n_sites) of accumulation curves
//'
// [[Rcpp::export]]
IntegerMatrix cpp_knn_parallel(IntegerMatrix species_pa, NumericMatrix dist_mat, int n_seeds) {
  int n_sites = species_pa.nrow();
  IntegerMatrix curves(n_seeds, n_sites);

  // Generate random seeds
  IntegerVector seeds = Rcpp::sample(n_sites, n_seeds, true) - 1; // 0-based

  // TODO: Replace with RcppParallel Worker
  for (int s = 0; s < n_seeds; s++) {
    IntegerVector curve = cpp_knn_single(species_pa, dist_mat, seeds[s]);
    curves(s, _) = curve;
  }

  return curves;
}
