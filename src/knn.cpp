// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <set>
using namespace Rcpp;
using namespace RcppParallel;


//' Single kNN Accumulation Curve
//'
//' Traverse sites in nearest-neighbor order from a starting seed,
//' accumulating species counts.
//'
//' @param species_pa Integer matrix (sites x species), presence/absence
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


// Worker struct for parallel kNN
struct KnnWorker : public Worker {
  // Inputs (read-only)
  const RMatrix<int> species_pa;
  const RMatrix<double> dist_mat;
  const RVector<int> seeds;

  // Output
  RMatrix<int> curves;

  // Constructor
  KnnWorker(const IntegerMatrix& species_pa_,
            const NumericMatrix& dist_mat_,
            const IntegerVector& seeds_,
            IntegerMatrix& curves_)
    : species_pa(species_pa_), dist_mat(dist_mat_),
      seeds(seeds_), curves(curves_) {}

  // Parallel operator
  void operator()(std::size_t begin, std::size_t end) {
    int n_sites = species_pa.nrow();
    int n_species = species_pa.ncol();

    for (std::size_t s = begin; s < end; s++) {
      std::vector<bool> visited(n_sites, false);
      std::set<int> species_seen;

      int current = seeds[s];
      visited[current] = true;

      // Record species at starting site
      for (int sp = 0; sp < n_species; sp++) {
        if (species_pa(current, sp) > 0) {
          species_seen.insert(sp);
        }
      }
      curves(s, 0) = species_seen.size();

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

        current = nearest;
        visited[current] = true;

        // Add new species
        for (int sp = 0; sp < n_species; sp++) {
          if (species_pa(current, sp) > 0) {
            species_seen.insert(sp);
          }
        }
        curves(s, step) = species_seen.size();
      }
    }
  }
};


//' Parallel kNN Accumulation
//'
//' Run kNN accumulation from multiple random starting points in parallel.
//'
//' @param species_pa Integer matrix (sites x species)
//' @param dist_mat Numeric matrix of pairwise distances
//' @param n_seeds Number of random starting points
//' @param n_cores Number of cores to use
//' @param progress Show progress (currently ignored in C++)
//' @return Integer matrix (n_seeds x n_sites) of accumulation curves
//'
// [[Rcpp::export]]
IntegerMatrix cpp_knn_parallel(IntegerMatrix species_pa,
                               NumericMatrix dist_mat,
                               int n_seeds,
                               int n_cores = 1,
                               bool progress = false) {
  int n_sites = species_pa.nrow();
  IntegerMatrix curves(n_seeds, n_sites);

  // Generate random seeds (0-based indices)
  IntegerVector seeds = Rcpp::sample(n_sites, n_seeds, true) - 1;

  // Set number of threads
  if (n_cores > 1) {
    // Create worker and run in parallel
    KnnWorker worker(species_pa, dist_mat, seeds, curves);
    parallelFor(0, n_seeds, worker);
  } else {
    // Sequential execution
    for (int s = 0; s < n_seeds; s++) {
      IntegerVector curve = cpp_knn_single(species_pa, dist_mat, seeds[s]);
      curves(s, _) = curve;
    }
  }

  return curves;
}


//' Parallel kNN Accumulation with Explicit Seeds
//'
//' Run kNN accumulation from specified starting points in parallel.
//'
//' @param species_pa Integer matrix (sites x species)
//' @param dist_mat Numeric matrix of pairwise distances
//' @param seeds Integer vector of starting point indices (0-based)
//' @param n_cores Number of cores to use
//' @param progress Show progress (currently ignored in C++)
//' @return Integer matrix (n_seeds x n_sites) of accumulation curves
//'
// [[Rcpp::export]]
IntegerMatrix cpp_knn_parallel_seeds(IntegerMatrix species_pa,
                                      NumericMatrix dist_mat,
                                      IntegerVector seeds,
                                      int n_cores = 1,
                                      bool progress = false) {
  int n_sites = species_pa.nrow();
  int n_seeds = seeds.size();
  IntegerMatrix curves(n_seeds, n_sites);

  if (n_cores > 1) {
    KnnWorker worker(species_pa, dist_mat, seeds, curves);
    parallelFor(0, n_seeds, worker);
  } else {
    for (int s = 0; s < n_seeds; s++) {
      IntegerVector curve = cpp_knn_single(species_pa, dist_mat, seeds[s]);
      curves(s, _) = curve;
    }
  }

  return curves;
}
