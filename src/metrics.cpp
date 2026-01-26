// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <set>
using namespace Rcpp;
using namespace RcppParallel;


// Worker struct for parallel kNN metrics (each site is its own seed)
struct KnnMetricsWorker : public Worker {
  // Inputs (read-only)
  const RMatrix<int> species_pa;
  const RMatrix<double> dist_mat;

  // Output
  RMatrix<int> curves;

  // Constructor
  KnnMetricsWorker(const IntegerMatrix& species_pa_,
                   const NumericMatrix& dist_mat_,
                   IntegerMatrix& curves_)
    : species_pa(species_pa_), dist_mat(dist_mat_), curves(curves_) {}

  // Parallel operator - each site is its own seed
  void operator()(std::size_t begin, std::size_t end) {
    int n_sites = species_pa.nrow();
    int n_species = species_pa.ncol();

    for (std::size_t seed = begin; seed < end; seed++) {
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
      curves(seed, 0) = species_seen.size();

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
        curves(seed, step) = species_seen.size();
      }
    }
  }
};


//' Parallel kNN Metrics Accumulation
//'
//' Run kNN accumulation from each site as its own starting point.
//' Returns one curve per site for extracting per-site metrics.
//'
//' @param species_pa Integer matrix (sites x species)
//' @param dist_mat Numeric matrix of pairwise distances
//' @param n_cores Number of cores to use
//' @param progress Show progress (currently ignored in C++)
//' @return Integer matrix (n_sites x n_sites) of accumulation curves
//'
// [[Rcpp::export]]
IntegerMatrix cpp_knn_metrics_parallel(IntegerMatrix species_pa,
                                        NumericMatrix dist_mat,
                                        int n_cores = 1,
                                        bool progress = false) {
  int n_sites = species_pa.nrow();
  IntegerMatrix curves(n_sites, n_sites);

  if (n_cores > 1) {
    // Create worker and run in parallel
    KnnMetricsWorker worker(species_pa, dist_mat, curves);
    parallelFor(0, n_sites, worker);
  } else {
    // Sequential execution
    KnnMetricsWorker worker(species_pa, dist_mat, curves);
    worker(0, n_sites);
  }

  return curves;
}


// Worker struct for parallel kNCN metrics
struct KncnMetricsWorker : public Worker {
  // Inputs (read-only)
  const RMatrix<int> species_pa;
  const RVector<double> x_coords;
  const RVector<double> y_coords;

  // Output
  RMatrix<int> curves;

  // Constructor
  KncnMetricsWorker(const IntegerMatrix& species_pa_,
                    const NumericVector& x_,
                    const NumericVector& y_,
                    IntegerMatrix& curves_)
    : species_pa(species_pa_), x_coords(x_), y_coords(y_), curves(curves_) {}

  // Parallel operator
  void operator()(std::size_t begin, std::size_t end) {
    int n_sites = species_pa.nrow();
    int n_species = species_pa.ncol();

    for (std::size_t seed = begin; seed < end; seed++) {
      std::vector<bool> visited(n_sites, false);
      std::set<int> species_seen;

      // Centroid tracking
      double cx = x_coords[seed];
      double cy = y_coords[seed];
      int n_visited = 1;

      visited[seed] = true;

      // Record species at starting site
      for (int sp = 0; sp < n_species; sp++) {
        if (species_pa(seed, sp) > 0) {
          species_seen.insert(sp);
        }
      }
      curves(seed, 0) = species_seen.size();

      // Visit remaining sites
      for (int step = 1; step < n_sites; step++) {
        // Find nearest unvisited to centroid
        double min_dist = R_PosInf;
        int nearest = -1;

        for (int j = 0; j < n_sites; j++) {
          if (!visited[j]) {
            double dx = x_coords[j] - cx;
            double dy = y_coords[j] - cy;
            double d = std::sqrt(dx * dx + dy * dy);
            if (d < min_dist) {
              min_dist = d;
              nearest = j;
            }
          }
        }

        // Update centroid
        cx = (cx * n_visited + x_coords[nearest]) / (n_visited + 1);
        cy = (cy * n_visited + y_coords[nearest]) / (n_visited + 1);
        n_visited++;

        visited[nearest] = true;

        // Add new species
        for (int sp = 0; sp < n_species; sp++) {
          if (species_pa(nearest, sp) > 0) {
            species_seen.insert(sp);
          }
        }
        curves(seed, step) = species_seen.size();
      }
    }
  }
};


//' Parallel kNCN Metrics Accumulation
//'
//' Run kNCN accumulation from each site as its own starting point.
//'
//' @param species_pa Integer matrix (sites x species)
//' @param x_coords Numeric vector of x coordinates
//' @param y_coords Numeric vector of y coordinates
//' @param n_cores Number of cores to use
//' @param progress Show progress (currently ignored in C++)
//' @return Integer matrix (n_sites x n_sites) of accumulation curves
//'
// [[Rcpp::export]]
IntegerMatrix cpp_kncn_metrics_parallel(IntegerMatrix species_pa,
                                         NumericVector x_coords,
                                         NumericVector y_coords,
                                         int n_cores = 1,
                                         bool progress = false) {
  int n_sites = species_pa.nrow();
  IntegerMatrix curves(n_sites, n_sites);

  if (n_cores > 1) {
    KncnMetricsWorker worker(species_pa, x_coords, y_coords, curves);
    parallelFor(0, n_sites, worker);
  } else {
    KncnMetricsWorker worker(species_pa, x_coords, y_coords, curves);
    worker(0, n_sites);
  }

  return curves;
}
