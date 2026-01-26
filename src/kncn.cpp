// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <set>
#include <cmath>
using namespace Rcpp;
using namespace RcppParallel;


//' Single kNCN Accumulation Curve
//'
//' k-Nearest Centroid Neighbor: after each step, recalculate the centroid
//' of all visited sites, then find the unvisited site nearest to that centroid.
//'
//' @param species_pa Integer matrix (sites x species), presence/absence
//' @param x Numeric vector of x coordinates
//' @param y Numeric vector of y coordinates
//' @param seed Starting site index (0-based)
//' @return Integer vector of cumulative species counts
//'
// [[Rcpp::export]]
IntegerVector cpp_kncn_single(IntegerMatrix species_pa,
                              NumericVector x, NumericVector y,
                              int seed) {
  int n_sites = species_pa.nrow();
  int n_species = species_pa.ncol();

  IntegerVector curve(n_sites);
  std::vector<bool> visited(n_sites, false);
  std::set<int> species_seen;
  std::vector<int> visited_sites;

  // Start at seed
  int current = seed;
  visited[current] = true;
  visited_sites.push_back(current);

  // Record species at starting site
  for (int sp = 0; sp < n_species; sp++) {
    if (species_pa(current, sp) > 0) {
      species_seen.insert(sp);
    }
  }
  curve[0] = species_seen.size();

  // Visit remaining sites
  for (int step = 1; step < n_sites; step++) {
    // Calculate centroid of visited sites
    double cx = 0.0, cy = 0.0;
    for (size_t v = 0; v < visited_sites.size(); v++) {
      cx += x[visited_sites[v]];
      cy += y[visited_sites[v]];
    }
    cx /= visited_sites.size();
    cy /= visited_sites.size();

    // Find nearest unvisited to centroid
    double min_dist = R_PosInf;
    int nearest = -1;

    for (int j = 0; j < n_sites; j++) {
      if (!visited[j]) {
        double dx = x[j] - cx;
        double dy = y[j] - cy;
        double d = std::sqrt(dx * dx + dy * dy);
        if (d < min_dist) {
          min_dist = d;
          nearest = j;
        }
      }
    }

    // Move to nearest
    current = nearest;
    visited[current] = true;
    visited_sites.push_back(current);

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


// Worker struct for parallel kNCN
struct KncnWorker : public Worker {
  // Inputs (read-only)
  const RMatrix<int> species_pa;
  const RVector<double> x;
  const RVector<double> y;
  const RVector<int> seeds;

  // Output
  RMatrix<int> curves;

  // Constructor
  KncnWorker(const IntegerMatrix& species_pa_,
             const NumericVector& x_,
             const NumericVector& y_,
             const IntegerVector& seeds_,
             IntegerMatrix& curves_)
    : species_pa(species_pa_), x(x_), y(y_),
      seeds(seeds_), curves(curves_) {}

  // Parallel operator
  void operator()(std::size_t begin, std::size_t end) {
    int n_sites = species_pa.nrow();
    int n_species = species_pa.ncol();

    for (std::size_t s = begin; s < end; s++) {
      std::vector<bool> visited(n_sites, false);
      std::set<int> species_seen;
      std::vector<int> visited_sites;

      int current = seeds[s];
      visited[current] = true;
      visited_sites.push_back(current);

      // Record species at starting site
      for (int sp = 0; sp < n_species; sp++) {
        if (species_pa(current, sp) > 0) {
          species_seen.insert(sp);
        }
      }
      curves(s, 0) = species_seen.size();

      // Visit remaining sites
      for (int step = 1; step < n_sites; step++) {
        // Calculate centroid
        double cx = 0.0, cy = 0.0;
        for (size_t v = 0; v < visited_sites.size(); v++) {
          cx += x[visited_sites[v]];
          cy += y[visited_sites[v]];
        }
        cx /= visited_sites.size();
        cy /= visited_sites.size();

        // Find nearest unvisited to centroid
        double min_dist = R_PosInf;
        int nearest = -1;

        for (int j = 0; j < n_sites; j++) {
          if (!visited[j]) {
            double dx = x[j] - cx;
            double dy = y[j] - cy;
            double d = std::sqrt(dx * dx + dy * dy);
            if (d < min_dist) {
              min_dist = d;
              nearest = j;
            }
          }
        }

        current = nearest;
        visited[current] = true;
        visited_sites.push_back(current);

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


//' Parallel kNCN Accumulation
//'
//' Run kNCN accumulation from multiple random starting points in parallel.
//'
//' @param species_pa Integer matrix (sites x species)
//' @param x Numeric vector of x coordinates
//' @param y Numeric vector of y coordinates
//' @param n_seeds Number of random starting points
//' @param n_cores Number of cores to use
//' @param progress Show progress (currently ignored)
//' @return Integer matrix (n_seeds x n_sites) of accumulation curves
//'
// [[Rcpp::export]]
IntegerMatrix cpp_kncn_parallel(IntegerMatrix species_pa,
                                NumericVector x, NumericVector y,
                                int n_seeds,
                                int n_cores = 1,
                                bool progress = false) {
  int n_sites = species_pa.nrow();
  IntegerMatrix curves(n_seeds, n_sites);

  // Generate random seeds (0-based indices)
  IntegerVector seeds = Rcpp::sample(n_sites, n_seeds, true) - 1;

  if (n_cores > 1) {
    KncnWorker worker(species_pa, x, y, seeds, curves);
    parallelFor(0, n_seeds, worker);
  } else {
    for (int s = 0; s < n_seeds; s++) {
      IntegerVector curve = cpp_kncn_single(species_pa, x, y, seeds[s]);
      curves(s, _) = curve;
    }
  }

  return curves;
}
