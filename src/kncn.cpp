#include <Rcpp.h>
#include <vector>
#include <set>
using namespace Rcpp;

//' Single kNCN Accumulation Curve
//'
//' k-Nearest Centroid Neighbor: after each step, recalculate the centroid
//' of all visited sites, then find the unvisited site nearest to that centroid.
//'
//' @param species_pa Integer matrix (sites × species), presence/absence
//' @param dist_mat Numeric matrix of pairwise distances (unused, compute from coords)
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
    for (int v : visited_sites) {
      cx += x[v];
      cy += y[v];
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


//' Parallel kNCN Accumulation (multiple seeds)
//'
//' @param species_pa Integer matrix (sites × species)
//' @param x Numeric vector of x coordinates
//' @param y Numeric vector of y coordinates
//' @param n_seeds Number of random starting points
//' @return Integer matrix (n_seeds × n_sites) of accumulation curves
//'
// [[Rcpp::export]]
IntegerMatrix cpp_kncn_parallel(IntegerMatrix species_pa,
                                NumericVector x, NumericVector y,
                                int n_seeds) {
  int n_sites = species_pa.nrow();
  IntegerMatrix curves(n_seeds, n_sites);

  // Generate random seeds
  IntegerVector seeds = Rcpp::sample(n_sites, n_seeds, true) - 1;

  // TODO: Replace with RcppParallel Worker
  for (int s = 0; s < n_seeds; s++) {
    IntegerVector curve = cpp_kncn_single(species_pa, x, y, seeds[s]);
    curves(s, _) = curve;
  }

  return curves;
}
