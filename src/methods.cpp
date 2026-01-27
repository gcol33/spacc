// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <set>
#include <cmath>
#include <algorithm>
#include <random>
using namespace Rcpp;
using namespace RcppParallel;


// ============================================================================
// COLLECTOR METHOD - sites in data order
// ============================================================================

// [[Rcpp::export]]
IntegerVector cpp_collector_single(IntegerMatrix species_pa) {
  int n_sites = species_pa.nrow();
  int n_species = species_pa.ncol();

  IntegerVector curve(n_sites);
  std::set<int> species_seen;

  for (int i = 0; i < n_sites; i++) {
    for (int sp = 0; sp < n_species; sp++) {
      if (species_pa(i, sp) > 0) {
        species_seen.insert(sp);
      }
    }
    curve[i] = species_seen.size();
  }

  return curve;
}


// ============================================================================
// GAUSSIAN-WEIGHTED SELECTION
// ============================================================================

// [[Rcpp::export]]
IntegerVector cpp_gaussian_single(IntegerMatrix species_pa,
                                  NumericMatrix dist_mat,
                                  int seed,
                                  double sigma) {
  int n_sites = species_pa.nrow();
  int n_species = species_pa.ncol();

  IntegerVector curve(n_sites);
  std::vector<bool> visited(n_sites, false);
  std::set<int> species_seen;

  // Random generator
  std::random_device rd;
  std::mt19937 gen(rd());

  int current = seed;
  visited[current] = true;

  // Record species at starting site
  for (int sp = 0; sp < n_species; sp++) {
    if (species_pa(current, sp) > 0) {
      species_seen.insert(sp);
    }
  }
  curve[0] = species_seen.size();

  double sigma2 = 2.0 * sigma * sigma;

  for (int step = 1; step < n_sites; step++) {
    // Compute Gaussian weights for unvisited sites
    std::vector<double> weights;
    std::vector<int> candidates;

    for (int j = 0; j < n_sites; j++) {
      if (!visited[j]) {
        double d = dist_mat(current, j);
        double w = std::exp(-(d * d) / sigma2);
        weights.push_back(w);
        candidates.push_back(j);
      }
    }

    // Normalize and sample
    double sum_w = 0.0;
    for (size_t i = 0; i < weights.size(); i++) sum_w += weights[i];
    for (size_t i = 0; i < weights.size(); i++) weights[i] /= sum_w;

    std::discrete_distribution<> dist(weights.begin(), weights.end());
    int next_idx = dist(gen);
    current = candidates[next_idx];

    visited[current] = true;

    for (int sp = 0; sp < n_species; sp++) {
      if (species_pa(current, sp) > 0) {
        species_seen.insert(sp);
      }
    }
    curve[step] = species_seen.size();
  }

  return curve;
}


// Parallel worker for Gaussian
struct GaussianWorker : public Worker {
  const RMatrix<int> species_pa;
  const RMatrix<double> dist_mat;
  const RVector<int> seeds;
  const double sigma;
  RMatrix<int> curves;

  GaussianWorker(const IntegerMatrix& sp, const NumericMatrix& dm,
                 const IntegerVector& s, double sig, IntegerMatrix& c)
    : species_pa(sp), dist_mat(dm), seeds(s), sigma(sig), curves(c) {}

  void operator()(std::size_t begin, std::size_t end) {
    int n_sites = species_pa.nrow();
    int n_species = species_pa.ncol();
    double sigma2 = 2.0 * sigma * sigma;

    std::random_device rd;
    std::mt19937 gen(rd());

    for (std::size_t s = begin; s < end; s++) {
      std::vector<bool> visited(n_sites, false);
      std::set<int> species_seen;

      int current = seeds[s];
      visited[current] = true;

      for (int sp = 0; sp < n_species; sp++) {
        if (species_pa(current, sp) > 0) species_seen.insert(sp);
      }
      curves(s, 0) = species_seen.size();

      for (int step = 1; step < n_sites; step++) {
        std::vector<double> weights;
        std::vector<int> candidates;

        for (int j = 0; j < n_sites; j++) {
          if (!visited[j]) {
            double d = dist_mat(current, j);
            weights.push_back(std::exp(-(d * d) / sigma2));
            candidates.push_back(j);
          }
        }

        double sum_w = 0.0;
        for (size_t i = 0; i < weights.size(); i++) sum_w += weights[i];
        for (size_t i = 0; i < weights.size(); i++) weights[i] /= sum_w;

        std::discrete_distribution<> dist(weights.begin(), weights.end());
        current = candidates[dist(gen)];
        visited[current] = true;

        for (int sp = 0; sp < n_species; sp++) {
          if (species_pa(current, sp) > 0) species_seen.insert(sp);
        }
        curves(s, step) = species_seen.size();
      }
    }
  }
};


// [[Rcpp::export]]
IntegerMatrix cpp_gaussian_parallel(IntegerMatrix species_pa,
                                    NumericMatrix dist_mat,
                                    int n_seeds,
                                    double sigma,
                                    int n_cores = 1,
                                    bool progress = false) {
  int n_sites = species_pa.nrow();
  IntegerMatrix curves(n_seeds, n_sites);
  IntegerVector seeds = Rcpp::sample(n_sites, n_seeds, true) - 1;

  if (n_cores > 1) {
    GaussianWorker worker(species_pa, dist_mat, seeds, sigma, curves);
    parallelFor(0, n_seeds, worker);
  } else {
    for (int s = 0; s < n_seeds; s++) {
      IntegerVector curve = cpp_gaussian_single(species_pa, dist_mat, seeds[s], sigma);
      curves(s, _) = curve;
    }
  }

  return curves;
}


// ============================================================================
// EXPANDING RADIUS (WAVEFRONT) METHOD
// ============================================================================

// [[Rcpp::export]]
List cpp_wavefront_single(IntegerMatrix species_pa,
                          NumericMatrix dist_mat,
                          int seed,
                          double r0,
                          double dr,
                          int n_steps) {
  int n_sites = species_pa.nrow();
  int n_species = species_pa.ncol();

  IntegerVector species_count(n_steps);
  NumericVector radius_vals(n_steps);
  IntegerVector sites_included(n_steps);

  for (int step = 0; step < n_steps; step++) {
    double radius = r0 + step * dr;
    radius_vals[step] = radius;

    std::set<int> species_seen;
    int n_included = 0;

    for (int i = 0; i < n_sites; i++) {
      if (dist_mat(seed, i) <= radius) {
        n_included++;
        for (int sp = 0; sp < n_species; sp++) {
          if (species_pa(i, sp) > 0) species_seen.insert(sp);
        }
      }
    }

    species_count[step] = species_seen.size();
    sites_included[step] = n_included;
  }

  return List::create(
    Named("species") = species_count,
    Named("radius") = radius_vals,
    Named("n_sites") = sites_included
  );
}


// [[Rcpp::export]]
List cpp_wavefront_parallel(IntegerMatrix species_pa,
                            NumericMatrix dist_mat,
                            int n_seeds,
                            double r0,
                            double dr,
                            int n_steps,
                            int n_cores = 1,
                            bool progress = false) {
  int n_sites = species_pa.nrow();
  IntegerMatrix curves(n_seeds, n_steps);
  IntegerMatrix sites_mat(n_seeds, n_steps);
  IntegerVector seeds = Rcpp::sample(n_sites, n_seeds, true) - 1;

  NumericVector radius_vals(n_steps);
  for (int step = 0; step < n_steps; step++) {
    radius_vals[step] = r0 + step * dr;
  }

  for (int s = 0; s < n_seeds; s++) {
    List result = cpp_wavefront_single(species_pa, dist_mat, seeds[s], r0, dr, n_steps);
    IntegerVector sp_count = result["species"];
    IntegerVector n_sites_inc = result["n_sites"];
    curves(s, _) = sp_count;
    sites_mat(s, _) = n_sites_inc;
  }

  return List::create(
    Named("curves") = curves,
    Named("radius") = radius_vals,
    Named("sites_included") = sites_mat
  );
}


// ============================================================================
// SIMPLE RADIUS EXPANSION (by distance order)
// ============================================================================

// [[Rcpp::export]]
IntegerVector cpp_radius_single(IntegerMatrix species_pa,
                                NumericMatrix dist_mat,
                                int seed) {
  int n_sites = species_pa.nrow();
  int n_species = species_pa.ncol();

  // Get distances from seed and sort
  std::vector<std::pair<double, int>> dist_idx(n_sites);
  for (int i = 0; i < n_sites; i++) {
    dist_idx[i] = std::make_pair(dist_mat(seed, i), i);
  }
  std::sort(dist_idx.begin(), dist_idx.end());

  // Accumulate in order of distance
  IntegerVector curve(n_sites);
  std::set<int> species_seen;

  for (int step = 0; step < n_sites; step++) {
    int site = dist_idx[step].second;
    for (int sp = 0; sp < n_species; sp++) {
      if (species_pa(site, sp) > 0) {
        species_seen.insert(sp);
      }
    }
    curve[step] = species_seen.size();
  }

  return curve;
}


// Parallel worker for radius expansion
struct RadiusWorker : public Worker {
  const RMatrix<int> species_pa;
  const RMatrix<double> dist_mat;
  const RVector<int> seeds;
  RMatrix<int> curves;

  RadiusWorker(const IntegerMatrix& sp, const NumericMatrix& dm,
               const IntegerVector& s, IntegerMatrix& c)
    : species_pa(sp), dist_mat(dm), seeds(s), curves(c) {}

  void operator()(std::size_t begin, std::size_t end) {
    int n_sites = species_pa.nrow();
    int n_species = species_pa.ncol();

    for (std::size_t s = begin; s < end; s++) {
      int seed = seeds[s];

      // Sort by distance from seed
      std::vector<std::pair<double, int>> dist_idx(n_sites);
      for (int i = 0; i < n_sites; i++) {
        dist_idx[i] = std::make_pair(dist_mat(seed, i), i);
      }
      std::sort(dist_idx.begin(), dist_idx.end());

      std::set<int> species_seen;
      for (int step = 0; step < n_sites; step++) {
        int site = dist_idx[step].second;
        for (int sp = 0; sp < n_species; sp++) {
          if (species_pa(site, sp) > 0) species_seen.insert(sp);
        }
        curves(s, step) = species_seen.size();
      }
    }
  }
};


// [[Rcpp::export]]
IntegerMatrix cpp_radius_parallel(IntegerMatrix species_pa,
                                  NumericMatrix dist_mat,
                                  int n_seeds,
                                  int n_cores = 1,
                                  bool progress = false) {
  int n_sites = species_pa.nrow();
  IntegerMatrix curves(n_seeds, n_sites);
  IntegerVector seeds = Rcpp::sample(n_sites, n_seeds, true) - 1;

  if (n_cores > 1) {
    RadiusWorker worker(species_pa, dist_mat, seeds, curves);
    parallelFor(0, n_seeds, worker);
  } else {
    for (int s = 0; s < n_seeds; s++) {
      IntegerVector curve = cpp_radius_single(species_pa, dist_mat, seeds[s]);
      curves(s, _) = curve;
    }
  }

  return curves;
}


// ============================================================================
// DIRECTIONAL CONE METHOD
// ============================================================================

// [[Rcpp::export]]
IntegerVector cpp_cone_single(IntegerMatrix species_pa,
                              NumericVector x, NumericVector y,
                              int seed,
                              double angle,
                              double width) {
  int n_sites = species_pa.nrow();
  int n_species = species_pa.ncol();

  double seed_x = x[seed];
  double seed_y = y[seed];

  // Compute distance and angle from seed for each site
  std::vector<std::tuple<double, double, int>> sites_info; // dist, angle_diff, idx

  for (int i = 0; i < n_sites; i++) {
    double dx = x[i] - seed_x;
    double dy = y[i] - seed_y;
    double dist = std::sqrt(dx * dx + dy * dy);
    double site_angle = std::atan2(dy, dx);

    // Angular difference (handle wraparound)
    double angle_diff = std::abs(site_angle - angle);
    if (angle_diff > M_PI) angle_diff = 2.0 * M_PI - angle_diff;

    sites_info.push_back(std::make_tuple(dist, angle_diff, i));
  }

  // Sort by distance, but only include sites within cone
  std::vector<std::pair<double, int>> in_cone;
  std::vector<std::pair<double, int>> out_cone;

  for (const auto& si : sites_info) {
    double dist = std::get<0>(si);
    double angle_diff = std::get<1>(si);
    int idx = std::get<2>(si);

    if (angle_diff <= width || idx == seed) {
      in_cone.push_back(std::make_pair(dist, idx));
    } else {
      out_cone.push_back(std::make_pair(dist, idx));
    }
  }

  std::sort(in_cone.begin(), in_cone.end());
  std::sort(out_cone.begin(), out_cone.end());

  // Accumulate: first sites in cone, then rest
  IntegerVector curve(n_sites);
  std::set<int> species_seen;
  int step = 0;

  for (const auto& p : in_cone) {
    int site = p.second;
    for (int sp = 0; sp < n_species; sp++) {
      if (species_pa(site, sp) > 0) species_seen.insert(sp);
    }
    curve[step++] = species_seen.size();
  }

  for (const auto& p : out_cone) {
    int site = p.second;
    for (int sp = 0; sp < n_species; sp++) {
      if (species_pa(site, sp) > 0) species_seen.insert(sp);
    }
    curve[step++] = species_seen.size();
  }

  return curve;
}


// [[Rcpp::export]]
IntegerMatrix cpp_cone_parallel(IntegerMatrix species_pa,
                                NumericVector x, NumericVector y,
                                int n_seeds,
                                double width = 0.785398, // pi/4 = 45 degrees
                                int n_cores = 1,
                                bool progress = false) {
  int n_sites = species_pa.nrow();
  IntegerMatrix curves(n_seeds, n_sites);
  IntegerVector seeds = Rcpp::sample(n_sites, n_seeds, true) - 1;

  // Random angles for each seed
  NumericVector angles = Rcpp::runif(n_seeds, 0, 2.0 * M_PI);

  // Sequential for now (cone is fast enough)
  for (int s = 0; s < n_seeds; s++) {
    IntegerVector curve = cpp_cone_single(species_pa, x, y, seeds[s], angles[s], width);
    curves(s, _) = curve;
  }

  return curves;
}


// ============================================================================
// DISTANCE-DECAY METHOD
// ============================================================================

// [[Rcpp::export]]
IntegerVector cpp_distance_decay_single(IntegerMatrix species_pa,
                                        NumericMatrix dist_mat,
                                        int seed,
                                        NumericVector breaks) {
  int n_sites = species_pa.nrow();
  int n_species = species_pa.ncol();
  int n_breaks = breaks.size();

  IntegerVector curve(n_breaks);

  for (int b = 0; b < n_breaks; b++) {
    double threshold = breaks[b];
    std::set<int> species_seen;

    for (int i = 0; i < n_sites; i++) {
      if (dist_mat(seed, i) <= threshold) {
        for (int sp = 0; sp < n_species; sp++) {
          if (species_pa(i, sp) > 0) species_seen.insert(sp);
        }
      }
    }
    curve[b] = species_seen.size();
  }

  return curve;
}


// [[Rcpp::export]]
IntegerMatrix cpp_distance_decay_parallel(IntegerMatrix species_pa,
                                          NumericMatrix dist_mat,
                                          int n_seeds,
                                          NumericVector breaks,
                                          int n_cores = 1,
                                          bool progress = false) {
  int n_breaks = breaks.size();
  IntegerMatrix curves(n_seeds, n_breaks);
  int n_sites = species_pa.nrow();
  IntegerVector seeds = Rcpp::sample(n_sites, n_seeds, true) - 1;

  for (int s = 0; s < n_seeds; s++) {
    IntegerVector curve = cpp_distance_decay_single(species_pa, dist_mat, seeds[s], breaks);
    curves(s, _) = curve;
  }

  return curves;
}
