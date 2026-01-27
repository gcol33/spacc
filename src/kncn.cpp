// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <cmath>
#include "kdtree_adapter.h"
#include "balltree.h"
using namespace Rcpp;
using namespace RcppParallel;


// ============================================================================
// EXACT (brute-force) kNCN
// Uses: incremental centroid, squared distances, shrinking candidate list,
//       boolean vector for species tracking.
// ============================================================================

// [[Rcpp::export]]
IntegerVector cpp_kncn_single(IntegerMatrix species_pa,
                              NumericVector x, NumericVector y,
                              int seed) {
  int n_sites = species_pa.nrow();
  int n_species = species_pa.ncol();

  IntegerVector curve(n_sites);
  std::vector<bool> species_seen(n_species, false);
  int richness = 0;

  // Shrinking candidate list: indices of unvisited sites
  std::vector<int> candidates(n_sites);
  for (int i = 0; i < n_sites; i++) candidates[i] = i;

  // Remove seed from candidates
  for (int i = 0; i < n_sites; i++) {
    if (candidates[i] == seed) {
      candidates[i] = candidates.back();
      candidates.pop_back();
      break;
    }
  }

  // Incremental centroid
  double sum_x = x[seed];
  double sum_y = y[seed];
  int n_visited = 1;

  // Record species at seed
  int current = seed;
  for (int sp = 0; sp < n_species; sp++) {
    if (species_pa(current, sp) > 0 && !species_seen[sp]) {
      species_seen[sp] = true;
      richness++;
    }
  }
  curve[0] = richness;

  // Visit remaining sites
  for (int step = 1; step < n_sites; step++) {
    double cx = sum_x / n_visited;
    double cy = sum_y / n_visited;

    // Find nearest unvisited to centroid (squared distance, no sqrt)
    double min_dist_sq = R_PosInf;
    int best_idx = -1; // index into candidates

    int n_cand = static_cast<int>(candidates.size());
    for (int i = 0; i < n_cand; i++) {
      int j = candidates[i];
      double dx = x[j] - cx;
      double dy = y[j] - cy;
      double d_sq = dx * dx + dy * dy;
      if (d_sq < min_dist_sq) {
        min_dist_sq = d_sq;
        best_idx = i;
      }
    }

    // Move to nearest
    current = candidates[best_idx];
    // Swap-remove from candidates
    candidates[best_idx] = candidates.back();
    candidates.pop_back();

    sum_x += x[current];
    sum_y += y[current];
    n_visited++;

    // Add new species
    for (int sp = 0; sp < n_species; sp++) {
      if (species_pa(current, sp) > 0 && !species_seen[sp]) {
        species_seen[sp] = true;
        richness++;
      }
    }
    curve[step] = richness;
  }

  return curve;
}


// Worker struct for parallel exact kNCN
struct KncnWorker : public Worker {
  const RMatrix<int> species_pa;
  const RVector<double> x;
  const RVector<double> y;
  const RVector<int> seeds;
  RMatrix<int> curves;

  KncnWorker(const IntegerMatrix& species_pa_,
             const NumericVector& x_,
             const NumericVector& y_,
             const IntegerVector& seeds_,
             IntegerMatrix& curves_)
    : species_pa(species_pa_), x(x_), y(y_),
      seeds(seeds_), curves(curves_) {}

  void operator()(std::size_t begin, std::size_t end) {
    int n_sites = species_pa.nrow();
    int n_species = species_pa.ncol();

    for (std::size_t s = begin; s < end; s++) {
      std::vector<bool> species_seen(n_species, false);
      int richness = 0;

      // Shrinking candidate list
      std::vector<int> candidates(n_sites);
      for (int i = 0; i < n_sites; i++) candidates[i] = i;

      int current = seeds[s];
      // Remove seed from candidates
      for (int i = 0; i < n_sites; i++) {
        if (candidates[i] == current) {
          candidates[i] = candidates.back();
          candidates.pop_back();
          break;
        }
      }

      double sum_x = x[current];
      double sum_y = y[current];
      int n_visited = 1;

      for (int sp = 0; sp < n_species; sp++) {
        if (species_pa(current, sp) > 0) {
          species_seen[sp] = true;
          richness++;
        }
      }
      curves(s, 0) = richness;

      for (int step = 1; step < n_sites; step++) {
        double cx = sum_x / n_visited;
        double cy = sum_y / n_visited;

        double min_dist_sq = R_PosInf;
        int best_idx = -1;

        int n_cand = static_cast<int>(candidates.size());
        for (int i = 0; i < n_cand; i++) {
          int j = candidates[i];
          double dx = x[j] - cx;
          double dy = y[j] - cy;
          double d_sq = dx * dx + dy * dy;
          if (d_sq < min_dist_sq) {
            min_dist_sq = d_sq;
            best_idx = i;
          }
        }

        current = candidates[best_idx];
        candidates[best_idx] = candidates.back();
        candidates.pop_back();

        sum_x += x[current];
        sum_y += y[current];
        n_visited++;

        for (int sp = 0; sp < n_species; sp++) {
          if (species_pa(current, sp) > 0 && !species_seen[sp]) {
            species_seen[sp] = true;
            richness++;
          }
        }
        curves(s, step) = richness;
      }
    }
  }
};


// [[Rcpp::export]]
IntegerMatrix cpp_kncn_parallel(IntegerMatrix species_pa,
                                NumericVector x, NumericVector y,
                                int n_seeds,
                                int n_cores = 1,
                                bool progress = false) {
  int n_sites = species_pa.nrow();
  IntegerMatrix curves(n_seeds, n_sites);

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


// ============================================================================
// SPATIAL TREE kNCN — k-d tree (Euclidean) or ball tree (haversine)
// Queries tree with centroid each step.
// Uses boolean vector for species tracking + incremental centroid.
// ============================================================================

// Generic kNCN accumulation with a find-nearest functor
template<typename FindNN>
IntegerVector kncn_tree_single_impl(const IntegerMatrix& species_pa,
                                     const std::vector<double>& px,
                                     const std::vector<double>& py,
                                     int seed, FindNN find_nn) {
  int n_sites = species_pa.nrow();
  int n_species = species_pa.ncol();

  IntegerVector curve(n_sites);
  std::vector<bool> visited(n_sites, false);
  std::vector<bool> species_seen(n_species, false);
  int richness = 0;

  double sum_x = px[seed];
  double sum_y = py[seed];
  int n_visited = 1;

  int current = seed;
  visited[current] = true;

  for (int sp = 0; sp < n_species; sp++) {
    if (species_pa(current, sp) > 0) {
      species_seen[sp] = true;
      richness++;
    }
  }
  curve[0] = richness;

  for (int step = 1; step < n_sites; step++) {
    double cx = sum_x / n_visited;
    double cy = sum_y / n_visited;

    int nearest = find_nn(cx, cy, visited);

    current = nearest;
    visited[current] = true;
    sum_x += px[current];
    sum_y += py[current];
    n_visited++;

    for (int sp = 0; sp < n_species; sp++) {
      if (species_pa(current, sp) > 0 && !species_seen[sp]) {
        species_seen[sp] = true;
        richness++;
      }
    }
    curve[step] = richness;
  }

  return curve;
}


// [[Rcpp::export]]
IntegerVector cpp_kncn_kdtree_single(IntegerMatrix species_pa,
                                      NumericVector x, NumericVector y,
                                      int seed,
                                      std::string distance = "euclidean") {
  int n_sites = species_pa.nrow();
  bool use_haversine = (distance == "haversine");

  std::vector<double> vx(n_sites), vy(n_sites);
  for (int i = 0; i < n_sites; i++) { vx[i] = x[i]; vy[i] = y[i]; }

  if (use_haversine) {
    BallTree btree(vx, vy);
    return kncn_tree_single_impl(species_pa, vx, vy, seed,
      [&](double qx, double qy, const std::vector<bool>& vis) {
        return btree.find_nearest_unvisited(qx, qy, vis);
      });
  } else {
    PointCloud2D cloud;
    cloud.pts_x = vx;
    cloud.pts_y = vy;
    KDTree2D* tree = build_kdtree(cloud);
    auto result = kncn_tree_single_impl(species_pa, vx, vy, seed,
      [&](double qx, double qy, const std::vector<bool>& vis) {
        return find_nearest_unvisited(*tree, qx, qy, vis, n_sites);
      });
    delete tree;
    return result;
  }
}


// Parallel worker for k-d tree (Euclidean) kNCN
struct KncnKdtreeWorker : public Worker {
  const RMatrix<int> species_pa;
  const PointCloud2D& cloud;
  const RVector<double> x;
  const RVector<double> y;
  const RVector<int> seeds;
  RMatrix<int> curves;

  KncnKdtreeWorker(const IntegerMatrix& sp,
                   const PointCloud2D& cl,
                   const NumericVector& x_,
                   const NumericVector& y_,
                   const IntegerVector& s,
                   IntegerMatrix& c)
    : species_pa(sp), cloud(cl), x(x_), y(y_), seeds(s), curves(c) {}

  void operator()(std::size_t begin, std::size_t end) {
    int n_sites = species_pa.nrow();
    int n_species = species_pa.ncol();

    KDTree2D tree(2, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    tree.buildIndex();

    for (std::size_t s = begin; s < end; s++) {
      std::vector<bool> visited(n_sites, false);
      std::vector<bool> species_seen(n_species, false);
      int richness = 0;

      double sum_x = x[seeds[s]];
      double sum_y = y[seeds[s]];
      int n_visited = 1;

      int current = seeds[s];
      visited[current] = true;

      for (int sp = 0; sp < n_species; sp++) {
        if (species_pa(current, sp) > 0) {
          species_seen[sp] = true;
          richness++;
        }
      }
      curves(s, 0) = richness;

      for (int step = 1; step < n_sites; step++) {
        double cx = sum_x / n_visited;
        double cy = sum_y / n_visited;

        int nearest = find_nearest_unvisited(tree, cx, cy, visited, n_sites);

        current = nearest;
        visited[current] = true;
        sum_x += x[current];
        sum_y += y[current];
        n_visited++;

        for (int sp = 0; sp < n_species; sp++) {
          if (species_pa(current, sp) > 0 && !species_seen[sp]) {
            species_seen[sp] = true;
            richness++;
          }
        }
        curves(s, step) = richness;
      }
    }
  }
};


// Parallel worker for ball tree (haversine) kNCN
struct KncnBalltreeWorker : public Worker {
  const RMatrix<int> species_pa;
  const std::vector<double>& px;
  const std::vector<double>& py;
  const RVector<int> seeds;
  RMatrix<int> curves;

  KncnBalltreeWorker(const IntegerMatrix& sp,
                     const std::vector<double>& x,
                     const std::vector<double>& y,
                     const IntegerVector& s,
                     IntegerMatrix& c)
    : species_pa(sp), px(x), py(y), seeds(s), curves(c) {}

  void operator()(std::size_t begin, std::size_t end) {
    int n_sites = species_pa.nrow();
    int n_species = species_pa.ncol();

    BallTree btree(px, py);

    for (std::size_t s = begin; s < end; s++) {
      std::vector<bool> visited(n_sites, false);
      std::vector<bool> species_seen(n_species, false);
      int richness = 0;

      double sum_x = px[seeds[s]];
      double sum_y = py[seeds[s]];
      int n_visited = 1;

      int current = seeds[s];
      visited[current] = true;

      for (int sp = 0; sp < n_species; sp++) {
        if (species_pa(current, sp) > 0) {
          species_seen[sp] = true;
          richness++;
        }
      }
      curves(s, 0) = richness;

      for (int step = 1; step < n_sites; step++) {
        double cx = sum_x / n_visited;
        double cy = sum_y / n_visited;

        int nearest = btree.find_nearest_unvisited(cx, cy, visited);

        current = nearest;
        visited[current] = true;
        sum_x += px[current];
        sum_y += py[current];
        n_visited++;

        for (int sp = 0; sp < n_species; sp++) {
          if (species_pa(current, sp) > 0 && !species_seen[sp]) {
            species_seen[sp] = true;
            richness++;
          }
        }
        curves(s, step) = richness;
      }
    }
  }
};


// [[Rcpp::export]]
IntegerMatrix cpp_kncn_kdtree_parallel(IntegerMatrix species_pa,
                                        NumericVector x, NumericVector y,
                                        int n_seeds,
                                        int n_cores = 1,
                                        bool progress = false,
                                        std::string distance = "euclidean") {
  int n_sites = species_pa.nrow();
  IntegerMatrix curves(n_seeds, n_sites);
  IntegerVector seeds = Rcpp::sample(n_sites, n_seeds, true) - 1;
  bool use_haversine = (distance == "haversine");

  std::vector<double> vx(n_sites), vy(n_sites);
  for (int i = 0; i < n_sites; i++) { vx[i] = x[i]; vy[i] = y[i]; }

  if (use_haversine) {
    if (n_cores > 1) {
      KncnBalltreeWorker worker(species_pa, vx, vy, seeds, curves);
      parallelFor(0, n_seeds, worker);
    } else {
      for (int s = 0; s < n_seeds; s++) {
        IntegerVector curve = cpp_kncn_kdtree_single(species_pa, x, y, seeds[s], distance);
        curves(s, _) = curve;
      }
    }
  } else {
    PointCloud2D cloud;
    cloud.pts_x = vx;
    cloud.pts_y = vy;

    if (n_cores > 1) {
      KncnKdtreeWorker worker(species_pa, cloud, x, y, seeds, curves);
      parallelFor(0, n_seeds, worker);
    } else {
      for (int s = 0; s < n_seeds; s++) {
        IntegerVector curve = cpp_kncn_kdtree_single(species_pa, x, y, seeds[s], distance);
        curves(s, _) = curve;
      }
    }
  }

  return curves;
}
