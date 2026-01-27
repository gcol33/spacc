// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <set>
#include "kdtree_adapter.h"
#include "balltree.h"
using namespace Rcpp;
using namespace RcppParallel;


// ============================================================================
// EXACT (brute-force) kNN — uses precomputed distance matrix
// ============================================================================

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


// Worker struct for parallel exact kNN
struct KnnWorker : public Worker {
  const RMatrix<int> species_pa;
  const RMatrix<double> dist_mat;
  const RVector<int> seeds;
  RMatrix<int> curves;

  KnnWorker(const IntegerMatrix& species_pa_,
            const NumericMatrix& dist_mat_,
            const IntegerVector& seeds_,
            IntegerMatrix& curves_)
    : species_pa(species_pa_), dist_mat(dist_mat_),
      seeds(seeds_), curves(curves_) {}

  void operator()(std::size_t begin, std::size_t end) {
    int n_sites = species_pa.nrow();
    int n_species = species_pa.ncol();

    for (std::size_t s = begin; s < end; s++) {
      std::vector<bool> visited(n_sites, false);
      std::set<int> species_seen;

      int current = seeds[s];
      visited[current] = true;

      for (int sp = 0; sp < n_species; sp++) {
        if (species_pa(current, sp) > 0) {
          species_seen.insert(sp);
        }
      }
      curves(s, 0) = species_seen.size();

      for (int step = 1; step < n_sites; step++) {
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


// [[Rcpp::export]]
IntegerMatrix cpp_knn_parallel(IntegerMatrix species_pa,
                               NumericMatrix dist_mat,
                               int n_seeds,
                               int n_cores = 1,
                               bool progress = false) {
  int n_sites = species_pa.nrow();
  IntegerMatrix curves(n_seeds, n_sites);

  IntegerVector seeds = Rcpp::sample(n_sites, n_seeds, true) - 1;

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


// ============================================================================
// SPATIAL TREE kNN — k-d tree (Euclidean) or ball tree (haversine)
// No precomputed distance matrix needed.
// ============================================================================

// Helper: build coordinate vectors from R vectors
inline void fill_coords(std::vector<double>& vx, std::vector<double>& vy,
                         NumericVector x, NumericVector y, int n) {
  vx.resize(n);
  vy.resize(n);
  for (int i = 0; i < n; i++) {
    vx[i] = x[i];
    vy[i] = y[i];
  }
}

// Run a single kNN accumulation curve using a spatial tree.
// qx/qy for the query are taken from the coordinate arrays.
// find_nn is a lambda: (double qx, double qy, visited) -> int nearest
template<typename FindNN>
IntegerVector knn_tree_single_impl(const IntegerMatrix& species_pa,
                                    const std::vector<double>& px,
                                    const std::vector<double>& py,
                                    int seed, FindNN find_nn) {
  int n_sites = species_pa.nrow();
  int n_species = species_pa.ncol();

  IntegerVector curve(n_sites);
  std::vector<bool> visited(n_sites, false);
  std::set<int> species_seen;

  int current = seed;
  visited[current] = true;

  for (int sp = 0; sp < n_species; sp++) {
    if (species_pa(current, sp) > 0) species_seen.insert(sp);
  }
  curve[0] = species_seen.size();

  for (int step = 1; step < n_sites; step++) {
    int nearest = find_nn(px[current], py[current], visited);
    current = nearest;
    visited[current] = true;

    for (int sp = 0; sp < n_species; sp++) {
      if (species_pa(current, sp) > 0) species_seen.insert(sp);
    }
    curve[step] = species_seen.size();
  }

  return curve;
}


// [[Rcpp::export]]
IntegerVector cpp_knn_kdtree_single(IntegerMatrix species_pa,
                                     NumericVector x, NumericVector y,
                                     int seed,
                                     std::string distance = "euclidean") {
  int n_sites = species_pa.nrow();
  bool use_haversine = (distance == "haversine");

  std::vector<double> vx, vy;
  fill_coords(vx, vy, x, y, n_sites);

  if (use_haversine) {
    BallTree btree(vx, vy);
    return knn_tree_single_impl(species_pa, vx, vy, seed,
      [&](double qx, double qy, const std::vector<bool>& vis) {
        return btree.find_nearest_unvisited(qx, qy, vis);
      });
  } else {
    PointCloud2D cloud;
    cloud.pts_x = vx;
    cloud.pts_y = vy;
    KDTree2D* tree = build_kdtree(cloud);
    auto result = knn_tree_single_impl(species_pa, vx, vy, seed,
      [&](double qx, double qy, const std::vector<bool>& vis) {
        return find_nearest_unvisited(*tree, qx, qy, vis, n_sites);
      });
    delete tree;
    return result;
  }
}


// Parallel worker for k-d tree (Euclidean) kNN
struct KnnKdtreeWorker : public Worker {
  const RMatrix<int> species_pa;
  const PointCloud2D& cloud;
  const RVector<int> seeds;
  RMatrix<int> curves;

  KnnKdtreeWorker(const IntegerMatrix& sp,
                  const PointCloud2D& cl,
                  const IntegerVector& s,
                  IntegerMatrix& c)
    : species_pa(sp), cloud(cl), seeds(s), curves(c) {}

  void operator()(std::size_t begin, std::size_t end) {
    int n_sites = species_pa.nrow();
    int n_species = species_pa.ncol();

    KDTree2D tree(2, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    tree.buildIndex();

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
        int nearest = find_nearest_unvisited(
          tree, cloud.pts_x[current], cloud.pts_y[current],
          visited, n_sites
        );
        current = nearest;
        visited[current] = true;

        for (int sp = 0; sp < n_species; sp++) {
          if (species_pa(current, sp) > 0) species_seen.insert(sp);
        }
        curves(s, step) = species_seen.size();
      }
    }
  }
};


// Parallel worker for ball tree (haversine) kNN
struct KnnBalltreeWorker : public Worker {
  const RMatrix<int> species_pa;
  const std::vector<double>& px;
  const std::vector<double>& py;
  const RVector<int> seeds;
  RMatrix<int> curves;

  KnnBalltreeWorker(const IntegerMatrix& sp,
                    const std::vector<double>& x,
                    const std::vector<double>& y,
                    const IntegerVector& s,
                    IntegerMatrix& c)
    : species_pa(sp), px(x), py(y), seeds(s), curves(c) {}

  void operator()(std::size_t begin, std::size_t end) {
    int n_sites = species_pa.nrow();
    int n_species = species_pa.ncol();

    // Each thread builds its own ball tree
    BallTree btree(px, py);

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
        int nearest = btree.find_nearest_unvisited(px[current], py[current], visited);
        current = nearest;
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
IntegerMatrix cpp_knn_kdtree_parallel(IntegerMatrix species_pa,
                                       NumericVector x, NumericVector y,
                                       int n_seeds,
                                       int n_cores = 1,
                                       bool progress = false,
                                       std::string distance = "euclidean") {
  int n_sites = species_pa.nrow();
  IntegerMatrix curves(n_seeds, n_sites);
  IntegerVector seeds = Rcpp::sample(n_sites, n_seeds, true) - 1;
  bool use_haversine = (distance == "haversine");

  std::vector<double> vx, vy;
  fill_coords(vx, vy, x, y, n_sites);

  if (use_haversine) {
    if (n_cores > 1) {
      KnnBalltreeWorker worker(species_pa, vx, vy, seeds, curves);
      parallelFor(0, n_seeds, worker);
    } else {
      for (int s = 0; s < n_seeds; s++) {
        IntegerVector curve = cpp_knn_kdtree_single(species_pa, x, y, seeds[s], distance);
        curves(s, _) = curve;
      }
    }
  } else {
    PointCloud2D cloud;
    cloud.pts_x = vx;
    cloud.pts_y = vy;

    if (n_cores > 1) {
      KnnKdtreeWorker worker(species_pa, cloud, seeds, curves);
      parallelFor(0, n_seeds, worker);
    } else {
      for (int s = 0; s < n_seeds; s++) {
        IntegerVector curve = cpp_knn_kdtree_single(species_pa, x, y, seeds[s], distance);
        curves(s, _) = curve;
      }
    }
  }

  return curves;
}


// [[Rcpp::export]]
IntegerMatrix cpp_knn_kdtree_parallel_seeds(IntegerMatrix species_pa,
                                             NumericVector x, NumericVector y,
                                             IntegerVector seeds,
                                             int n_cores = 1,
                                             bool progress = false,
                                             std::string distance = "euclidean") {
  int n_sites = species_pa.nrow();
  int n_seeds = seeds.size();
  IntegerMatrix curves(n_seeds, n_sites);
  bool use_haversine = (distance == "haversine");

  std::vector<double> vx, vy;
  fill_coords(vx, vy, x, y, n_sites);

  if (use_haversine) {
    if (n_cores > 1) {
      KnnBalltreeWorker worker(species_pa, vx, vy, seeds, curves);
      parallelFor(0, n_seeds, worker);
    } else {
      for (int s = 0; s < n_seeds; s++) {
        IntegerVector curve = cpp_knn_kdtree_single(species_pa, x, y, seeds[s], distance);
        curves(s, _) = curve;
      }
    }
  } else {
    PointCloud2D cloud;
    cloud.pts_x = vx;
    cloud.pts_y = vy;

    if (n_cores > 1) {
      KnnKdtreeWorker worker(species_pa, cloud, seeds, curves);
      parallelFor(0, n_seeds, worker);
    } else {
      for (int s = 0; s < n_seeds; s++) {
        IntegerVector curve = cpp_knn_kdtree_single(species_pa, x, y, seeds[s], distance);
        curves(s, _) = curve;
      }
    }
  }

  return curves;
}
