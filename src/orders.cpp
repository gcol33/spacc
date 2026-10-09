// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <string>
#include <vector>
#include "parallel_utils.h"
#include "kdtree_adapter.h"
#include "balltree.h"
using namespace Rcpp;
using namespace RcppParallel;

// Site-visitation orders, one row per seed (or focal point), 0-based.
// Every accumulation method in spacc is an order generator; the species data
// are only read afterwards, by cpp_order_parallel() and the metric workers.
//
// A generator provides
//   Context context() const;                       built once per worker chunk
//   void generate(Context&, std::size_t s, std::vector<int>& order) const;
// and the OrderWorker below writes each order into its row.

template<typename Gen>
struct OrderWorker : public Worker {
  const Gen gen;
  RMatrix<int> orders;

  OrderWorker(const Gen& gen_, IntegerMatrix& orders_) : gen(gen_), orders(orders_) {}

  void operator()(std::size_t begin, std::size_t end) {
    auto ctx = gen.context();
    std::vector<int> order;
    for (std::size_t s = begin; s < end; s++) {
      gen.generate(ctx, s, order);
      for (std::size_t k = 0; k < order.size(); k++) orders(s, k) = order[k];
    }
  }
};

template<typename Gen>
IntegerMatrix generate_orders(const Gen& gen, int n_orders, int n_sites, int n_cores) {
  IntegerMatrix orders(n_orders, n_sites);
  OrderWorker<Gen> worker(gen, orders);
  run_workers(static_cast<std::size_t>(n_orders), n_cores, worker);
  return orders;
}

struct NoContext {};


// ---------------------------------------------------------------------------
// Fixed focus: rank all sites by distance from one continuous focal point.
// ---------------------------------------------------------------------------

struct FocalGen {
  std::vector<double> x, y;
  RVector<double> fx, fy;
  spacc::PointDistance dist;

  typedef NoContext Context;
  Context context() const { return Context(); }

  void generate(Context&, std::size_t f, std::vector<int>& order) const {
    const double qx = fx[f], qy = fy[f];
    spacc::focal_order(static_cast<int>(x.size()),
      [&](int site) { return dist(qx, qy, x[site], y[site]); }, order);
  }
};

// [[Rcpp::export]]
IntegerMatrix cpp_focal_order(NumericVector x, NumericVector y,
                              NumericVector focus_x, NumericVector focus_y,
                              std::string distance = "euclidean",
                              int n_cores = 1) {
  FocalGen gen{to_std(x), to_std(y), RVector<double>(focus_x), RVector<double>(focus_y),
               spacc::PointDistance{distance == "haversine"}};
  return generate_orders(gen, focus_x.size(), x.size(), n_cores);
}


// ---------------------------------------------------------------------------
// Nearest-neighbour walk
// ---------------------------------------------------------------------------

struct NnWalkMatrixGen {
  RMatrix<double> d;
  RVector<int> seeds;

  typedef NoContext Context;
  Context context() const { return Context(); }

  void generate(Context&, std::size_t s, std::vector<int>& order) const {
    const int n = static_cast<int>(d.nrow());
    auto dist = [&](int i, int j) { return d(i, j); };
    spacc::nn_walk_order(n, seeds[s],
      [&](int current, const spacc::VisitedFlags& visited) {
        return spacc::nearest_unvisited(dist, current, n, visited);
      }, order);
  }
};

// [[Rcpp::export]]
IntegerMatrix cpp_nn_walk_order(NumericMatrix dist_mat, IntegerVector seeds,
                                int n_cores = 1) {
  NnWalkMatrixGen gen{RMatrix<double>(dist_mat), RVector<int>(seeds)};
  return generate_orders(gen, seeds.size(), dist_mat.nrow(), n_cores);
}

template<typename Tree>
struct NnWalkTreeGen {
  std::vector<double> x, y;
  RVector<int> seeds;

  typedef Tree Context;
  Context context() const { return Tree(x, y); }

  void generate(Context& tree, std::size_t s, std::vector<int>& order) const {
    spacc::nn_walk_order(static_cast<int>(x.size()), seeds[s],
      [&](int current, const spacc::VisitedFlags& visited) {
        return tree.nearest_unvisited(x[current], y[current], visited);
      }, order);
  }
};

// [[Rcpp::export]]
IntegerMatrix cpp_nn_walk_tree_order(NumericVector x, NumericVector y,
                                     IntegerVector seeds,
                                     std::string distance = "euclidean",
                                     int n_cores = 1) {
  if (distance == "haversine") {
    NnWalkTreeGen<BallTree> gen{to_std(x), to_std(y), RVector<int>(seeds)};
    return generate_orders(gen, seeds.size(), x.size(), n_cores);
  }
  NnWalkTreeGen<KdTree> gen{to_std(x), to_std(y), RVector<int>(seeds)};
  return generate_orders(gen, seeds.size(), x.size(), n_cores);
}


// ---------------------------------------------------------------------------
// Nearest-centroid walk (kNCN)
// ---------------------------------------------------------------------------

struct KncnGen {
  std::vector<double> x, y;
  RVector<int> seeds;
  spacc::PointDistance dist;

  typedef NoContext Context;
  Context context() const { return Context(); }

  void generate(Context&, std::size_t s, std::vector<int>& order) const {
    spacc::kncn_order(x, y, seeds[s],
      [&](double cx, double cy, const spacc::VisitedFlags& visited) {
        return spacc::nearest_unvisited_point(cx, cy, x, y, visited, dist);
      }, order);
  }
};

// [[Rcpp::export]]
IntegerMatrix cpp_kncn_order(NumericVector x, NumericVector y, IntegerVector seeds,
                             std::string distance = "euclidean", int n_cores = 1) {
  KncnGen gen{to_std(x), to_std(y), RVector<int>(seeds),
              spacc::PointDistance{distance == "haversine"}};
  return generate_orders(gen, seeds.size(), x.size(), n_cores);
}

template<typename Tree>
struct KncnTreeGen {
  std::vector<double> x, y;
  RVector<int> seeds;

  typedef Tree Context;
  Context context() const { return Tree(x, y); }

  void generate(Context& tree, std::size_t s, std::vector<int>& order) const {
    spacc::kncn_order(x, y, seeds[s],
      [&](double cx, double cy, const spacc::VisitedFlags& visited) {
        return tree.nearest_unvisited(cx, cy, visited);
      }, order);
  }
};

// [[Rcpp::export]]
IntegerMatrix cpp_kncn_tree_order(NumericVector x, NumericVector y, IntegerVector seeds,
                                  std::string distance = "euclidean", int n_cores = 1) {
  if (distance == "haversine") {
    KncnTreeGen<BallTree> gen{to_std(x), to_std(y), RVector<int>(seeds)};
    return generate_orders(gen, seeds.size(), x.size(), n_cores);
  }
  KncnTreeGen<KdTree> gen{to_std(x), to_std(y), RVector<int>(seeds)};
  return generate_orders(gen, seeds.size(), x.size(), n_cores);
}


// ---------------------------------------------------------------------------
// Gaussian-weighted walk. Draws come from R's RNG, so set.seed() reproduces
// the orders; this runs serially because R's RNG is not thread-safe.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
IntegerMatrix cpp_gaussian_order(NumericMatrix dist_mat, IntegerVector seeds, double sigma) {
  const int n = dist_mat.nrow();
  IntegerMatrix orders(seeds.size(), n);
  std::vector<int> order;
  auto dist = [&](int i, int j) { return dist_mat(i, j); };
  auto unif = []() { return R::unif_rand(); };
  for (int s = 0; s < seeds.size(); s++) {
    spacc::gaussian_order(dist, n, seeds[s], sigma, unif, order);
    for (int k = 0; k < n; k++) orders(s, k) = order[k];
  }
  return orders;
}


// ---------------------------------------------------------------------------
// Directional cone from each seed along a given angle
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
IntegerMatrix cpp_cone_order(NumericVector x, NumericVector y, IntegerVector seeds,
                             NumericVector angles, double width) {
  const std::vector<double> vx = to_std(x), vy = to_std(y);
  const int n = x.size();
  IntegerMatrix orders(seeds.size(), n);
  std::vector<int> order;
  for (int s = 0; s < seeds.size(); s++) {
    spacc::cone_order(vx, vy, seeds[s], angles[s], width, order);
    for (int k = 0; k < n; k++) orders(s, k) = order[k];
  }
  return orders;
}
