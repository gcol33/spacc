// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <utility>
#include "parallel_utils.h"
#include "core/geometry_core.h"
using namespace Rcpp;
using namespace RcppParallel;

// Species richness accumulated along site orders (one curve per order).
struct RichnessWorker : public Worker {
  const RMatrix<int> species;
  const RMatrix<int> orders;
  RMatrix<int> curves;

  RichnessWorker(const IntegerMatrix& species_, const IntegerMatrix& orders_,
                 IntegerMatrix& curves_)
    : species(species_), orders(orders_), curves(curves_) {}

  void operator()(std::size_t begin, std::size_t end) {
    const int n_sites = species.nrow();
    const int n_species = species.ncol();
    auto species_at = [&](int site, int sp) { return species(site, sp); };
    for (std::size_t s = begin; s < end; s++) {
      spacc::richness_along_order(
        species_at, n_sites, n_species,
        [&](int step) { return orders(s, step); },
        [&](int step, int richness) { curves(s, step) = richness; });
    }
  }
};

// [[Rcpp::export]]
IntegerMatrix cpp_order_parallel(IntegerMatrix species_pa, IntegerMatrix orders,
                                 int n_cores = 1) {
  IntegerMatrix curves(orders.nrow(), species_pa.nrow());
  RichnessWorker worker(species_pa, orders, curves);
  run_workers(static_cast<std::size_t>(orders.nrow()), n_cores, worker);
  return curves;
}


// Species richness within expanding radii around each seed (wavefront and
// distance-decay curves).
struct RadiusWorker : public Worker {
  const RMatrix<int> species;
  const RMatrix<double> dist;
  const RVector<int> seeds;
  const RVector<double> radii;
  RMatrix<int> curves;
  RMatrix<int> sites_included;

  RadiusWorker(const IntegerMatrix& species_, const NumericMatrix& dist_,
               const IntegerVector& seeds_, const NumericVector& radii_,
               IntegerMatrix& curves_, IntegerMatrix& sites_included_)
    : species(species_), dist(dist_), seeds(seeds_), radii(radii_),
      curves(curves_), sites_included(sites_included_) {}

  void operator()(std::size_t begin, std::size_t end) {
    const int n_sites = species.nrow();
    const int n_species = species.ncol();
    auto species_at = [&](int site, int sp) { return species(site, sp); };
    auto dist_at = [&](int i, int j) { return dist(i, j); };
    for (std::size_t s = begin; s < end; s++) {
      for (std::size_t r = 0; r < radii.size(); r++) {
        std::pair<int, int> counts = spacc::richness_within_radius(
          species_at, dist_at, n_sites, n_species, seeds[s], radii[r]);
        curves(s, r) = counts.first;
        sites_included(s, r) = counts.second;
      }
    }
  }
};

// [[Rcpp::export]]
List cpp_radius_parallel(IntegerMatrix species_pa, NumericMatrix dist_mat,
                         IntegerVector seeds, NumericVector radii, int n_cores = 1) {
  IntegerMatrix curves(seeds.size(), radii.size());
  IntegerMatrix sites_included(seeds.size(), radii.size());
  RadiusWorker worker(species_pa, dist_mat, seeds, radii, curves, sites_included);
  run_workers(static_cast<std::size_t>(seeds.size()), n_cores, worker);
  return List::create(Named("curves") = curves,
                      Named("sites_included") = sites_included);
}


// Cumulative area of the sites visited along each order (n_orders x n_steps,
// orders 0-based): the convex hull of the first k sites, or the summed
// Voronoi cells of the first k sites within the bounding box of all sites.
// [[Rcpp::export]]
NumericMatrix cpp_area_along_orders(NumericVector x, NumericVector y, IntegerMatrix orders,
                                    std::string method) {
  const int n_orders = orders.nrow(), n_steps = orders.ncol();
  NumericMatrix area(n_orders, n_steps);
  auto xs = [&](int i) { return x[i]; };
  auto ys = [&](int i) { return y[i]; };
  if (method == "convex_hull") {
    std::vector<double> row(n_steps);
    for (int s = 0; s < n_orders; s++) {
      spacc::hull_area_along_order(xs, ys, [&](int k) { return orders(s, k); }, n_steps,
                                   row.data());
      for (int k = 0; k < n_steps; k++) area(s, k) = row[k];
    }
  } else if (method == "voronoi") {
    const std::vector<double> cell = spacc::voronoi_cell_areas(
      xs, ys, x.size(), Rcpp::min(x), Rcpp::max(x), Rcpp::min(y), Rcpp::max(y));
    for (int s = 0; s < n_orders; s++) {
      double total = 0.0;
      for (int k = 0; k < n_steps; k++) {
        total += cell[orders(s, k)];
        area(s, k) = total;
      }
    }
  } else {
    Rcpp::stop("unknown area method");
  }
  return area;
}
