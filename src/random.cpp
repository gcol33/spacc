// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <set>
#include <algorithm>
#include "parallel_utils.h"
using namespace Rcpp;
using namespace RcppParallel;


// [[Rcpp::export]]
IntegerVector cpp_random_single(IntegerMatrix species_pa, IntegerVector order) {
  int n_sites = species_pa.nrow();
  int n_species = species_pa.ncol();

  IntegerVector curve(n_sites);
  std::set<int> species_seen;

  for (int step = 0; step < n_sites; step++) {
    int site = order[step];

    // Add species from this site
    for (int sp = 0; sp < n_species; sp++) {
      if (species_pa(site, sp) > 0) {
        species_seen.insert(sp);
      }
    }
    curve[step] = species_seen.size();
  }

  return curve;
}


// Worker struct for parallel random accumulation
struct RandomWorker : public Worker {
  // Inputs (read-only)
  const RMatrix<int> species_pa;
  const RMatrix<int> orders;  // n_seeds x n_sites

  // Output
  RMatrix<int> curves;

  // Constructor
  RandomWorker(const IntegerMatrix& species_pa_,
               const IntegerMatrix& orders_,
               IntegerMatrix& curves_)
    : species_pa(species_pa_), orders(orders_), curves(curves_) {}

  // Parallel operator
  void operator()(std::size_t begin, std::size_t end) {
    int n_sites = species_pa.nrow();
    int n_species = species_pa.ncol();

    for (std::size_t s = begin; s < end; s++) {
      std::set<int> species_seen;

      for (int step = 0; step < n_sites; step++) {
        int site = orders(s, step);

        // Add species from this site
        for (int sp = 0; sp < n_species; sp++) {
          if (species_pa(site, sp) > 0) {
            species_seen.insert(sp);
          }
        }
        curves(s, step) = species_seen.size();
      }
    }
  }
};


// [[Rcpp::export]]
IntegerMatrix cpp_random_parallel(IntegerMatrix species_pa,
                                  int n_seeds,
                                  int n_cores = 1,
                                  bool progress = false) {
  int n_sites = species_pa.nrow();
  IntegerMatrix curves(n_seeds, n_sites);

  // Generate random orders for each seed
  IntegerMatrix orders(n_seeds, n_sites);
  for (int s = 0; s < n_seeds; s++) {
    IntegerVector order = Rcpp::sample(n_sites, n_sites, false) - 1;
    orders(s, _) = order;
  }

  RandomWorker worker(species_pa, orders, curves);
  dispatch_parallel(n_seeds, n_cores, worker, curves,
    [&](int s) { IntegerVector order = orders(s, _); return cpp_random_single(species_pa, order); });
  return curves;
}


// [[Rcpp::export]]
IntegerMatrix cpp_order_parallel(IntegerMatrix species_pa,
                                 IntegerMatrix orders,
                                 int n_cores = 1,
                                 bool progress = false) {
  int n_sites = species_pa.nrow();
  int n_orders = orders.nrow();
  IntegerMatrix curves(n_orders, n_sites);

  // Reuse the random-accumulation worker with caller-supplied orderings
  // (0-indexed, one ordering per row) instead of internally generated ones.
  RandomWorker worker(species_pa, orders, curves);
  dispatch_parallel(n_orders, n_cores, worker, curves,
    [&](int s) { IntegerVector order = orders(s, _); return cpp_random_single(species_pa, order); });
  return curves;
}
