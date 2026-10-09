// parallel_utils.h
// Shared RcppParallel driver for the accumulation workers.
//
// Every per-seed computation in spacc is a Worker whose operator()(begin, end)
// handles a range of seeds (or orders). run_workers() runs it on the
// RcppParallel pool when more than one core is requested and serially
// otherwise, so a single implementation serves both paths.

#ifndef SPACC_PARALLEL_UTILS_H
#define SPACC_PARALLEL_UTILS_H

#include <Rcpp.h>
#include <RcppParallel.h>
#include <cstdio>
#include <vector>
#include "core/accumulation_core.h"

template<typename Worker>
inline void run_workers(std::size_t n, int n_cores, Worker& worker) {
  if (n_cores > 1 && n > 1) {
    RcppParallel::parallelFor(0, n, worker);
  } else {
    worker(0, n);
  }
}

// Split a flat (n_seeds * k) x n_sites matrix, row s * k + i holding seed s
// and component i, into a named list of k matrices (n_seeds x n_sites).
inline Rcpp::List split_components(const Rcpp::NumericMatrix& flat,
                                   int n_seeds, int k,
                                   const Rcpp::CharacterVector& names) {
  const int n_sites = flat.ncol();
  Rcpp::List out(k);
  for (int i = 0; i < k; i++) {
    Rcpp::NumericMatrix mat(n_seeds, n_sites);
    for (int s = 0; s < n_seeds; s++) {
      for (int st = 0; st < n_sites; st++) mat(s, st) = flat(s * k + i, st);
    }
    out[i] = mat;
  }
  out.attr("names") = names;
  return out;
}

// Names "q0.0", "q1.0", ... for a vector of Hill orders.
inline Rcpp::CharacterVector q_names(const Rcpp::NumericVector& q_values) {
  Rcpp::CharacterVector names(q_values.size());
  for (int qi = 0; qi < q_values.size(); qi++) {
    char name[16];
    snprintf(name, sizeof(name), "q%.1f", q_values[qi]);
    names[qi] = name;
  }
  return names;
}

// Copy R coordinate vectors into std::vector (the spatial indices and the
// pure core routines read plain vectors).
inline std::vector<double> to_std(const Rcpp::NumericVector& v) {
  return std::vector<double>(v.begin(), v.end());
}

// ---------------------------------------------------------------------------
// Generic "metric along an order" worker.
//
// Metric must provide
//   struct Scratch;  Scratch scratch(int n_species) const;   per-order state
//   void at_step(std::size_t s, int step, int site,
//                const std::vector<T>& cumulative, Scratch& scratch);
// and own (as RMatrix/RVector views) the output matrices it writes into; the
// only mutation in at_step is through those views, each order writing its own
// rows, so the shared metric object is safe across threads.
// The worker walks each order, accumulating site rows into `cumulative`
// (type T) and calling at_step after each site.
// ---------------------------------------------------------------------------

template<typename T, typename SpeciesView, typename Metric>
struct OrderMetricWorker : public RcppParallel::Worker {
  const SpeciesView species;
  const RcppParallel::RMatrix<int> orders;
  Metric metric;

  template<typename SpeciesMat>
  OrderMetricWorker(const SpeciesMat& species_, const Rcpp::IntegerMatrix& orders_,
                    const Metric& metric_)
    : species(species_), orders(orders_), metric(metric_) {}

  void operator()(std::size_t begin, std::size_t end) {
    const int n_sites = static_cast<int>(species.nrow());
    const int n_species = static_cast<int>(species.ncol());
    auto species_at = [&](int site, int sp) { return species(site, sp); };
    for (std::size_t s = begin; s < end; s++) {
      typename Metric::Scratch scratch = metric.scratch(n_species);
      std::vector<T> cumulative(n_species, T(0));
      auto order_at = [&](int step) { return orders(s, step); };
      spacc::accumulate_along_order(
        species_at, n_sites, n_species, order_at, cumulative,
        [&](int step, int site) { metric.at_step(s, step, site, cumulative, scratch); });
    }
  }
};

template<typename T, typename SpeciesMat, typename Metric>
void run_order_metric(const SpeciesMat& species, const Rcpp::IntegerMatrix& orders,
                      const Metric& metric, int n_cores) {
  typedef RcppParallel::RMatrix<typename SpeciesMat::stored_type> View;
  OrderMetricWorker<T, View, Metric> worker(species, orders, metric);
  run_workers(static_cast<std::size_t>(orders.nrow()), n_cores, worker);
}

#endif // SPACC_PARALLEL_UTILS_H
