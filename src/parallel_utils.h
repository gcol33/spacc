// parallel_utils.h
// Shared parallel dispatch logic for spacc accumulation methods
// Eliminates boilerplate in *_parallel functions

#ifndef SPACC_PARALLEL_UTILS_H
#define SPACC_PARALLEL_UTILS_H

#include <Rcpp.h>
#include <RcppParallel.h>

// Sample random seed sites (0-indexed)
inline Rcpp::IntegerVector sample_seeds(int n_sites, int n_seeds) {
  return Rcpp::sample(n_sites, n_seeds, true) - 1;
}

// Run seeds in parallel (via Worker) or sequentially (via single function).
// Worker must already be constructed with a reference to the output matrix.
// SingleFn: callable(int s) -> IntegerVector, where s is the seed index.
template<typename Worker, typename SingleFn>
void dispatch_parallel(int n_seeds, int n_cores,
                       Worker& worker,
                       Rcpp::IntegerMatrix& curves,
                       SingleFn run_single) {
  if (n_cores > 1) {
    RcppParallel::parallelFor(0, n_seeds, worker);
  } else {
    for (int s = 0; s < n_seeds; s++) {
      curves(s, Rcpp::_) = run_single(s);
    }
  }
}

#endif // SPACC_PARALLEL_UTILS_H
