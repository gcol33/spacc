// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include "core/beta_core.h"
#include "parallel_utils.h"
using namespace Rcpp;
using namespace RcppParallel;


// Baselga (2010) partition of the dissimilarity between the accumulated
// species pool and each newly added site, along each order. Step k compares
// the pool of the first k sites with site k + 1, so outputs have n_sites - 1
// columns.
struct BetaMetric {
  RMatrix<int> species;
  RMatrix<double> dist;
  RMatrix<double> beta_total;
  RMatrix<double> beta_turnover;
  RMatrix<double> beta_nestedness;
  RMatrix<double> distance;
  bool use_jaccard;

  BetaMetric(const IntegerMatrix& species_, const NumericMatrix& dist_,
             NumericMatrix& total_, NumericMatrix& turnover_,
             NumericMatrix& nestedness_, NumericMatrix& distance_, bool jaccard)
    : species(species_), dist(dist_), beta_total(total_), beta_turnover(turnover_),
      beta_nestedness(nestedness_), distance(distance_), use_jaccard(jaccard) {}

  struct Scratch {
    std::vector<bool> pool;   // species present in the accumulated sites
    int current;
    double cum_distance;
  };
  Scratch scratch(int n_species) const {
    Scratch sc;
    sc.pool.assign(n_species, false);
    sc.current = -1;
    sc.cum_distance = 0.0;
    return sc;
  }

  void at_step(std::size_t s, int step, int site, const std::vector<int>&,
               Scratch& sc) {
    const int n_species = species.ncol();
    if (step == 0) {
      for (int sp = 0; sp < n_species; sp++) sc.pool[sp] = species(site, sp) > 0;
      sc.current = site;
      return;
    }

    sc.cum_distance += dist(sc.current, site);
    distance(s, step - 1) = sc.cum_distance;

    int a = 0, b = 0, c = 0;
    for (int sp = 0; sp < n_species; sp++) {
      const bool in_pool = sc.pool[sp];
      const bool in_site = species(site, sp) > 0;
      if (in_pool && in_site) a++;
      else if (in_pool) b++;
      else if (in_site) c++;
      if (in_site) sc.pool[sp] = true;
    }

    spacc::BetaComponents beta = use_jaccard
      ? spacc::calc_beta_jaccard(a, b, c)
      : spacc::calc_beta_sorensen(a, b, c);
    beta_total(s, step - 1) = beta.total;
    beta_turnover(s, step - 1) = beta.turnover;
    beta_nestedness(s, step - 1) = beta.nestedness;
    sc.current = site;
  }
};

// [[Rcpp::export]]
List cpp_beta_order_parallel(IntegerMatrix species_pa, NumericMatrix dist_mat,
                             IntegerMatrix orders, bool use_jaccard = false,
                             int n_cores = 1) {
  const int n_seeds = orders.nrow();
  const int n_steps = species_pa.nrow() - 1;
  NumericMatrix beta_total(n_seeds, n_steps);
  NumericMatrix beta_turnover(n_seeds, n_steps);
  NumericMatrix beta_nestedness(n_seeds, n_steps);
  NumericMatrix distance(n_seeds, n_steps);
  BetaMetric metric(species_pa, dist_mat, beta_total, beta_turnover,
                    beta_nestedness, distance, use_jaccard);
  run_order_metric<int>(species_pa, orders, metric, n_cores);
  return List::create(
    Named("beta_total") = beta_total,
    Named("beta_turnover") = beta_turnover,
    Named("beta_nestedness") = beta_nestedness,
    Named("distance") = distance
  );
}
