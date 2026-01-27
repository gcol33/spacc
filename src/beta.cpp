// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <set>
#include <algorithm>
using namespace Rcpp;
using namespace RcppParallel;


// ============================================================================
// BETA DIVERSITY CALCULATIONS (Baselga 2010)
// ============================================================================

// Structure to hold beta diversity components
struct BetaComponents {
  double total;
  double turnover;
  double nestedness;
};


// Calculate Sorensen-based beta diversity and its components
// a = shared species, b = unique to set1, c = unique to set2
BetaComponents calc_beta_sorensen(int a, int b, int c) {
  BetaComponents result = {0.0, 0.0, 0.0};

  if (a + b + c == 0) return result;

  // Total beta (Sorensen dissimilarity)
  result.total = (double)(b + c) / (2.0 * a + b + c);

  // Turnover component (Simpson)
  int min_bc = std::min(b, c);
  if (a + min_bc > 0) {
    result.turnover = (double)min_bc / (a + min_bc);
  }

  // Nestedness component (total - turnover)
  result.nestedness = result.total - result.turnover;

  return result;
}


// Calculate Jaccard-based beta diversity and its components
BetaComponents calc_beta_jaccard(int a, int b, int c) {
  BetaComponents result = {0.0, 0.0, 0.0};

  if (a + b + c == 0) return result;

  // Total beta (Jaccard dissimilarity)
  result.total = (double)(b + c) / (a + b + c);

  // Turnover component
  int min_bc = std::min(b, c);
  if (a + min_bc > 0) {
    result.turnover = 2.0 * (double)min_bc / (a + 2.0 * min_bc);
  }

  // Nestedness component
  result.nestedness = result.total - result.turnover;

  return result;
}


// Count shared and unique species between two sets
void count_abc(const std::set<int>& set1, const std::set<int>& set2,
               int& a, int& b, int& c) {
  a = 0; b = 0; c = 0;

  for (int sp : set1) {
    if (set2.count(sp) > 0) {
      a++;
    } else {
      b++;
    }
  }

  for (int sp : set2) {
    if (set1.count(sp) == 0) {
      c++;
    }
  }
}


// [[Rcpp::export]]
List cpp_beta_knn_single(IntegerMatrix species_pa,
                         NumericMatrix dist_mat,
                         int seed,
                         bool use_jaccard = false) {
  int n_sites = species_pa.nrow();
  int n_species = species_pa.ncol();

  // Output vectors (n_sites - 1 comparisons)
  NumericVector beta_total(n_sites - 1);
  NumericVector beta_turn(n_sites - 1);
  NumericVector beta_nest(n_sites - 1);
  NumericVector cum_distance(n_sites - 1);
  IntegerVector n_species_acc(n_sites);

  std::vector<bool> visited(n_sites, false);
  std::set<int> accumulated_species;
  std::set<int> new_site_species;

  int current = seed;
  visited[current] = true;

  // Initialize with first site
  for (int sp = 0; sp < n_species; sp++) {
    if (species_pa(current, sp) > 0) {
      accumulated_species.insert(sp);
    }
  }
  n_species_acc[0] = accumulated_species.size();

  double total_dist = 0.0;

  for (int step = 0; step < n_sites - 1; step++) {
    // Find nearest unvisited
    double min_dist = R_PosInf;
    int next = -1;
    for (int j = 0; j < n_sites; j++) {
      if (!visited[j] && dist_mat(current, j) < min_dist) {
        min_dist = dist_mat(current, j);
        next = j;
      }
    }

    total_dist += min_dist;
    cum_distance[step] = total_dist;

    // Get species at new site
    new_site_species.clear();
    for (int sp = 0; sp < n_species; sp++) {
      if (species_pa(next, sp) > 0) {
        new_site_species.insert(sp);
      }
    }

    // Calculate beta between accumulated and new site
    int a, b, c;
    count_abc(accumulated_species, new_site_species, a, b, c);

    BetaComponents beta;
    if (use_jaccard) {
      beta = calc_beta_jaccard(a, b, c);
    } else {
      beta = calc_beta_sorensen(a, b, c);
    }

    beta_total[step] = beta.total;
    beta_turn[step] = beta.turnover;
    beta_nest[step] = beta.nestedness;

    // Add new species to accumulated set
    for (int sp : new_site_species) {
      accumulated_species.insert(sp);
    }
    n_species_acc[step + 1] = accumulated_species.size();

    current = next;
    visited[current] = true;
  }

  return List::create(
    Named("beta_total") = beta_total,
    Named("beta_turnover") = beta_turn,
    Named("beta_nestedness") = beta_nest,
    Named("distance") = cum_distance,
    Named("richness") = n_species_acc
  );
}


// Worker for parallel beta accumulation
struct BetaKnnWorker : public Worker {
  const RMatrix<int> species_pa;
  const RMatrix<double> dist_mat;
  const RVector<int> seeds;
  const bool use_jaccard;

  RMatrix<double> beta_total;
  RMatrix<double> beta_turn;
  RMatrix<double> beta_nest;
  RMatrix<double> distances;

  BetaKnnWorker(const IntegerMatrix& sp, const NumericMatrix& dm,
                const IntegerVector& s, bool jac,
                NumericMatrix& bt, NumericMatrix& btu,
                NumericMatrix& bn, NumericMatrix& d)
    : species_pa(sp), dist_mat(dm), seeds(s), use_jaccard(jac),
      beta_total(bt), beta_turn(btu), beta_nest(bn), distances(d) {}

  void operator()(std::size_t begin, std::size_t end) {
    int n_sites = species_pa.nrow();
    int n_species = species_pa.ncol();

    for (std::size_t s = begin; s < end; s++) {
      std::vector<bool> visited(n_sites, false);
      std::set<int> accumulated;
      std::set<int> new_species;

      int current = seeds[s];
      visited[current] = true;

      for (int sp = 0; sp < n_species; sp++) {
        if (species_pa(current, sp) > 0) {
          accumulated.insert(sp);
        }
      }

      double total_dist = 0.0;

      for (int step = 0; step < n_sites - 1; step++) {
        double min_dist = R_PosInf;
        int next = -1;
        for (int j = 0; j < n_sites; j++) {
          if (!visited[j] && dist_mat(current, j) < min_dist) {
            min_dist = dist_mat(current, j);
            next = j;
          }
        }

        total_dist += min_dist;
        distances(s, step) = total_dist;

        new_species.clear();
        for (int sp = 0; sp < n_species; sp++) {
          if (species_pa(next, sp) > 0) {
            new_species.insert(sp);
          }
        }

        int a = 0, b = 0, c = 0;
        for (int sp : accumulated) {
          if (new_species.count(sp) > 0) a++;
          else b++;
        }
        for (int sp : new_species) {
          if (accumulated.count(sp) == 0) c++;
        }

        BetaComponents beta;
        if (use_jaccard) {
          beta = calc_beta_jaccard(a, b, c);
        } else {
          beta = calc_beta_sorensen(a, b, c);
        }

        beta_total(s, step) = beta.total;
        beta_turn(s, step) = beta.turnover;
        beta_nest(s, step) = beta.nestedness;

        for (int sp : new_species) {
          accumulated.insert(sp);
        }

        current = next;
        visited[current] = true;
      }
    }
  }
};


// [[Rcpp::export]]
List cpp_beta_knn_parallel(IntegerMatrix species_pa,
                           NumericMatrix dist_mat,
                           int n_seeds,
                           bool use_jaccard = false,
                           int n_cores = 1,
                           bool progress = false) {
  int n_sites = species_pa.nrow();
  int n_steps = n_sites - 1;

  IntegerVector seeds = Rcpp::sample(n_sites, n_seeds, true) - 1;

  NumericMatrix beta_total(n_seeds, n_steps);
  NumericMatrix beta_turn(n_seeds, n_steps);
  NumericMatrix beta_nest(n_seeds, n_steps);
  NumericMatrix distances(n_seeds, n_steps);

  if (n_cores > 1) {
    BetaKnnWorker worker(species_pa, dist_mat, seeds, use_jaccard,
                         beta_total, beta_turn, beta_nest, distances);
    parallelFor(0, n_seeds, worker);
  } else {
    for (int s = 0; s < n_seeds; s++) {
      List single = cpp_beta_knn_single(species_pa, dist_mat,
                                        seeds[s], use_jaccard);
      NumericVector bt = single["beta_total"];
      NumericVector btu = single["beta_turnover"];
      NumericVector bn = single["beta_nestedness"];
      NumericVector d = single["distance"];

      for (int step = 0; step < n_steps; step++) {
        beta_total(s, step) = bt[step];
        beta_turn(s, step) = btu[step];
        beta_nest(s, step) = bn[step];
        distances(s, step) = d[step];
      }
    }
  }

  return List::create(
    Named("beta_total") = beta_total,
    Named("beta_turnover") = beta_turn,
    Named("beta_nestedness") = beta_nest,
    Named("distance") = distances
  );
}
