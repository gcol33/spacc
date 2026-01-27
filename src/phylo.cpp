// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <set>
#include <cmath>
using namespace Rcpp;
using namespace RcppParallel;


// ============================================================================
// PHYLOGENETIC DIVERSITY CALCULATIONS
// ============================================================================

// [[Rcpp::export]]
double calc_mpd(NumericMatrix dist_mat,
                LogicalVector species_present,
                bool abundance_weighted = false,
                NumericVector abundances = NumericVector()) {
  int n_species = dist_mat.nrow();

  // Get indices of present species
  std::vector<int> present;
  for (int i = 0; i < n_species; i++) {
    if (species_present[i]) {
      present.push_back(i);
    }
  }

  if (present.size() < 2) return 0.0;

  double sum_dist = 0.0;
  double sum_weight = 0.0;

  if (abundance_weighted && abundances.size() == n_species) {
    // Abundance-weighted MPD
    for (size_t i = 0; i < present.size(); i++) {
      for (size_t j = i + 1; j < present.size(); j++) {
        int sp_i = present[i];
        int sp_j = present[j];
        double w = abundances[sp_i] * abundances[sp_j];
        sum_dist += dist_mat(sp_i, sp_j) * w;
        sum_weight += w;
      }
    }
  } else {
    // Unweighted MPD
    for (size_t i = 0; i < present.size(); i++) {
      for (size_t j = i + 1; j < present.size(); j++) {
        sum_dist += dist_mat(present[i], present[j]);
        sum_weight += 1.0;
      }
    }
  }

  if (sum_weight == 0) return 0.0;
  return sum_dist / sum_weight;
}


// [[Rcpp::export]]
double calc_mntd(NumericMatrix dist_mat,
                 LogicalVector species_present,
                 bool abundance_weighted = false,
                 NumericVector abundances = NumericVector()) {
  int n_species = dist_mat.nrow();

  std::vector<int> present;
  for (int i = 0; i < n_species; i++) {
    if (species_present[i]) {
      present.push_back(i);
    }
  }

  if (present.size() < 2) return 0.0;

  double sum_nnd = 0.0;
  double sum_weight = 0.0;

  for (size_t i = 0; i < present.size(); i++) {
    int sp_i = present[i];

    // Find nearest neighbor
    double min_dist = R_PosInf;
    for (size_t j = 0; j < present.size(); j++) {
      if (i != j) {
        double d = dist_mat(sp_i, present[j]);
        if (d < min_dist) min_dist = d;
      }
    }

    double w = 1.0;
    if (abundance_weighted && abundances.size() == n_species) {
      w = abundances[sp_i];
    }

    sum_nnd += min_dist * w;
    sum_weight += w;
  }

  if (sum_weight == 0) return 0.0;
  return sum_nnd / sum_weight;
}


// [[Rcpp::export]]
double calc_faith_pd(IntegerMatrix edge,
                     NumericVector edge_length,
                     int n_tips,
                     LogicalVector species_present) {
  int n_edges = edge.nrow();

  // Mark which edges are used
  std::vector<bool> edge_used(n_edges, false);

  // For each present species, trace path to root and mark edges
  for (int tip = 0; tip < n_tips; tip++) {
    if (!species_present[tip]) continue;

    int current_node = tip + 1;  // R uses 1-based indexing

    // Trace up to root
    while (true) {
      // Find edge leading to current_node
      int found_edge = -1;
      for (int e = 0; e < n_edges; e++) {
        if (edge(e, 1) == current_node) {
          found_edge = e;
          break;
        }
      }

      if (found_edge < 0) break;  // Reached root

      edge_used[found_edge] = true;
      current_node = edge(found_edge, 0);  // Move to parent
    }
  }

  // Sum used edge lengths
  double pd = 0.0;
  for (int e = 0; e < n_edges; e++) {
    if (edge_used[e]) {
      pd += edge_length[e];
    }
  }

  return pd;
}


// [[Rcpp::export]]
List cpp_phylo_knn_single(IntegerMatrix species_pa,
                          NumericMatrix site_dist_mat,
                          NumericMatrix phylo_dist_mat,
                          int seed,
                          CharacterVector metrics) {
  int n_sites = species_pa.nrow();
  int n_species = species_pa.ncol();
  int n_metrics = metrics.size();

  // Output matrices: one row per metric
  NumericMatrix results(n_metrics, n_sites);

  std::vector<bool> visited(n_sites, false);
  std::vector<int> cumulative(n_species, 0);
  LogicalVector species_present(n_species);

  int current = seed;
  visited[current] = true;

  // Add first site
  for (int sp = 0; sp < n_species; sp++) {
    cumulative[sp] += species_pa(current, sp);
    species_present[sp] = cumulative[sp] > 0;
  }

  // Calculate metrics for first site
  for (int m = 0; m < n_metrics; m++) {
    String metric = metrics[m];
    if (metric == "mpd") {
      results(m, 0) = calc_mpd(phylo_dist_mat, species_present, false, NumericVector());
    } else if (metric == "mntd") {
      results(m, 0) = calc_mntd(phylo_dist_mat, species_present, false, NumericVector());
    } else {
      results(m, 0) = NA_REAL;
    }
  }

  for (int step = 1; step < n_sites; step++) {
    // Find nearest unvisited
    double min_dist = R_PosInf;
    int next = -1;
    for (int j = 0; j < n_sites; j++) {
      if (!visited[j] && site_dist_mat(current, j) < min_dist) {
        min_dist = site_dist_mat(current, j);
        next = j;
      }
    }

    current = next;
    visited[current] = true;

    // Accumulate
    for (int sp = 0; sp < n_species; sp++) {
      cumulative[sp] += species_pa(current, sp);
      species_present[sp] = cumulative[sp] > 0;
    }

    // Calculate metrics
    for (int m = 0; m < n_metrics; m++) {
      String metric = metrics[m];
      if (metric == "mpd") {
        results(m, step) = calc_mpd(phylo_dist_mat, species_present, false, NumericVector());
      } else if (metric == "mntd") {
        results(m, step) = calc_mntd(phylo_dist_mat, species_present, false, NumericVector());
      } else {
        results(m, step) = NA_REAL;
      }
    }
  }

  // Convert to named list
  List out(n_metrics);
  for (int m = 0; m < n_metrics; m++) {
    NumericVector curve(n_sites);
    for (int s = 0; s < n_sites; s++) {
      curve[s] = results(m, s);
    }
    out[m] = curve;
  }
  out.attr("names") = metrics;

  return out;
}


// [[Rcpp::export]]
List cpp_phylo_knn_parallel(IntegerMatrix species_pa,
                            NumericMatrix site_dist_mat,
                            NumericMatrix phylo_dist_mat,
                            int n_seeds,
                            CharacterVector metrics,
                            int n_cores = 1,
                            bool progress = false) {
  int n_sites = species_pa.nrow();
  int n_metrics = metrics.size();

  IntegerVector seeds = Rcpp::sample(n_sites, n_seeds, true) - 1;

  // Allocate output: list of matrices, one per metric
  List out(n_metrics);
  for (int m = 0; m < n_metrics; m++) {
    out[m] = NumericMatrix(n_seeds, n_sites);
  }
  out.attr("names") = metrics;

  // Sequential for now (parallelization would require more complex setup)
  for (int s = 0; s < n_seeds; s++) {
    List single = cpp_phylo_knn_single(species_pa, site_dist_mat,
                                        phylo_dist_mat, seeds[s], metrics);

    for (int m = 0; m < n_metrics; m++) {
      NumericMatrix mat = out[m];
      NumericVector curve = single[m];
      for (int st = 0; st < n_sites; st++) {
        mat(s, st) = curve[st];
      }
      out[m] = mat;
    }
  }

  return out;
}


// ============================================================================
// FUNCTIONAL DIVERSITY CALCULATIONS
// ============================================================================

// [[Rcpp::export]]
double calc_fdis(NumericMatrix traits,
                 LogicalVector species_present,
                 NumericVector abundances) {
  int n_species = traits.nrow();
  int n_traits = traits.ncol();

  // Get present species
  std::vector<int> present;
  double total_abund = 0.0;

  for (int i = 0; i < n_species; i++) {
    if (species_present[i]) {
      present.push_back(i);
      total_abund += abundances[i];
    }
  }

  if (present.size() < 2 || total_abund == 0) return 0.0;

  // Calculate weighted centroid
  std::vector<double> centroid(n_traits, 0.0);
  for (size_t i = 0; i < present.size(); i++) {
    int sp = present[i];
    double w = abundances[sp] / total_abund;
    for (int t = 0; t < n_traits; t++) {
      centroid[t] += traits(sp, t) * w;
    }
  }

  // Calculate weighted mean distance to centroid
  double sum_dist = 0.0;
  for (size_t i = 0; i < present.size(); i++) {
    int sp = present[i];
    double dist_sq = 0.0;
    for (int t = 0; t < n_traits; t++) {
      double diff = traits(sp, t) - centroid[t];
      dist_sq += diff * diff;
    }
    sum_dist += std::sqrt(dist_sq) * (abundances[sp] / total_abund);
  }

  return sum_dist;
}


// [[Rcpp::export]]
double calc_fric_approx(NumericMatrix traits,
                        LogicalVector species_present) {
  int n_species = traits.nrow();
  int n_traits = traits.ncol();

  std::vector<int> present;
  for (int i = 0; i < n_species; i++) {
    if (species_present[i]) {
      present.push_back(i);
    }
  }

  if (present.size() <= (size_t)n_traits) return 0.0;

  // Calculate range for each trait and multiply (hypervolume approximation)
  double volume = 1.0;
  for (int t = 0; t < n_traits; t++) {
    double min_val = R_PosInf;
    double max_val = R_NegInf;

    for (size_t i = 0; i < present.size(); i++) {
      double val = traits(present[i], t);
      if (val < min_val) min_val = val;
      if (val > max_val) max_val = val;
    }

    double range = max_val - min_val;
    if (range > 0) {
      volume *= range;
    } else {
      volume = 0.0;
      break;
    }
  }

  return volume;
}


// [[Rcpp::export]]
List cpp_func_knn_single(IntegerMatrix species_mat,
                         NumericMatrix site_dist_mat,
                         NumericMatrix traits,
                         int seed,
                         CharacterVector metrics) {
  int n_sites = species_mat.nrow();
  int n_species = species_mat.ncol();
  int n_metrics = metrics.size();

  NumericMatrix results(n_metrics, n_sites);

  std::vector<bool> visited(n_sites, false);
  std::vector<double> cumulative(n_species, 0.0);
  LogicalVector species_present(n_species);
  NumericVector abundances(n_species);

  int current = seed;
  visited[current] = true;

  for (int sp = 0; sp < n_species; sp++) {
    cumulative[sp] += species_mat(current, sp);
    species_present[sp] = cumulative[sp] > 0;
    abundances[sp] = cumulative[sp];
  }

  for (int m = 0; m < n_metrics; m++) {
    String metric = metrics[m];
    if (metric == "fdis") {
      results(m, 0) = calc_fdis(traits, species_present, abundances);
    } else if (metric == "fric") {
      results(m, 0) = calc_fric_approx(traits, species_present);
    } else {
      results(m, 0) = NA_REAL;
    }
  }

  for (int step = 1; step < n_sites; step++) {
    double min_dist = R_PosInf;
    int next = -1;
    for (int j = 0; j < n_sites; j++) {
      if (!visited[j] && site_dist_mat(current, j) < min_dist) {
        min_dist = site_dist_mat(current, j);
        next = j;
      }
    }

    current = next;
    visited[current] = true;

    for (int sp = 0; sp < n_species; sp++) {
      cumulative[sp] += species_mat(current, sp);
      species_present[sp] = cumulative[sp] > 0;
      abundances[sp] = cumulative[sp];
    }

    for (int m = 0; m < n_metrics; m++) {
      String metric = metrics[m];
      if (metric == "fdis") {
        results(m, step) = calc_fdis(traits, species_present, abundances);
      } else if (metric == "fric") {
        results(m, step) = calc_fric_approx(traits, species_present);
      } else {
        results(m, step) = NA_REAL;
      }
    }
  }

  List out(n_metrics);
  for (int m = 0; m < n_metrics; m++) {
    NumericVector curve(n_sites);
    for (int s = 0; s < n_sites; s++) {
      curve[s] = results(m, s);
    }
    out[m] = curve;
  }
  out.attr("names") = metrics;

  return out;
}


// [[Rcpp::export]]
List cpp_func_knn_parallel(IntegerMatrix species_mat,
                           NumericMatrix site_dist_mat,
                           NumericMatrix traits,
                           int n_seeds,
                           CharacterVector metrics,
                           int n_cores = 1,
                           bool progress = false) {
  int n_sites = species_mat.nrow();
  int n_metrics = metrics.size();

  IntegerVector seeds = Rcpp::sample(n_sites, n_seeds, true) - 1;

  List out(n_metrics);
  for (int m = 0; m < n_metrics; m++) {
    out[m] = NumericMatrix(n_seeds, n_sites);
  }
  out.attr("names") = metrics;

  for (int s = 0; s < n_seeds; s++) {
    List single = cpp_func_knn_single(species_mat, site_dist_mat,
                                       traits, seeds[s], metrics);

    for (int m = 0; m < n_metrics; m++) {
      NumericMatrix mat = out[m];
      NumericVector curve = single[m];
      for (int st = 0; st < n_sites; st++) {
        mat(s, st) = curve[st];
      }
      out[m] = mat;
    }
  }

  return out;
}
