// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <set>
#include <cmath>
#include <limits>
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


// Rao's quadratic entropy: Q = sum_i sum_j p_i p_j d_ij, with
// p_i the relative abundance of species i and d_ij a pairwise distance
// (phylogenetic or functional). Absent species contribute zero.
// [[Rcpp::export]]
double calc_rao(NumericMatrix dist_mat, NumericVector abundances) {
  int n_species = dist_mat.nrow();

  double total = 0.0;
  for (int i = 0; i < n_species; i++) total += abundances[i];
  if (total <= 0.0) return 0.0;

  double q = 0.0;
  for (int i = 0; i < n_species; i++) {
    if (abundances[i] <= 0.0) continue;
    double pi = abundances[i] / total;
    for (int j = i + 1; j < n_species; j++) {
      if (abundances[j] <= 0.0) continue;
      q += pi * (abundances[j] / total) * dist_mat(i, j);
    }
  }
  return 2.0 * q;  // symmetric pairs (i, j) and (j, i); diagonal is zero
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
List cpp_phylo_knn_single(NumericMatrix species_pa,
                          NumericMatrix site_dist_mat,
                          NumericMatrix phylo_dist_mat,
                          int seed,
                          CharacterVector metrics,
                          Rcpp::Nullable<IntegerMatrix> tree_edge = R_NilValue,
                          Rcpp::Nullable<NumericVector> tree_edge_length = R_NilValue,
                          int tree_n_tips = 0) {
  int n_sites = species_pa.nrow();
  int n_species = species_pa.ncol();
  int n_metrics = metrics.size();

  // Check if PD is requested and tree data is available
  bool has_tree = tree_edge.isNotNull() && tree_edge_length.isNotNull();
  IntegerMatrix edge_mat;
  NumericVector edge_len;
  if (has_tree) {
    edge_mat = Rcpp::as<IntegerMatrix>(tree_edge);
    edge_len = Rcpp::as<NumericVector>(tree_edge_length);
  }

  // Output matrices: one row per metric
  NumericMatrix results(n_metrics, n_sites);

  std::vector<bool> visited(n_sites, false);
  std::vector<double> cumulative(n_species, 0.0);
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
    } else if (metric == "rao") {
      results(m, 0) = calc_rao(phylo_dist_mat, NumericVector(cumulative.begin(), cumulative.end()));
    } else if (metric == "pd" && has_tree) {
      results(m, 0) = calc_faith_pd(edge_mat, edge_len, tree_n_tips, species_present);
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
      } else if (metric == "rao") {
        results(m, step) = calc_rao(phylo_dist_mat, NumericVector(cumulative.begin(), cumulative.end()));
      } else if (metric == "pd" && has_tree) {
        results(m, step) = calc_faith_pd(edge_mat, edge_len, tree_n_tips, species_present);
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


// Worker struct for parallel phylo kNN
// Uses a single packed RMatrix<double> for output (n_seeds * n_metrics rows, n_sites cols)
// to match the working KnnWorker pattern exactly.
struct PhyloKnnWorker : public Worker {
  const RMatrix<double> species_pa;
  const RMatrix<double> site_dist_mat;
  const RMatrix<double> phylo_dist_mat;
  const RVector<int> seeds;
  RMatrix<double> results;  // packed: row = s * n_metrics + m

  bool do_mpd;
  bool do_mntd;
  bool do_pd;
  bool do_rao;
  int mpd_idx;
  int mntd_idx;
  int pd_idx;
  int rao_idx;
  int n_metrics;

  const bool has_tree;
  const RMatrix<int> edge_mat;
  const RVector<double> edge_len;
  const int tree_n_tips;

  PhyloKnnWorker(const NumericMatrix& species_pa_,
                 const NumericMatrix& site_dist_mat_,
                 const NumericMatrix& phylo_dist_mat_,
                 const IntegerVector& seeds_,
                 const std::vector<std::string>& metric_names_,
                 NumericMatrix& results_,
                 bool has_tree_,
                 const IntegerMatrix& edge_mat_,
                 const NumericVector& edge_len_,
                 int tree_n_tips_)
    : species_pa(species_pa_), site_dist_mat(site_dist_mat_),
      phylo_dist_mat(phylo_dist_mat_), seeds(seeds_),
      results(results_),
      do_mpd(false), do_mntd(false), do_pd(false), do_rao(false),
      mpd_idx(-1), mntd_idx(-1), pd_idx(-1), rao_idx(-1),
      n_metrics(metric_names_.size()),
      has_tree(has_tree_),
      edge_mat(edge_mat_), edge_len(edge_len_), tree_n_tips(tree_n_tips_) {
    for (size_t m = 0; m < metric_names_.size(); m++) {
      if (metric_names_[m] == "mpd") { do_mpd = true; mpd_idx = m; }
      else if (metric_names_[m] == "mntd") { do_mntd = true; mntd_idx = m; }
      else if (metric_names_[m] == "pd") { do_pd = true; pd_idx = m; }
      else if (metric_names_[m] == "rao") { do_rao = true; rao_idx = m; }
    }
  }

  void operator()(std::size_t begin, std::size_t end) {
    int n_sites = species_pa.nrow();
    int n_species = species_pa.ncol();
    const double INF = std::numeric_limits<double>::infinity();

    for (std::size_t s = begin; s < end; s++) {
      std::vector<bool> visited(n_sites, false);
      std::vector<bool> species_present(n_species, false);
      std::vector<double> cumulative(n_species, 0.0);

      int current = seeds[s];
      visited[current] = true;

      for (int sp = 0; sp < n_species; sp++) {
        cumulative[sp] += species_pa(current, sp);
        species_present[sp] = cumulative[sp] > 0;
      }

      if (do_mpd)  results(s * n_metrics + mpd_idx, 0) = calc_mpd_internal(species_present, n_species);
      if (do_mntd) results(s * n_metrics + mntd_idx, 0) = calc_mntd_internal(species_present, n_species);
      if (do_pd && has_tree) results(s * n_metrics + pd_idx, 0) = calc_pd_internal(species_present);
      if (do_rao) results(s * n_metrics + rao_idx, 0) = calc_rao_internal(cumulative, n_species);

      for (int step = 1; step < n_sites; step++) {
        double min_dist = INF;
        int next_site = -1;
        for (int j = 0; j < n_sites; j++) {
          if (!visited[j] && site_dist_mat(current, j) < min_dist) {
            min_dist = site_dist_mat(current, j);
            next_site = j;
          }
        }

        current = next_site;
        visited[current] = true;

        for (int sp = 0; sp < n_species; sp++) {
          cumulative[sp] += species_pa(current, sp);
          species_present[sp] = cumulative[sp] > 0;
        }

        if (do_mpd)  results(s * n_metrics + mpd_idx, step) = calc_mpd_internal(species_present, n_species);
        if (do_mntd) results(s * n_metrics + mntd_idx, step) = calc_mntd_internal(species_present, n_species);
        if (do_pd && has_tree) results(s * n_metrics + pd_idx, step) = calc_pd_internal(species_present);
        if (do_rao) results(s * n_metrics + rao_idx, step) = calc_rao_internal(cumulative, n_species);
      }
    }
  }

private:
  double calc_mpd_internal(const std::vector<bool>& species_present,
                           int n_species) const {
    std::vector<int> present;
    for (int i = 0; i < n_species; i++) {
      if (species_present[i]) present.push_back(i);
    }
    if (present.size() < 2) return 0.0;

    double sum_dist = 0.0;
    double count = 0.0;
    for (size_t i = 0; i < present.size(); i++) {
      for (size_t j = i + 1; j < present.size(); j++) {
        sum_dist += phylo_dist_mat(present[i], present[j]);
        count += 1.0;
      }
    }
    return count == 0 ? 0.0 : sum_dist / count;
  }

  double calc_mntd_internal(const std::vector<bool>& species_present,
                            int n_species) const {
    std::vector<int> present;
    for (int i = 0; i < n_species; i++) {
      if (species_present[i]) present.push_back(i);
    }
    if (present.size() < 2) return 0.0;

    double sum_nnd = 0.0;
    const double INF = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < present.size(); i++) {
      int sp_i = present[i];
      double min_d = INF;
      for (size_t j = 0; j < present.size(); j++) {
        if (i != j) {
          double d = phylo_dist_mat(sp_i, present[j]);
          if (d < min_d) min_d = d;
        }
      }
      sum_nnd += min_d;
    }
    return sum_nnd / present.size();
  }

  double calc_rao_internal(const std::vector<double>& abundances,
                           int n_species) const {
    double total = 0.0;
    for (int i = 0; i < n_species; i++) total += abundances[i];
    if (total <= 0.0) return 0.0;

    double q = 0.0;
    for (int i = 0; i < n_species; i++) {
      if (abundances[i] <= 0.0) continue;
      double pi = abundances[i] / total;
      for (int j = i + 1; j < n_species; j++) {
        if (abundances[j] <= 0.0) continue;
        q += pi * (abundances[j] / total) * phylo_dist_mat(i, j);
      }
    }
    return 2.0 * q;
  }

  double calc_pd_internal(const std::vector<bool>& species_present) const {
    int n_edges = edge_mat.nrow();
    std::vector<bool> edge_used(n_edges, false);

    for (int tip = 0; tip < tree_n_tips; tip++) {
      if (!species_present[tip]) continue;
      int current_node = tip + 1;

      while (true) {
        int found_edge = -1;
        for (int e = 0; e < n_edges; e++) {
          if (edge_mat(e, 1) == current_node) {
            found_edge = e;
            break;
          }
        }
        if (found_edge < 0) break;
        edge_used[found_edge] = true;
        current_node = edge_mat(found_edge, 0);
      }
    }

    double pd = 0.0;
    for (int e = 0; e < n_edges; e++) {
      if (edge_used[e]) pd += edge_len[e];
    }
    return pd;
  }
};


// [[Rcpp::export]]
List cpp_phylo_knn_parallel(NumericMatrix species_pa,
                            NumericMatrix site_dist_mat,
                            NumericMatrix phylo_dist_mat,
                            int n_seeds,
                            CharacterVector metrics,
                            int n_cores = 1,
                            bool progress = false,
                            Rcpp::Nullable<IntegerMatrix> tree_edge = R_NilValue,
                            Rcpp::Nullable<NumericVector> tree_edge_length = R_NilValue,
                            int tree_n_tips = 0) {
  int n_sites = species_pa.nrow();
  int n_metrics = metrics.size();

  IntegerVector seeds = Rcpp::sample(n_sites, n_seeds, true) - 1;

  // Check if tree data is available
  bool has_tree = tree_edge.isNotNull() && tree_edge_length.isNotNull();
  IntegerMatrix edge_mat;
  NumericVector edge_len;
  if (has_tree) {
    edge_mat = Rcpp::as<IntegerMatrix>(tree_edge);
    edge_len = Rcpp::as<NumericVector>(tree_edge_length);
  } else {
    // Create dummy matrices for worker
    edge_mat = IntegerMatrix(1, 2);
    edge_len = NumericVector(1);
  }

  // Convert metrics to std::vector<std::string>
  std::vector<std::string> metric_names(n_metrics);
  for (int m = 0; m < n_metrics; m++) {
    metric_names[m] = Rcpp::as<std::string>(metrics[m]);
  }

  // Single packed output matrix: (n_seeds * n_metrics) rows x n_sites cols
  NumericMatrix packed_results(n_seeds * n_metrics, n_sites);

  if (n_cores > 1) {
    PhyloKnnWorker worker(species_pa, site_dist_mat, phylo_dist_mat,
                          seeds, metric_names, packed_results,
                          has_tree, edge_mat, edge_len, tree_n_tips);
    parallelFor(0, n_seeds, worker, 1, n_cores);
  } else {
    // Sequential fallback
    for (int s = 0; s < n_seeds; s++) {
      List single = cpp_phylo_knn_single(species_pa, site_dist_mat,
                                          phylo_dist_mat, seeds[s], metrics,
                                          tree_edge, tree_edge_length, tree_n_tips);
      for (int m = 0; m < n_metrics; m++) {
        NumericVector curve = single[m];
        for (int st = 0; st < n_sites; st++) {
          packed_results(s * n_metrics + m, st) = curve[st];
        }
      }
    }
  }

  // Unpack into per-metric matrices
  List out(n_metrics);
  for (int m = 0; m < n_metrics; m++) {
    NumericMatrix mat(n_seeds, n_sites);
    for (int s = 0; s < n_seeds; s++) {
      for (int st = 0; st < n_sites; st++) {
        mat(s, st) = packed_results(s * n_metrics + m, st);
      }
    }
    out[m] = mat;
  }
  out.attr("names") = metrics;

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


// Functional Rao's quadratic entropy using Euclidean trait distance:
// Q = sum_i sum_j p_i p_j ||t_i - t_j||, with p the relative abundances.
// [[Rcpp::export]]
double calc_rao_traits(NumericMatrix traits, NumericVector abundances) {
  int n_species = traits.nrow();
  int n_traits = traits.ncol();

  double total = 0.0;
  for (int i = 0; i < n_species; i++) total += abundances[i];
  if (total <= 0.0) return 0.0;

  double q = 0.0;
  for (int i = 0; i < n_species; i++) {
    if (abundances[i] <= 0.0) continue;
    double pi = abundances[i] / total;
    for (int j = i + 1; j < n_species; j++) {
      if (abundances[j] <= 0.0) continue;
      double d2 = 0.0;
      for (int t = 0; t < n_traits; t++) {
        double diff = traits(i, t) - traits(j, t);
        d2 += diff * diff;
      }
      q += pi * (abundances[j] / total) * std::sqrt(d2);
    }
  }
  return 2.0 * q;
}


// [[Rcpp::export]]
List cpp_func_knn_single(NumericMatrix species_mat,
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
    } else if (metric == "rao") {
      results(m, 0) = calc_rao_traits(traits, abundances);
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
      } else if (metric == "rao") {
        results(m, step) = calc_rao_traits(traits, abundances);
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


// Worker struct for parallel functional kNN
// Uses a single packed RMatrix<double> for output to match KnnWorker pattern.
struct FuncKnnWorker : public Worker {
  const RMatrix<double> species_mat;
  const RMatrix<double> site_dist_mat;
  const RMatrix<double> traits;
  const RVector<int> seeds;
  RMatrix<double> results;  // packed: row = s * n_metrics + m

  bool do_fdis;
  bool do_fric;
  bool do_rao;
  int fdis_idx;
  int fric_idx;
  int rao_idx;
  int n_metrics;

  FuncKnnWorker(const NumericMatrix& species_mat_,
                const NumericMatrix& site_dist_mat_,
                const NumericMatrix& traits_,
                const IntegerVector& seeds_,
                const std::vector<std::string>& metric_names_,
                NumericMatrix& results_)
    : species_mat(species_mat_), site_dist_mat(site_dist_mat_),
      traits(traits_), seeds(seeds_), results(results_),
      do_fdis(false), do_fric(false), do_rao(false),
      fdis_idx(-1), fric_idx(-1), rao_idx(-1),
      n_metrics(metric_names_.size()) {
    for (size_t m = 0; m < metric_names_.size(); m++) {
      if (metric_names_[m] == "fdis") { do_fdis = true; fdis_idx = m; }
      else if (metric_names_[m] == "fric") { do_fric = true; fric_idx = m; }
      else if (metric_names_[m] == "rao") { do_rao = true; rao_idx = m; }
    }
  }

  void operator()(std::size_t begin, std::size_t end) {
    int n_sites = species_mat.nrow();
    int n_species = species_mat.ncol();
    int n_traits = traits.ncol();
    const double INF = std::numeric_limits<double>::infinity();

    for (std::size_t s = begin; s < end; s++) {
      std::vector<bool> visited(n_sites, false);
      std::vector<double> cumulative(n_species, 0.0);
      std::vector<bool> species_present(n_species, false);
      std::vector<double> abundances(n_species, 0.0);

      int current = seeds[s];
      visited[current] = true;

      for (int sp = 0; sp < n_species; sp++) {
        cumulative[sp] += species_mat(current, sp);
        species_present[sp] = cumulative[sp] > 0;
        abundances[sp] = cumulative[sp];
      }

      if (do_fdis) results(s * n_metrics + fdis_idx, 0) = calc_fdis_internal(species_present, abundances, n_species, n_traits);
      if (do_fric) results(s * n_metrics + fric_idx, 0) = calc_fric_internal(species_present, n_species, n_traits);
      if (do_rao) results(s * n_metrics + rao_idx, 0) = calc_rao_traits_internal(abundances, n_species, n_traits);

      for (int step = 1; step < n_sites; step++) {
        double min_dist = INF;
        int next_site = -1;
        for (int j = 0; j < n_sites; j++) {
          if (!visited[j] && site_dist_mat(current, j) < min_dist) {
            min_dist = site_dist_mat(current, j);
            next_site = j;
          }
        }

        current = next_site;
        visited[current] = true;

        for (int sp = 0; sp < n_species; sp++) {
          cumulative[sp] += species_mat(current, sp);
          species_present[sp] = cumulative[sp] > 0;
          abundances[sp] = cumulative[sp];
        }

        if (do_fdis) results(s * n_metrics + fdis_idx, step) = calc_fdis_internal(species_present, abundances, n_species, n_traits);
        if (do_fric) results(s * n_metrics + fric_idx, step) = calc_fric_internal(species_present, n_species, n_traits);
        if (do_rao) results(s * n_metrics + rao_idx, step) = calc_rao_traits_internal(abundances, n_species, n_traits);
      }
    }
  }

private:
  double calc_fdis_internal(const std::vector<bool>& species_present,
                            const std::vector<double>& abundances,
                            int n_species, int n_traits) const {
    std::vector<int> present;
    double total_abund = 0.0;
    for (int i = 0; i < n_species; i++) {
      if (species_present[i]) {
        present.push_back(i);
        total_abund += abundances[i];
      }
    }
    if (present.size() < 2 || total_abund == 0) return 0.0;

    std::vector<double> centroid(n_traits, 0.0);
    for (size_t i = 0; i < present.size(); i++) {
      int sp = present[i];
      double w = abundances[sp] / total_abund;
      for (int t = 0; t < n_traits; t++) {
        centroid[t] += traits(sp, t) * w;
      }
    }

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

  double calc_rao_traits_internal(const std::vector<double>& abundances,
                                  int n_species, int n_traits) const {
    double total = 0.0;
    for (int i = 0; i < n_species; i++) total += abundances[i];
    if (total <= 0.0) return 0.0;

    double q = 0.0;
    for (int i = 0; i < n_species; i++) {
      if (abundances[i] <= 0.0) continue;
      double pi = abundances[i] / total;
      for (int j = i + 1; j < n_species; j++) {
        if (abundances[j] <= 0.0) continue;
        double d2 = 0.0;
        for (int t = 0; t < n_traits; t++) {
          double diff = traits(i, t) - traits(j, t);
          d2 += diff * diff;
        }
        q += pi * (abundances[j] / total) * std::sqrt(d2);
      }
    }
    return 2.0 * q;
  }

  double calc_fric_internal(const std::vector<bool>& species_present,
                            int n_species, int n_traits) const {
    std::vector<int> present;
    for (int i = 0; i < n_species; i++) {
      if (species_present[i]) present.push_back(i);
    }
    if (present.size() <= (size_t)n_traits) return 0.0;

    const double INF = std::numeric_limits<double>::infinity();
    double volume = 1.0;
    for (int t = 0; t < n_traits; t++) {
      double min_val = INF;
      double max_val = -INF;
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
};


// [[Rcpp::export]]
List cpp_func_knn_parallel(NumericMatrix species_mat,
                           NumericMatrix site_dist_mat,
                           NumericMatrix traits,
                           int n_seeds,
                           CharacterVector metrics,
                           int n_cores = 1,
                           bool progress = false) {
  int n_sites = species_mat.nrow();
  int n_metrics = metrics.size();

  IntegerVector seeds = Rcpp::sample(n_sites, n_seeds, true) - 1;

  std::vector<std::string> metric_names(n_metrics);
  for (int m = 0; m < n_metrics; m++) {
    metric_names[m] = Rcpp::as<std::string>(metrics[m]);
  }

  // Single packed output matrix: (n_seeds * n_metrics) rows x n_sites cols
  NumericMatrix packed_results(n_seeds * n_metrics, n_sites);

  if (n_cores > 1) {
    FuncKnnWorker worker(species_mat, site_dist_mat, traits,
                         seeds, metric_names, packed_results);
    parallelFor(0, n_seeds, worker, 1, n_cores);
  } else {
    for (int s = 0; s < n_seeds; s++) {
      List single = cpp_func_knn_single(species_mat, site_dist_mat,
                                         traits, seeds[s], metrics);
      for (int m = 0; m < n_metrics; m++) {
        NumericVector curve = single[m];
        for (int st = 0; st < n_sites; st++) {
          packed_results(s * n_metrics + m, st) = curve[st];
        }
      }
    }
  }

  // Unpack into per-metric matrices
  List out(n_metrics);
  for (int m = 0; m < n_metrics; m++) {
    NumericMatrix mat(n_seeds, n_sites);
    for (int s = 0; s < n_seeds; s++) {
      for (int st = 0; st < n_sites; st++) {
        mat(s, st) = packed_results(s * n_metrics + m, st);
      }
    }
    out[m] = mat;
  }
  out.attr("names") = metrics;

  return out;
}
