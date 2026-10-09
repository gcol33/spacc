// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <algorithm>
#include <memory>
#include <string>
#include <vector>
#include "core/phylo_core.h"
#include "core/func_core.h"
#include "parallel_utils.h"
using namespace Rcpp;
using namespace RcppParallel;


// ---------------------------------------------------------------------------
// Helpers shared by the exported single-community functions and the workers
// ---------------------------------------------------------------------------

static std::vector<int> present_from_flags(const LogicalVector& flags) {
  std::vector<int> present;
  for (int i = 0; i < flags.size(); i++) {
    if (flags[i]) present.push_back(i);
  }
  return present;
}

// Position of `key` in `names`, or -1.
static int index_of(const std::vector<std::string>& names, const char* key) {
  for (std::size_t i = 0; i < names.size(); i++) {
    if (names[i] == key) return static_cast<int>(i);
  }
  return -1;
}

static std::vector<std::string> metric_names(const CharacterVector& metrics) {
  std::vector<std::string> names(metrics.size());
  for (int m = 0; m < metrics.size(); m++) names[m] = Rcpp::as<std::string>(metrics[m]);
  return names;
}

static spacc::PhyloTree build_tree(const IntegerMatrix& edge, const NumericVector& edge_length,
                                   int n_tips) {
  return spacc::PhyloTree(
    edge.nrow(),
    [&](int e) { return edge(e, 0); },
    [&](int e) { return edge(e, 1); },
    [&](int e) { return edge_length[e]; },
    n_tips);
}

// Unpack a flat (n_seeds * n_metrics) x n_sites matrix into a named list of
// n_seeds x n_sites matrices, one per metric.
static List run_metric_list(const NumericMatrix& flat, int n_seeds,
                            const CharacterVector& metrics) {
  return split_components(flat, n_seeds, metrics.size(), metrics);
}


// ===========================================================================
// PHYLOGENETIC DIVERSITY
// ===========================================================================

// [[Rcpp::export]]
double calc_mpd(NumericMatrix dist_mat, LogicalVector species_present,
                bool abundance_weighted = false,
                NumericVector abundances = NumericVector()) {
  const bool weighted = abundance_weighted && abundances.size() == dist_mat.nrow();
  return spacc::mean_pairwise_distance(
    [&](int i, int j) { return dist_mat(i, j); },
    present_from_flags(species_present),
    [&](int i) { return weighted ? abundances[i] : 1.0; });
}


// [[Rcpp::export]]
double calc_mntd(NumericMatrix dist_mat, LogicalVector species_present,
                 bool abundance_weighted = false,
                 NumericVector abundances = NumericVector()) {
  const bool weighted = abundance_weighted && abundances.size() == dist_mat.nrow();
  return spacc::mean_nearest_taxon_distance(
    [&](int i, int j) { return dist_mat(i, j); },
    present_from_flags(species_present),
    [&](int i) { return weighted ? abundances[i] : 1.0; });
}


// Rao's quadratic entropy on a pairwise (phylogenetic or functional) distance
// matrix with relative abundances; absent species contribute zero.
// [[Rcpp::export]]
double calc_rao(NumericMatrix dist_mat, NumericVector abundances) {
  return spacc::rao_quadratic_entropy(
    [&](int i, int j) { return dist_mat(i, j); },
    [&](int i) { return abundances[i]; },
    dist_mat.nrow());
}


// [[Rcpp::export]]
double calc_faith_pd(IntegerMatrix edge, NumericVector edge_length, int n_tips,
                     LogicalVector species_present) {
  spacc::PhyloTree tree = build_tree(edge, edge_length, n_tips);
  return tree.faith_pd([&](int tip) { return species_present[tip] != 0; });
}


struct PhyloMetric {
  RMatrix<double> phylo_dist;
  RMatrix<double> out;       // (n_seeds * n_metrics) x n_sites
  int n_metrics, mpd_idx, mntd_idx, pd_idx, rao_idx;
  std::shared_ptr<const spacc::PhyloTree> tree;   // null when no tree supplied

  PhyloMetric(const NumericMatrix& phylo_dist_, NumericMatrix& out_,
              const std::vector<std::string>& names,
              std::shared_ptr<const spacc::PhyloTree> tree_)
    : phylo_dist(phylo_dist_), out(out_),
      n_metrics(static_cast<int>(names.size())),
      mpd_idx(index_of(names, "mpd")), mntd_idx(index_of(names, "mntd")),
      pd_idx(tree_ ? index_of(names, "pd") : -1), rao_idx(index_of(names, "rao")),
      tree(tree_) {}

  struct Scratch {};
  Scratch scratch(int) const { return Scratch(); }

  void at_step(std::size_t s, int step, int, const std::vector<double>& cumulative,
               Scratch&) {
    auto d = [&](int i, int j) { return phylo_dist(i, j); };
    auto one = [](int) { return 1.0; };
    const std::size_t row = s * n_metrics;
    std::vector<int> present;
    if (mpd_idx >= 0 || mntd_idx >= 0) present = spacc::present_species(cumulative);
    if (mpd_idx >= 0)  out(row + mpd_idx, step) = spacc::mean_pairwise_distance(d, present, one);
    if (mntd_idx >= 0) out(row + mntd_idx, step) = spacc::mean_nearest_taxon_distance(d, present, one);
    if (pd_idx >= 0)   out(row + pd_idx, step) = tree->faith_pd([&](int tip) { return cumulative[tip] > 0; });
    if (rao_idx >= 0)  out(row + rao_idx, step) = spacc::rao_quadratic_entropy(
                         d, [&](int i) { return cumulative[i]; }, static_cast<int>(cumulative.size()));
  }
};

// [[Rcpp::export]]
List cpp_phylo_order_parallel(NumericMatrix species_pa, IntegerMatrix orders,
                              NumericMatrix phylo_dist_mat, CharacterVector metrics,
                              int n_cores = 1,
                              Rcpp::Nullable<IntegerMatrix> tree_edge = R_NilValue,
                              Rcpp::Nullable<NumericVector> tree_edge_length = R_NilValue,
                              int tree_n_tips = 0) {
  std::shared_ptr<const spacc::PhyloTree> tree;
  if (tree_edge.isNotNull() && tree_edge_length.isNotNull()) {
    tree = std::make_shared<const spacc::PhyloTree>(
      build_tree(Rcpp::as<IntegerMatrix>(tree_edge),
                 Rcpp::as<NumericVector>(tree_edge_length), tree_n_tips));
  }
  const int n_seeds = orders.nrow();
  NumericMatrix flat(n_seeds * metrics.size(), species_pa.nrow());
  std::fill(flat.begin(), flat.end(), NA_REAL);
  PhyloMetric metric(phylo_dist_mat, flat, metric_names(metrics), tree);
  run_order_metric<double>(species_pa, orders, metric, n_cores);
  return run_metric_list(flat, n_seeds, metrics);
}


// ===========================================================================
// FUNCTIONAL DIVERSITY
// ===========================================================================

// [[Rcpp::export]]
double calc_fdis(NumericMatrix traits, LogicalVector species_present,
                 NumericVector abundances) {
  return spacc::functional_dispersion(
    [&](int i, int t) { return traits(i, t); }, traits.ncol(),
    present_from_flags(species_present),
    [&](int i) { return abundances[i]; });
}


// [[Rcpp::export]]
double calc_fric_approx(NumericMatrix traits, LogicalVector species_present) {
  return spacc::functional_range_volume(
    [&](int i, int t) { return traits(i, t); }, traits.ncol(),
    present_from_flags(species_present));
}


// Rao's quadratic entropy on Euclidean trait distances.
// [[Rcpp::export]]
double calc_rao_traits(NumericMatrix traits, NumericVector abundances) {
  return spacc::rao_trait_entropy(
    [&](int i, int t) { return traits(i, t); }, traits.ncol(), traits.nrow(),
    [&](int i) { return abundances[i]; });
}


struct FuncMetric {
  RMatrix<double> traits;
  RMatrix<double> out;       // (n_seeds * n_metrics) x n_sites
  int n_metrics, fdis_idx, fric_idx, rao_idx;

  FuncMetric(const NumericMatrix& traits_, NumericMatrix& out_,
             const std::vector<std::string>& names)
    : traits(traits_), out(out_), n_metrics(static_cast<int>(names.size())),
      fdis_idx(index_of(names, "fdis")), fric_idx(index_of(names, "fric")),
      rao_idx(index_of(names, "rao")) {}

  struct Scratch {};
  Scratch scratch(int) const { return Scratch(); }

  void at_step(std::size_t s, int step, int, const std::vector<double>& cumulative,
               Scratch&) {
    auto trait = [&](int i, int t) { return traits(i, t); };
    auto abundance = [&](int i) { return cumulative[i]; };
    const int n_traits = traits.ncol();
    const std::size_t row = s * n_metrics;
    std::vector<int> present;
    if (fdis_idx >= 0 || fric_idx >= 0) present = spacc::present_species(cumulative);
    if (fdis_idx >= 0) out(row + fdis_idx, step) = spacc::functional_dispersion(trait, n_traits, present, abundance);
    if (fric_idx >= 0) out(row + fric_idx, step) = spacc::functional_range_volume(trait, n_traits, present);
    if (rao_idx >= 0)  out(row + rao_idx, step) = spacc::rao_trait_entropy(
                         trait, n_traits, static_cast<int>(cumulative.size()), abundance);
  }
};

// [[Rcpp::export]]
List cpp_func_order_parallel(NumericMatrix species_mat, IntegerMatrix orders,
                             NumericMatrix traits, CharacterVector metrics,
                             int n_cores = 1) {
  const int n_seeds = orders.nrow();
  NumericMatrix flat(n_seeds * metrics.size(), species_mat.nrow());
  std::fill(flat.begin(), flat.end(), NA_REAL);
  FuncMetric metric(traits, flat, metric_names(metrics));
  run_order_metric<double>(species_mat, orders, metric, n_cores);
  return run_metric_list(flat, n_seeds, metrics);
}
