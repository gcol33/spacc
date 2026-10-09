// core/phylo_core.h
// Phylogenetic diversity of a species set: mean pairwise distance (MPD), mean
// nearest taxon distance (MNTD), Rao's quadratic entropy, Faith's PD and
// EvoHeritage. Pure C++: distance matrices and abundance vectors are accessed
// through d(i, j) / abundance(i) callables.

#ifndef SPACC_CORE_PHYLO_H
#define SPACC_CORE_PHYLO_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace spacc {

// Indices of the species whose cumulative abundance is positive.
template<class T>
std::vector<int> present_species(const std::vector<T>& cumulative) {
  std::vector<int> present;
  present.reserve(cumulative.size());
  for (size_t i = 0; i < cumulative.size(); i++) {
    if (cumulative[i] > 0) present.push_back(static_cast<int>(i));
  }
  return present;
}

// Mean pairwise distance among `present`; `weight(i)` gives an optional
// abundance weight per species (pairs weighted by w_i * w_j). Pass a callable
// returning 1.0 for the unweighted mean.
template<class DistFn, class WeightFn>
double mean_pairwise_distance(const DistFn& d, const std::vector<int>& present,
                              const WeightFn& weight) {
  if (present.size() < 2) return 0.0;
  double sum_dist = 0.0, sum_weight = 0.0;
  for (size_t i = 0; i < present.size(); i++) {
    for (size_t j = i + 1; j < present.size(); j++) {
      double w = weight(present[i]) * weight(present[j]);
      sum_dist += d(present[i], present[j]) * w;
      sum_weight += w;
    }
  }
  return sum_weight == 0.0 ? 0.0 : sum_dist / sum_weight;
}

// Mean distance from each present species to its nearest present neighbour,
// weighted by `weight(i)`.
template<class DistFn, class WeightFn>
double mean_nearest_taxon_distance(const DistFn& d, const std::vector<int>& present,
                                   const WeightFn& weight) {
  if (present.size() < 2) return 0.0;
  double sum_nnd = 0.0, sum_weight = 0.0;
  for (size_t i = 0; i < present.size(); i++) {
    double min_d = std::numeric_limits<double>::infinity();
    for (size_t j = 0; j < present.size(); j++) {
      if (i == j) continue;
      double dij = d(present[i], present[j]);
      if (dij < min_d) min_d = dij;
    }
    double w = weight(present[i]);
    sum_nnd += min_d * w;
    sum_weight += w;
  }
  return sum_weight == 0.0 ? 0.0 : sum_nnd / sum_weight;
}

// Rao's quadratic entropy Q = sum_i sum_j p_i p_j d_ij with p the relative
// abundances over all n species (absent species contribute zero).
template<class DistFn, class AbundFn>
double rao_quadratic_entropy(const DistFn& d, const AbundFn& abundance, int n_species) {
  double total = 0.0;
  for (int i = 0; i < n_species; i++) total += abundance(i);
  if (total <= 0.0) return 0.0;
  double q = 0.0;
  for (int i = 0; i < n_species; i++) {
    double ai = abundance(i);
    if (ai <= 0.0) continue;
    double pi = ai / total;
    for (int j = i + 1; j < n_species; j++) {
      double aj = abundance(j);
      if (aj <= 0.0) continue;
      q += pi * (aj / total) * d(i, j);
    }
  }
  return 2.0 * q;  // (i, j) and (j, i); the diagonal is zero
}

// Rooted tree in ape's edge-list encoding: nodes are numbered 1..n_nodes with
// tips 1..n_tips; edge e joins parent(e) to child(e).
struct PhyloTree {
  std::vector<int> edge_into;      // edge index leading into each node, -1 at the root
  std::vector<int> parent_of;      // parent node of each node, 0 at the root
  std::vector<double> edge_length;
  std::vector<int> postorder;      // non-root nodes, every child before its parent
  int n_tips;

  template<class ParentFn, class ChildFn, class LengthFn>
  PhyloTree(int n_edges, const ParentFn& parent, const ChildFn& child,
            const LengthFn& length, int n_tips_)
    : n_tips(n_tips_) {
    int n_nodes = 0;
    for (int e = 0; e < n_edges; e++) {
      n_nodes = std::max(n_nodes, std::max(parent(e), child(e)));
    }
    edge_into.assign(n_nodes + 1, -1);
    parent_of.assign(n_nodes + 1, 0);
    edge_length.resize(n_edges);
    for (int e = 0; e < n_edges; e++) {
      edge_into[child(e)] = e;
      parent_of[child(e)] = parent(e);
      edge_length[e] = length(e);
    }

    // Sorting nodes by decreasing depth (edge count to the root) puts every
    // child ahead of its parent.
    std::vector<int> depth(n_nodes + 1, -1);
    for (int v = 1; v <= n_nodes; v++) {
      std::vector<int> path;
      int u = v;
      while (depth[u] < 0 && edge_into[u] >= 0) {
        path.push_back(u);
        u = parent_of[u];
      }
      if (depth[u] < 0) depth[u] = 0;
      for (auto it = path.rbegin(); it != path.rend(); ++it) depth[*it] = depth[parent_of[*it]] + 1;
    }
    for (int v = 1; v <= n_nodes; v++) {
      if (edge_into[v] >= 0) postorder.push_back(v);
    }
    std::stable_sort(postorder.begin(), postorder.end(),
                     [&](int a, int b) { return depth[a] > depth[b]; });
  }

  // Faith's PD: total length of the edges on the paths from the present tips
  // to the root. present[tip] for tips 0..n_tips-1 (0-based).
  template<class PresentFn>
  double faith_pd(const PresentFn& present) const {
    std::vector<bool> used(edge_length.size(), false);
    double pd = 0.0;
    for (int tip = 0; tip < n_tips; tip++) {
      if (!present(tip)) continue;
      int node = tip + 1;
      while (true) {
        int e = edge_into[node];
        if (e < 0 || used[e]) break;  // root reached, or path already counted
        used[e] = true;
        pd += edge_length[e];
        node = parent_of[node];
      }
    }
    return pd;
  }

  // Phylogenetic Hill number of order q, the mean phylogenetic diversity of
  // Chao et al. (2010): with a_e the relative abundance descending from edge
  // e and T = sum_e L_e a_e the abundance-weighted mean root-to-tip depth (the
  // tree depth for an ultrametric tree),
  //   qD(T) = (sum_e (L_e / T) a_e^q)^(1/(1-q)),  exp(-sum_e (L_e/T) a_e log a_e) at q = 1,
  // in effective numbers of lineages; it equals the taxonomic Hill number on
  // a star tree. abundance(tip) for tips 0..n_tips-1.
  template<class AbundanceFn>
  double hill(const AbundanceFn& abundance, double q) const {
    std::vector<double> below(parent_of.size(), 0.0);
    double total = 0.0;
    for (int tip = 0; tip < n_tips; tip++) {
      const double a = abundance(tip);
      if (a > 0) { below[tip + 1] = a; total += a; }
    }
    if (total <= 0.0) return 0.0;
    for (int v : postorder) below[parent_of[v]] += below[v];

    double depth = 0.0;
    for (int v : postorder) depth += edge_length[edge_into[v]] * below[v] / total;
    if (depth <= 0.0) return 0.0;

    const bool shannon = std::abs(q - 1.0) < 1e-10;
    double sum = 0.0;
    for (int v : postorder) {
      const double a = below[v] / total;
      if (a <= 0.0) continue;
      const double w = edge_length[edge_into[v]] / depth;
      sum += shannon ? -w * a * std::log(a) : w * std::pow(a, q);
    }
    return shannon ? std::exp(sum) : std::pow(sum, 1.0 / (1.0 - q));
  }

  // EvoHeritage phi_rho of the present tips (Rosindell et al. 2024, Eqns. 5
  // and 9). Heritage is gained along each edge and attrited at rate rho per
  // unit length; p(v) is the fraction of the heritage at node v that survives
  // to at least one present tip:
  //   p(tip) = present(tip),  p(v) = 1 - prod_children [1 - p(c) exp(-rho L_c)].
  // phi_rho = sum_e (1 - exp(-rho L_e)) p(child_e) / (1 - exp(-rho)), which is
  // Faith's PD at rho = 0 and tends to the number of present tips as rho grows.
  template<class PresentFn>
  double evo_heritage(const PresentFn& present, double rho) const {
    // (1 - exp(-rho L)) / (1 - exp(-rho)) in standardised units; L at rho = 0.
    const double unit = rho > 0.0 ? -std::expm1(-rho) : 1.0;
    auto gain = [&](double len) { return rho > 0.0 ? -std::expm1(-rho * len) / unit : len; };

    std::vector<double> lost_below(parent_of.size(), 1.0);  // prod_children [1 - p(c) R(c)]
    double phi = 0.0;
    for (int v : postorder) {
      const double p = v <= n_tips ? (present(v - 1) ? 1.0 : 0.0) : 1.0 - lost_below[v];
      if (p <= 0.0) continue;
      const double len = edge_length[edge_into[v]];
      phi += gain(len) * p;
      lost_below[parent_of[v]] *= 1.0 - p * std::exp(-rho * len);
    }
    return phi;
  }
};

} // namespace spacc

#endif // SPACC_CORE_PHYLO_H
