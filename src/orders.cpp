#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include "core/distance_core.h"
using namespace Rcpp;

// Site-visitation orders for custom-metric accumulation (spaccDiversity()).
// These return the 0-based order in which sites are visited, one row per seed,
// so an arbitrary R diversity function can be applied along the same spatial
// traversal used by the built-in accumulation methods.


// Nearest-neighbour walk: from each seed, repeatedly move to the nearest
// unvisited site using a precomputed distance matrix.
// [[Rcpp::export]]
IntegerMatrix cpp_nn_walk_order(NumericMatrix dist_mat, IntegerVector seeds) {
  int n_sites = dist_mat.nrow();
  int n_seeds = seeds.size();
  IntegerMatrix orders(n_seeds, n_sites);

  for (int s = 0; s < n_seeds; s++) {
    std::vector<bool> visited(n_sites, false);
    int current = seeds[s];
    visited[current] = true;
    orders(s, 0) = current;

    for (int step = 1; step < n_sites; step++) {
      double min_dist = R_PosInf;
      int nearest = -1;
      for (int j = 0; j < n_sites; j++) {
        if (!visited[j] && dist_mat(current, j) < min_dist) {
          min_dist = dist_mat(current, j);
          nearest = j;
        }
      }
      current = nearest;
      visited[current] = true;
      orders(s, step) = current;
    }
  }

  return orders;
}


// Fixed-focus spatial ordering. Each row ranks every sampling site by its
// distance from one continuous focal point.
// [[Rcpp::export]]
IntegerMatrix cpp_focal_order(NumericVector x, NumericVector y,
                              NumericVector focus_x, NumericVector focus_y,
                              std::string distance = "euclidean") {
  int n_sites = x.size();
  int n_foci = focus_x.size();
  IntegerMatrix orders(n_foci, n_sites);
  bool use_haversine = distance == "haversine";

  for (int focus = 0; focus < n_foci; focus++) {
    std::vector<std::pair<double, int>> ranked(n_sites);
    for (int site = 0; site < n_sites; site++) {
      double separation = use_haversine
        ? spacc::haversine_distance(focus_y[focus], focus_x[focus], y[site], x[site])
        : spacc::euclidean_distance(focus_x[focus], focus_y[focus], x[site], y[site]);
      ranked[site] = std::make_pair(separation, site);
    }
    std::stable_sort(ranked.begin(), ranked.end(),
      [](const std::pair<double, int>& left, const std::pair<double, int>& right) {
        if (left.first == right.first) return left.second < right.second;
        return left.first < right.first;
      });
    for (int step = 0; step < n_sites; step++) {
      orders(focus, step) = ranked[step].second;
    }
  }
  return orders;
}


// kNCN order: from each seed, repeatedly move to the site nearest the centroid
// of the already-visited sites. Mirrors the traversal in cpp_kncn_*.
// [[Rcpp::export]]
IntegerMatrix cpp_kncn_order(NumericVector x, NumericVector y, IntegerVector seeds) {
  int n_sites = x.size();
  int n_seeds = seeds.size();
  IntegerMatrix orders(n_seeds, n_sites);

  for (int s = 0; s < n_seeds; s++) {
    std::vector<bool> visited(n_sites, false);
    int current = seeds[s];
    visited[current] = true;
    orders(s, 0) = current;

    double cx = x[current];
    double cy = y[current];
    int n_visited = 1;

    for (int step = 1; step < n_sites; step++) {
      double min_dist = R_PosInf;
      int nearest = -1;
      for (int j = 0; j < n_sites; j++) {
        if (!visited[j]) {
          double dx = cx - x[j];
          double dy = cy - y[j];
          double d = dx * dx + dy * dy;
          if (d < min_dist) {
            min_dist = d;
            nearest = j;
          }
        }
      }
      current = nearest;
      visited[current] = true;
      cx = (cx * n_visited + x[current]) / (n_visited + 1);
      cy = (cy * n_visited + y[current]) / (n_visited + 1);
      n_visited++;
      orders(s, step) = current;
    }
  }

  return orders;
}
