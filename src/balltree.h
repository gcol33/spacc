#ifndef SPACC_BALLTREE_H
#define SPACC_BALLTREE_H

#include <vector>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include "core/distance_core.h"

// ============================================================================
// Ball tree for haversine nearest-neighbour queries on (lon, lat) coordinates.
//
// Points are recursively partitioned into nested bounding spheres. Each node
// stores a pivot point and a radius enclosing every point in its subtree; the
// search prunes subtrees whose closest possible point (pivot distance minus
// radius) is farther than the current best.
// ============================================================================

class BallTree {
public:
  // Coordinate arrays: x = longitude, y = latitude (degrees)
  const std::vector<double>& px;
  const std::vector<double>& py;
  int n_points;

  struct Node {
    int pivot;          // index of the pivot point
    double radius;      // max distance from pivot to any point in the subtree
    int left;           // left child node (-1 if leaf)
    int right;          // right child node (-1 if leaf)
    int start;          // start index in `indices` (leaves)
    int end;            // end index in `indices` (leaves)
  };

  std::vector<Node> nodes;
  std::vector<int> indices;  // permuted point indices

  static constexpr int LEAF_SIZE = 16;

  BallTree(const std::vector<double>& x, const std::vector<double>& y)
    : px(x), py(y), n_points(static_cast<int>(x.size()))
  {
    indices.resize(n_points);
    for (int i = 0; i < n_points; i++) indices[i] = i;
    nodes.reserve(2 * n_points / LEAF_SIZE + 1);
    build(0, n_points);
  }

  // Nearest unvisited point to the query (qx, qy); -1 when all are visited.
  int nearest_unvisited(double qx, double qy,
                        const std::vector<bool>& visited) const {
    double best_dist = std::numeric_limits<double>::infinity();
    int best_idx = -1;
    search_nn(0, qx, qy, visited, best_dist, best_idx);
    return best_idx;
  }

private:
  double dist(int i, int j) const {
    return spacc::haversine_distance(py[i], px[i], py[j], px[j]);
  }

  double dist_to(double qx, double qy, int j) const {
    return spacc::haversine_distance(qy, qx, py[j], px[j]);
  }

  // Build the subtree over indices[start..end)
  int build(int start, int end) {
    int node_idx = static_cast<int>(nodes.size());
    nodes.push_back(Node());

    int count = end - start;

    if (count <= LEAF_SIZE) {
      // Leaf: the medoid (point minimising its max distance to the others) is the pivot
      int best_pivot = indices[start];
      double best_max_d = std::numeric_limits<double>::infinity();
      for (int i = start; i < end; i++) {
        double max_d = 0.0;
        for (int j = start; j < end; j++) {
          double d = dist(indices[i], indices[j]);
          if (d > max_d) max_d = d;
        }
        if (max_d < best_max_d) {
          best_max_d = max_d;
          best_pivot = indices[i];
        }
      }
      nodes[node_idx].pivot = best_pivot;
      nodes[node_idx].radius = best_max_d;
      nodes[node_idx].left = -1;
      nodes[node_idx].right = -1;
      nodes[node_idx].start = start;
      nodes[node_idx].end = end;
      return node_idx;
    }

    // Two far-apart points define the split axis
    int p1 = indices[start];
    std::vector<double> d_from_p1(count);
    int far1_local = 0;
    for (int i = 0; i < count; i++) {
      d_from_p1[i] = dist(p1, indices[start + i]);
      if (d_from_p1[i] > d_from_p1[far1_local]) far1_local = i;
    }
    int far1 = indices[start + far1_local];

    std::vector<double> d_far1(count);
    int far2_local = 0;
    for (int i = 0; i < count; i++) {
      d_far1[i] = dist(far1, indices[start + i]);
      if (d_far1[i] > d_far1[far2_local]) far2_local = i;
    }
    int far2 = indices[start + far2_local];

    std::vector<double> d_far2(count);
    for (int i = 0; i < count; i++) {
      d_far2[i] = dist(far2, indices[start + i]);
    }

    // Projection on the far1-far2 axis drives the partition
    std::vector<double> proj(count);
    for (int i = 0; i < count; i++) proj[i] = d_far1[i] - d_far2[i];

    std::vector<int> order(count);
    std::iota(order.begin(), order.end(), 0);
    int mid_local = count / 2;
    std::nth_element(order.begin(), order.begin() + mid_local, order.end(),
      [&](int a, int b) { return proj[a] < proj[b]; }
    );

    std::vector<int> tmp(count);
    for (int i = 0; i < count; i++) tmp[i] = indices[start + order[i]];
    for (int i = 0; i < count; i++) indices[start + i] = tmp[i];

    int mid = start + mid_local;

    // Pivot: the point closest to the midpoint of the axis (proj ~ 0)
    std::vector<double> proj_reordered(count);
    for (int i = 0; i < count; i++) proj_reordered[i] = proj[order[i]];

    int best_pivot_local = 0;
    double best_abs_proj = std::abs(proj_reordered[0]);
    for (int i = 1; i < count; i++) {
      double ap = std::abs(proj_reordered[i]);
      if (ap < best_abs_proj) {
        best_abs_proj = ap;
        best_pivot_local = i;
      }
    }
    int pivot_pt = indices[start + best_pivot_local];

    double radius = 0.0;
    for (int i = start; i < end; i++) {
      double d = dist(pivot_pt, indices[i]);
      if (d > radius) radius = d;
    }

    nodes[node_idx].pivot = pivot_pt;
    nodes[node_idx].radius = radius;
    nodes[node_idx].start = start;
    nodes[node_idx].end = end;

    // Index into nodes[] after each recursion: the vector may reallocate
    int left_idx = build(start, mid);
    nodes[node_idx].left = left_idx;
    int right_idx = build(mid, end);
    nodes[node_idx].right = right_idx;

    return node_idx;
  }

  void search_nn(int node_idx, double qx, double qy,
                 const std::vector<bool>& visited,
                 double& best_dist, int& best_idx) const {
    const Node& node = nodes[node_idx];

    double d_pivot = dist_to(qx, qy, node.pivot);
    double lower_bound = d_pivot - node.radius;
    if (lower_bound >= best_dist) return;

    if (!visited[node.pivot] && d_pivot < best_dist) {
      best_dist = d_pivot;
      best_idx = node.pivot;
    }

    if (node.left == -1) {
      for (int i = node.start; i < node.end; i++) {
        int idx = indices[i];
        if (idx == node.pivot) continue;
        if (!visited[idx]) {
          double d = dist_to(qx, qy, idx);
          if (d < best_dist) {
            best_dist = d;
            best_idx = idx;
          }
        }
      }
      return;
    }

    // Visit the closer child first for better pruning; a child sharing this
    // node's pivot reuses d_pivot.
    int first = node.left;
    int second = node.right;
    const Node& lnode = nodes[node.left];
    const Node& rnode = nodes[node.right];

    double d_left = (lnode.pivot == node.pivot) ? d_pivot : dist_to(qx, qy, lnode.pivot);
    double d_right = (rnode.pivot == node.pivot) ? d_pivot : dist_to(qx, qy, rnode.pivot);

    if (d_right < d_left) {
      first = node.right;
      second = node.left;
    }

    search_nn(first, qx, qy, visited, best_dist, best_idx);
    search_nn(second, qx, qy, visited, best_dist, best_idx);
  }
};

#endif // SPACC_BALLTREE_H
