#ifndef SPACC_BALLTREE_H
#define SPACC_BALLTREE_H

#include <vector>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>

// ============================================================================
// Ball tree supporting arbitrary distance metrics
// Designed for haversine nearest-neighbor queries on geographic coordinates.
//
// A ball tree recursively partitions points into nested bounding spheres.
// Each node stores a center (pivot) and a radius that encloses all points
// in the subtree. NN search prunes subtrees whose closest possible point
// (center distance minus radius) exceeds the current best distance.
// ============================================================================

static constexpr double BALLTREE_EARTH_R = 6371.0;
static constexpr double BALLTREE_D2R = 3.14159265358979323846 / 180.0;

// Haversine distance in km between two (lon, lat) points in degrees
inline double haversine(double lon1, double lat1, double lon2, double lat2) {
  double dlat = (lat2 - lat1) * BALLTREE_D2R;
  double dlon = (lon2 - lon1) * BALLTREE_D2R;
  double la1 = lat1 * BALLTREE_D2R;
  double la2 = lat2 * BALLTREE_D2R;
  double a = std::sin(dlat * 0.5) * std::sin(dlat * 0.5) +
             std::cos(la1) * std::cos(la2) *
             std::sin(dlon * 0.5) * std::sin(dlon * 0.5);
  double c = 2.0 * std::atan2(std::sqrt(a), std::sqrt(1.0 - a));
  return BALLTREE_EARTH_R * c;
}


struct BallTreeNode {
  int pivot;          // Index of the pivot point
  double radius;      // Bounding radius (max distance from pivot to any point in subtree)
  int left;           // Index of left child node (-1 if leaf)
  int right;          // Index of right child node (-1 if leaf)
  int start;          // Start index in the index array (for leaves)
  int end;            // End index in the index array (for leaves)
};


class BallTree {
public:
  // Coordinate arrays (x = longitude, y = latitude for haversine)
  const std::vector<double>& px;
  const std::vector<double>& py;
  int n_points;

  // Tree structure
  std::vector<BallTreeNode> nodes;
  std::vector<int> indices;  // Permuted point indices

  static constexpr int LEAF_SIZE = 16;

  BallTree(const std::vector<double>& x, const std::vector<double>& y)
    : px(x), py(y), n_points(static_cast<int>(x.size()))
  {
    indices.resize(n_points);
    for (int i = 0; i < n_points; i++) indices[i] = i;
    nodes.reserve(2 * n_points / LEAF_SIZE + 1);
    build(0, n_points);
  }

  // Find the nearest unvisited point to a query (qx, qy).
  // Returns the index in the original point array, or -1 if all visited.
  int find_nearest_unvisited(double qx, double qy,
                             const std::vector<bool>& visited) const {
    double best_dist = std::numeric_limits<double>::infinity();
    int best_idx = -1;
    search_nn(0, qx, qy, visited, best_dist, best_idx);
    return best_idx;
  }

private:
  // Distance between two points by index
  double dist(int i, int j) const {
    return haversine(px[i], py[i], px[j], py[j]);
  }

  // Distance from query to point by index
  double dist_to(double qx, double qy, int j) const {
    return haversine(qx, qy, px[j], py[j]);
  }

  // Build subtree over indices[start..end)
  int build(int start, int end) {
    int node_idx = static_cast<int>(nodes.size());
    nodes.push_back(BallTreeNode());

    int count = end - start;

    if (count <= LEAF_SIZE) {
      // Leaf: pick medoid (point minimizing max distance to others) as pivot
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

    // Find two far-apart points to define the split axis.
    // Precompute distances to avoid redundant haversine calls.
    int p1 = indices[start];
    std::vector<double> d_from_p1(count);
    int far1_local = 0;
    for (int i = 0; i < count; i++) {
      d_from_p1[i] = dist(p1, indices[start + i]);
      if (d_from_p1[i] > d_from_p1[far1_local]) far1_local = i;
    }
    int far1 = indices[start + far1_local];

    // Distances from far1 — reused for partitioning
    std::vector<double> d_far1(count);
    int far2_local = 0;
    for (int i = 0; i < count; i++) {
      d_far1[i] = dist(far1, indices[start + i]);
      if (d_far1[i] > d_far1[far2_local]) far2_local = i;
    }
    int far2 = indices[start + far2_local];

    // Distances from far2 — precompute once
    std::vector<double> d_far2(count);
    for (int i = 0; i < count; i++) {
      d_far2[i] = dist(far2, indices[start + i]);
    }

    // Projection values for partitioning (no haversine calls in comparator)
    std::vector<double> proj(count);
    for (int i = 0; i < count; i++) {
      proj[i] = d_far1[i] - d_far2[i];
    }

    // Partition using nth_element on precomputed projections.
    // We need indices and proj to stay synchronized, so sort a local
    // index array and then reorder indices[start..end).
    std::vector<int> order(count);
    std::iota(order.begin(), order.end(), 0);
    int mid_local = count / 2;
    std::nth_element(order.begin(), order.begin() + mid_local, order.end(),
      [&](int a, int b) { return proj[a] < proj[b]; }
    );

    // Reorder indices according to the partition
    std::vector<int> tmp(count);
    for (int i = 0; i < count; i++) tmp[i] = indices[start + order[i]];
    for (int i = 0; i < count; i++) indices[start + i] = tmp[i];

    int mid = start + mid_local;

    // Pick medoid of the full set as pivot (point minimizing sum of distances).
    // Reorder d_far1 to match new index order and use it as a proxy:
    // the point closest to the midpoint of the far1-far2 axis (proj ~ 0)
    // is a good centroid approximation.
    // Reorder proj to match
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

    // Compute radius from pivot (single pass of haversine calls)
    double radius = 0.0;
    for (int i = start; i < end; i++) {
      double d = dist(pivot_pt, indices[i]);
      if (d > radius) radius = d;
    }

    nodes[node_idx].pivot = pivot_pt;
    nodes[node_idx].radius = radius;
    nodes[node_idx].start = start;
    nodes[node_idx].end = end;

    // Recurse — must index into nodes[] since vector may reallocate
    int left_idx = build(start, mid);
    nodes[node_idx].left = left_idx;
    int right_idx = build(mid, end);
    nodes[node_idx].right = right_idx;

    return node_idx;
  }

  // Recursive NN search with pruning
  void search_nn(int node_idx, double qx, double qy,
                 const std::vector<bool>& visited,
                 double& best_dist, int& best_idx) const {
    const BallTreeNode& node = nodes[node_idx];

    // Prune: if the closest possible point in this ball is farther than best
    double d_pivot = dist_to(qx, qy, node.pivot);
    double lower_bound = d_pivot - node.radius;
    if (lower_bound >= best_dist) return;

    // Check pivot itself (it's a real point, avoid recomputing its distance)
    if (!visited[node.pivot] && d_pivot < best_dist) {
      best_dist = d_pivot;
      best_idx = node.pivot;
    }

    // Leaf node: check all points
    if (node.left == -1) {
      for (int i = node.start; i < node.end; i++) {
        int idx = indices[i];
        if (idx == node.pivot) continue;  // Already checked above
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

    // Internal node: visit closer child first for better pruning.
    // Reuse d_pivot for the child whose pivot matches this node's pivot
    // to avoid a redundant haversine call.
    int first = node.left;
    int second = node.right;
    const BallTreeNode& lnode = nodes[node.left];
    const BallTreeNode& rnode = nodes[node.right];

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
