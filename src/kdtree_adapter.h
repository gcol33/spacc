#ifndef SPACC_KDTREE_ADAPTER_H
#define SPACC_KDTREE_ADAPTER_H

#include "nanoflann.hpp"
#include <vector>
#include <cstddef>
#include <cstdint>

// Point cloud adapter for nanoflann (2D Euclidean)
struct PointCloud2D {
  std::vector<double> pts_x;
  std::vector<double> pts_y;

  inline size_t kdtree_get_point_count() const { return pts_x.size(); }

  inline double kdtree_get_pt(const size_t idx, const size_t dim) const {
    return (dim == 0) ? pts_x[idx] : pts_y[idx];
  }

  template <class BBOX>
  bool kdtree_get_bbox(BBOX&) const { return false; }
};

// Use uint32_t as index type to match nanoflann default
typedef uint32_t IndexType;

// Convenience typedefs
typedef nanoflann::KDTreeSingleIndexAdaptor<
  nanoflann::L2_Simple_Adaptor<double, PointCloud2D>,
  PointCloud2D, 2, IndexType
> KDTree2D;

// Build a k-d tree from coordinate vectors
// Caller owns the returned pointer (use delete)
inline KDTree2D* build_kdtree(const PointCloud2D& cloud) {
  auto* tree = new KDTree2D(2, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  tree->buildIndex();
  return tree;
}

// Find the nearest unvisited site to a query point.
// Searches expanding k-NN queries until an unvisited site is found.
// Returns the index into the original point cloud.
inline int find_nearest_unvisited(KDTree2D& tree, double qx, double qy,
                                   const std::vector<bool>& visited,
                                   int n_sites) {
  double query[2] = {qx, qy};

  // Start with k=8, expand if all results are visited
  size_t k = 8;

  while (true) {
    if (k > static_cast<size_t>(n_sites)) k = static_cast<size_t>(n_sites);

    std::vector<IndexType> indices(k);
    std::vector<double> dists_sq(k);

    size_t found = tree.knnSearch(query, k, indices.data(), dists_sq.data());

    for (size_t i = 0; i < found; i++) {
      if (!visited[indices[i]]) {
        return static_cast<int>(indices[i]);
      }
    }

    // All k results were visited; double k and retry
    if (k >= static_cast<size_t>(n_sites)) {
      return -1;
    }
    k = std::min(k * 2, static_cast<size_t>(n_sites));
  }
}

#endif // SPACC_KDTREE_ADAPTER_H
