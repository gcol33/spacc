#ifndef SPACC_KDTREE_ADAPTER_H
#define SPACC_KDTREE_ADAPTER_H

#include <nanoflann.hpp>
#include <algorithm>
#include <memory>
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

typedef uint32_t IndexType;

typedef nanoflann::KDTreeSingleIndexAdaptor<
  nanoflann::L2_Simple_Adaptor<double, PointCloud2D>,
  PointCloud2D, 2, IndexType
> KDTree2D;

// k-d tree over planar coordinates with the same nearest_unvisited() interface
// as BallTree, so the order generators are written once against either index.
class KdTree {
public:
  // The cloud is heap-allocated because the nanoflann index keeps a reference
  // to it; the KdTree object itself can then be moved safely.
  KdTree(const std::vector<double>& x, const std::vector<double>& y)
    : n_sites(static_cast<int>(x.size())), cloud(new PointCloud2D()) {
    cloud->pts_x = x;
    cloud->pts_y = y;
    tree.reset(new KDTree2D(2, *cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10)));
    tree->buildIndex();
  }

  // Nearest unvisited site to (qx, qy): expanding k-NN queries until an
  // unvisited site appears; -1 when all are visited.
  int nearest_unvisited(double qx, double qy,
                        const std::vector<bool>& visited) const {
    double query[2] = {qx, qy};
    size_t k = 8;
    while (true) {
      if (k > static_cast<size_t>(n_sites)) k = static_cast<size_t>(n_sites);
      std::vector<IndexType> indices(k);
      std::vector<double> dists_sq(k);
      size_t found = tree->knnSearch(query, k, indices.data(), dists_sq.data());
      for (size_t i = 0; i < found; i++) {
        if (!visited[indices[i]]) return static_cast<int>(indices[i]);
      }
      if (k >= static_cast<size_t>(n_sites)) return -1;
      k = std::min(k * 2, static_cast<size_t>(n_sites));
    }
  }

private:
  int n_sites;
  std::unique_ptr<PointCloud2D> cloud;
  std::unique_ptr<KDTree2D> tree;
};

#endif // SPACC_KDTREE_ADAPTER_H
