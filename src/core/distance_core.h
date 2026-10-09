// core/distance_core.h
// Pure C++ distance calculations for spacc
// No Rcpp dependencies - testable with Catch2

#ifndef SPACC_CORE_DISTANCE_H
#define SPACC_CORE_DISTANCE_H

#include "types.h"
#include <cmath>
#include <string>
#include <limits>

namespace spacc {

inline double to_radians(double degrees) {
  return degrees * PI / 180.0;
}

inline double euclidean_distance(double x1, double y1, double x2, double y2) {
  double dx = x1 - x2;
  double dy = y1 - y2;
  return std::sqrt(dx * dx + dy * dy);
}

// Great-circle distance in kilometres between two (lat, lon) points in degrees.
inline double haversine_distance(double lat1, double lon1, double lat2, double lon2) {
  double dlat = to_radians(lat2 - lat1);
  double dlon = to_radians(lon2 - lon1);

  double rlat1 = to_radians(lat1);
  double rlat2 = to_radians(lat2);

  double a = std::sin(dlat / 2) * std::sin(dlat / 2) +
             std::cos(rlat1) * std::cos(rlat2) *
             std::sin(dlon / 2) * std::sin(dlon / 2);

  double c = 2 * std::atan2(std::sqrt(a), std::sqrt(1 - a));

  return EARTH_RADIUS_KM * c;
}

// Point distance in the metric selected by `haversine`. Coordinates are
// (x = lon, y = lat) in degrees for haversine, planar units otherwise.
struct PointDistance {
  bool haversine;
  double operator()(double x1, double y1, double x2, double y2) const {
    return haversine ? haversine_distance(y1, x1, y2, x2)
                     : euclidean_distance(x1, y1, x2, y2);
  }
};

// Full symmetric distance matrix. For haversine: x = longitude, y = latitude.
inline DistanceMatrix compute_distance_matrix(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::string& method = "euclidean") {

  size_t n = x.size();
  DistanceMatrix dist(n, std::vector<double>(n, 0.0));

  const PointDistance point_dist{method == "haversine"};

  for (size_t i = 0; i < n; i++) {
    for (size_t j = i + 1; j < n; j++) {
      double d = point_dist(x[i], y[i], x[j], y[j]);
      dist[i][j] = d;
      dist[j][i] = d;
    }
  }

  return dist;
}

// Nearest unvisited site to `current` under a precomputed distance matrix.
inline int find_nearest_unvisited(
    const DistanceMatrix& dist_mat,
    int current,
    const VisitedFlags& visited) {

  double min_dist = std::numeric_limits<double>::infinity();
  int nearest = -1;
  size_t n = visited.size();

  for (size_t j = 0; j < n; j++) {
    if (!visited[j] && dist_mat[current][j] < min_dist) {
      min_dist = dist_mat[current][j];
      nearest = static_cast<int>(j);
    }
  }

  return nearest;
}

} // namespace spacc

#endif // SPACC_CORE_DISTANCE_H
