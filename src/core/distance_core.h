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

// Convert degrees to radians
inline double to_radians(double degrees) {
  return degrees * PI / 180.0;
}

// Euclidean distance between two points
inline double euclidean_distance(double x1, double y1, double x2, double y2) {
  double dx = x1 - x2;
  double dy = y1 - y2;
  return std::sqrt(dx * dx + dy * dy);
}

// Haversine distance between two lat/lon points (in degrees)
// Returns distance in kilometers
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

// Compute full distance matrix
// For haversine: x = longitude, y = latitude
inline DistanceMatrix compute_distance_matrix(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::string& method = "euclidean") {

  size_t n = x.size();
  DistanceMatrix dist(n, std::vector<double>(n, 0.0));

  bool use_haversine = (method == "haversine");

  for (size_t i = 0; i < n; i++) {
    for (size_t j = i + 1; j < n; j++) {
      double d;
      if (use_haversine) {
        // For haversine: x = longitude, y = latitude
        d = haversine_distance(y[i], x[i], y[j], x[j]);
      } else {
        d = euclidean_distance(x[i], y[i], x[j], y[j]);
      }
      dist[i][j] = d;
      dist[j][i] = d;
    }
  }

  return dist;
}

// Find nearest unvisited site
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
