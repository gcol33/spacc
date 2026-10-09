#include <Rcpp.h>
#include <string>
#include "core/distance_core.h"
using namespace Rcpp;

// Full symmetric site distance matrix. For haversine: x = longitude,
// y = latitude in degrees, result in kilometres.
// [[Rcpp::export]]
NumericMatrix cpp_distance_matrix(NumericVector x, NumericVector y,
                                  std::string method = "euclidean") {
  const int n = x.size();
  NumericMatrix dist(n, n);
  const spacc::PointDistance point_dist{method == "haversine"};
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      double d = point_dist(x[i], y[i], x[j], y[j]);
      dist(i, j) = d;
      dist(j, i) = d;
    }
  }
  return dist;
}
