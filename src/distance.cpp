#include <Rcpp.h>
using namespace Rcpp;

//' Fast Euclidean Distance Matrix
//'
//' Compute pairwise Euclidean distances between points.
//' O(n²) but fully vectorized in C++.
//'
//' @param x Numeric vector of x coordinates
//' @param y Numeric vector of y coordinates
//' @return Symmetric n×n distance matrix
//'
// [[Rcpp::export]]
NumericMatrix cpp_distance_matrix(NumericVector x, NumericVector y) {
  int n = x.size();
  NumericMatrix dist(n, n);

  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      double dx = x[i] - x[j];
      double dy = y[i] - y[j];
      double d = std::sqrt(dx * dx + dy * dy);
      dist(i, j) = d;
      dist(j, i) = d;
    }
  }

  return dist;
}
