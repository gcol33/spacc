#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Earth radius in km for haversine
const double EARTH_RADIUS_KM = 6371.0;

// Convert degrees to radians
inline double to_radians(double degrees) {
 return degrees * M_PI / 180.0;
}

// Haversine distance between two lat/lon points (in degrees)
inline double haversine_dist(double lat1, double lon1, double lat2, double lon2) {
 double dlat = to_radians(lat2 - lat1);
 double dlon = to_radians(lon2 - lon1);

 lat1 = to_radians(lat1);
 lat2 = to_radians(lat2);

 double a = std::sin(dlat / 2) * std::sin(dlat / 2) +
            std::cos(lat1) * std::cos(lat2) *
            std::sin(dlon / 2) * std::sin(dlon / 2);

 double c = 2 * std::atan2(std::sqrt(a), std::sqrt(1 - a));

 return EARTH_RADIUS_KM * c;
}

// Euclidean distance
inline double euclidean_dist(double x1, double y1, double x2, double y2) {
 double dx = x1 - x2;
 double dy = y1 - y2;
 return std::sqrt(dx * dx + dy * dy);
}


// [[Rcpp::export]]
NumericMatrix cpp_distance_matrix(NumericVector x, NumericVector y, std::string method = "euclidean") {
 int n = x.size();
 NumericMatrix dist(n, n);

 bool use_haversine = (method == "haversine");

 for (int i = 0; i < n; i++) {
   for (int j = i + 1; j < n; j++) {
     double d;
     if (use_haversine) {
       // For haversine: x = longitude, y = latitude
       d = haversine_dist(y[i], x[i], y[j], x[j]);
     } else {
       d = euclidean_dist(x[i], y[i], x[j], y[j]);
     }
     dist(i, j) = d;
     dist(j, i) = d;
   }
 }

 return dist;
}
