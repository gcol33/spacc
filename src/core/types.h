// core/types.h
// Pure C++ type definitions for spacc
// No Rcpp dependencies - testable with Catch2

#ifndef SPACC_CORE_TYPES_H
#define SPACC_CORE_TYPES_H

#include <vector>
#include <set>
#include <cstddef>

namespace spacc {

// Type aliases for clarity
using Curve = std::vector<int>;
using CurveDouble = std::vector<double>;
using CurveMatrix = std::vector<std::vector<int>>;
using CurveMatrixDouble = std::vector<std::vector<double>>;

using SiteSpeciesMatrix = std::vector<std::vector<int>>;
using SiteSpeciesMatrixDouble = std::vector<std::vector<double>>;
using DistanceMatrix = std::vector<std::vector<double>>;

using SpeciesSet = std::set<int>;
using VisitedFlags = std::vector<bool>;

// Constants
constexpr double EARTH_RADIUS_KM = 6371.0;
constexpr double PI = 3.14159265358979323846;

} // namespace spacc

#endif // SPACC_CORE_TYPES_H
