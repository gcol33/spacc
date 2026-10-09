// test_accumulation.cpp
// Catch2 tests for the order generators and accumulation kernels

#include "catch.hpp"
#include "../core/accumulation_core.h"

using namespace spacc;

namespace {

// Accessors over std containers
struct Species {
  SiteSpeciesMatrix m;
  int operator()(int site, int sp) const { return m[site][sp]; }
};

struct Dist {
  DistanceMatrix d;
  double operator()(int i, int j) const { return d[i][j]; }
};

struct OrderOf {
  std::vector<int> v;
  int operator()(int step) const { return v[step]; }
};

std::vector<int> richness_curve(const SiteSpeciesMatrix& species, const std::vector<int>& order) {
  std::vector<int> curve(order.size());
  richness_along_order(Species{species}, static_cast<int>(species.size()),
                       static_cast<int>(species[0].size()), OrderOf{order},
                       [&](int step, int r) { curve[step] = r; });
  return curve;
}

std::vector<int> walk_by_matrix(const DistanceMatrix& d, int seed) {
  const int n = static_cast<int>(d.size());
  Dist dist{d};
  std::vector<int> order;
  nn_walk_order(n, seed,
    [&](int current, const VisitedFlags& visited) {
      return nearest_unvisited(dist, current, n, visited);
    }, order);
  return order;
}

bool is_permutation(const std::vector<int>& order, int n) {
  std::vector<int> sorted(order);
  std::sort(sorted.begin(), sorted.end());
  for (int i = 0; i < n; i++) if (sorted[i] != i) return false;
  return static_cast<int>(order.size()) == n;
}

} // namespace


TEST_CASE("richness_along_order", "[accumulation]") {
  SiteSpeciesMatrix species = {
    {1, 0, 0},  // site 0: species 0
    {1, 1, 0},  // site 1: species 0, 1
    {0, 1, 1}   // site 2: species 1, 2
  };

  SECTION("counts new species along the order") {
    std::vector<int> curve = richness_curve(species, {0, 1, 2});
    REQUIRE(curve == std::vector<int>{1, 2, 3});
  }

  SECTION("order changes the curve, not the endpoint") {
    std::vector<int> curve = richness_curve(species, {2, 0, 1});
    REQUIRE(curve == std::vector<int>{2, 3, 3});
  }

  SECTION("same species everywhere") {
    SiteSpeciesMatrix same = {{1, 1}, {1, 1}, {1, 1}};
    REQUIRE(richness_curve(same, {0, 1, 2}) == std::vector<int>{2, 2, 2});
  }

  SECTION("single site") {
    SiteSpeciesMatrix one = {{1, 1, 0, 1}};
    REQUIRE(richness_curve(one, {0}) == std::vector<int>{3});
  }
}


TEST_CASE("accumulate_along_order", "[accumulation]") {
  SiteSpeciesMatrix species = {{5, 0}, {1, 2}, {0, 3}};
  std::vector<int> cumulative(2, 0);
  std::vector<std::vector<int>> seen;
  accumulate_along_order(Species{species}, 3, 2, OrderOf{{2, 0, 1}}, cumulative,
    [&](int, int) { seen.push_back(cumulative); });

  REQUIRE(seen.size() == 3);
  REQUIRE(seen[0] == std::vector<int>{0, 3});
  REQUIRE(seen[1] == std::vector<int>{5, 3});
  REQUIRE(seen[2] == std::vector<int>{6, 5});
}


TEST_CASE("nn_walk_order", "[accumulation]") {
  // Sites on a line: 0 -- 1 -- 2 -- 3
  DistanceMatrix d = {
    {0.0, 1.0, 2.0, 3.0},
    {1.0, 0.0, 1.0, 2.0},
    {2.0, 1.0, 0.0, 1.0},
    {3.0, 2.0, 1.0, 0.0}
  };

  SECTION("walks away from the seed along the line") {
    REQUIRE(walk_by_matrix(d, 0) == std::vector<int>{0, 1, 2, 3});
    REQUIRE(walk_by_matrix(d, 3) == std::vector<int>{3, 2, 1, 0});
  }

  SECTION("from the middle, ties resolve to the lowest index then continue") {
    REQUIRE(walk_by_matrix(d, 1) == std::vector<int>{1, 0, 2, 3});
  }

  SECTION("visits every site once") {
    REQUIRE(is_permutation(walk_by_matrix(d, 2), 4));
  }
}


TEST_CASE("kncn_order", "[accumulation]") {
  // Sites at (0,0), (2,0), (1,2): from (0,0) the centroid is the seed, so
  // (2,0) at distance 2 beats (1,2) at sqrt(5).
  std::vector<double> x = {0, 2, 1};
  std::vector<double> y = {0, 0, 2};
  PointDistance euclid{false};
  std::vector<int> order;
  kncn_order(x, y, 0,
    [&](double cx, double cy, const VisitedFlags& visited) {
      return nearest_unvisited_point(cx, cy, x, y, visited, euclid);
    }, order);

  REQUIRE(order == std::vector<int>{0, 1, 2});
}


TEST_CASE("focal_order", "[accumulation]") {
  std::vector<double> x = {0.0, 3.0, 1.0, 1.0};
  std::vector<double> y = {0.0, 0.0, 0.0, 0.0};
  std::vector<int> order;

  SECTION("ranks sites by distance from the focus") {
    double fx = 2.9;
    focal_order(4, [&](int site) { return std::fabs(x[site] - fx); }, order);
    REQUIRE(order == std::vector<int>{1, 2, 3, 0});
  }

  SECTION("ties keep site index order") {
    double fx = 1.0;
    focal_order(4, [&](int site) { return std::fabs(x[site] - fx); }, order);
    REQUIRE(order == std::vector<int>{2, 3, 0, 1});
  }
}


TEST_CASE("gaussian_order", "[accumulation]") {
  DistanceMatrix d = {
    {0.0, 1.0, 10.0},
    {1.0, 0.0, 9.0},
    {10.0, 9.0, 0.0}
  };
  Dist dist{d};
  std::vector<int> order;

  SECTION("a draw of zero always takes the first unvisited candidate") {
    gaussian_order(dist, 3, 0, 1.0, []() { return 0.0; }, order);
    REQUIRE(order == std::vector<int>{0, 1, 2});
  }

  SECTION("a draw just below one takes the last candidate with positive weight") {
    gaussian_order(dist, 3, 0, 5.0, []() { return 0.999999; }, order);
    REQUIRE(order == std::vector<int>{0, 2, 1});
  }

  SECTION("underflowed weights fall back to the nearest unvisited site") {
    gaussian_order(dist, 3, 0, 1e-3, []() { return 0.5; }, order);
    REQUIRE(order == std::vector<int>{0, 1, 2});
  }

  SECTION("always a permutation") {
    double u = 0.37;
    gaussian_order(dist, 3, 2, 2.0, [&]() { u = std::fmod(u * 7.1, 1.0); return u; }, order);
    REQUIRE(is_permutation(order, 3));
  }
}


TEST_CASE("cone_order", "[accumulation]") {
  // Seed at the origin; sites east (in cone), far east, and north (out of cone)
  std::vector<double> x = {0.0, 1.0, 3.0, 0.0};
  std::vector<double> y = {0.0, 0.0, 0.0, 1.0};
  std::vector<int> order;
  cone_order(x, y, 0, 0.0, PI / 4, order);

  REQUIRE(order == std::vector<int>{0, 1, 2, 3});
  REQUIRE(is_permutation(order, 4));
}


TEST_CASE("richness_within_radius", "[accumulation]") {
  SiteSpeciesMatrix species = {{1, 0, 0}, {1, 1, 0}, {0, 1, 1}};
  DistanceMatrix d = {
    {0.0, 1.0, 2.0},
    {1.0, 0.0, 1.0},
    {2.0, 1.0, 0.0}
  };
  Species sp{species};
  Dist dist{d};

  auto r0 = richness_within_radius(sp, dist, 3, 3, 0, 0.0);
  REQUIRE(r0.first == 1);
  REQUIRE(r0.second == 1);

  auto r1 = richness_within_radius(sp, dist, 3, 3, 0, 1.0);
  REQUIRE(r1.first == 2);
  REQUIRE(r1.second == 2);

  auto r2 = richness_within_radius(sp, dist, 3, 3, 0, 2.0);
  REQUIRE(r2.first == 3);
  REQUIRE(r2.second == 3);
}
