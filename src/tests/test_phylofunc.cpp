// test_phylofunc.cpp
// Catch2 tests for the phylogenetic and functional diversity kernels

#include "catch.hpp"
#include "../core/phylo_core.h"
#include "../core/func_core.h"
#include <cmath>
#include <vector>

using namespace spacc;

namespace {

struct Mat {
  std::vector<std::vector<double>> m;
  double operator()(int i, int j) const { return m[i][j]; }
};

double one(int) { return 1.0; }

} // namespace


TEST_CASE("present_species", "[phylo]") {
  std::vector<double> cumulative = {0.0, 2.5, 0.0, 1.0};
  REQUIRE(present_species(cumulative) == std::vector<int>{1, 3});
  REQUIRE(present_species(std::vector<int>{0, 0}).empty());
}


TEST_CASE("mean_pairwise_distance", "[phylo]") {
  Mat d{{{0, 1, 2}, {1, 0, 3}, {2, 3, 0}}};

  SECTION("fewer than two species gives zero") {
    REQUIRE(mean_pairwise_distance(d, {0}, one) == Approx(0.0));
  }

  SECTION("unweighted mean over pairs") {
    REQUIRE(mean_pairwise_distance(d, {0, 1, 2}, one) == Approx(2.0));
  }

  SECTION("abundance weights") {
    std::vector<double> w = {1.0, 3.0, 0.0};
    auto weight = [&](int i) { return w[i]; };
    // pairs (0,1) weight 3 dist 1, (0,2) weight 0, (1,2) weight 0
    REQUIRE(mean_pairwise_distance(d, {0, 1, 2}, weight) == Approx(1.0));
  }
}


TEST_CASE("mean_nearest_taxon_distance", "[phylo]") {
  Mat d{{{0, 1, 5}, {1, 0, 4}, {5, 4, 0}}};
  // nearest: 0 -> 1 (1), 1 -> 0 (1), 2 -> 1 (4)
  REQUIRE(mean_nearest_taxon_distance(d, {0, 1, 2}, one) == Approx(2.0));
  REQUIRE(mean_nearest_taxon_distance(d, {2}, one) == Approx(0.0));
}


TEST_CASE("rao_quadratic_entropy", "[phylo]") {
  Mat d{{{0, 2}, {2, 0}}};
  std::vector<double> ab = {1.0, 1.0};
  auto abundance = [&](int i) { return ab[i]; };
  // 2 * p0 p1 d01 = 2 * 0.25 * 2 = 1
  REQUIRE(rao_quadratic_entropy(d, abundance, 2) == Approx(1.0));

  std::vector<double> single = {3.0, 0.0};
  REQUIRE(rao_quadratic_entropy(d, [&](int i) { return single[i]; }, 2) == Approx(0.0));
}


TEST_CASE("PhyloTree::faith_pd", "[phylo]") {
  // ((t1:1, t2:2):3, t3:4); in ape encoding: tips 1..3, root 4, inner 5
  // edges: 4->5 (3), 5->1 (1), 5->2 (2), 4->3 (4)
  std::vector<int> parent = {4, 5, 5, 4};
  std::vector<int> child = {5, 1, 2, 3};
  std::vector<double> len = {3, 1, 2, 4};
  PhyloTree tree(4,
    [&](int e) { return parent[e]; },
    [&](int e) { return child[e]; },
    [&](int e) { return len[e]; },
    3);

  auto pd = [&](std::vector<bool> present) {
    return tree.faith_pd([&](int tip) { return present[tip]; });
  };

  REQUIRE(pd({true, false, false}) == Approx(4.0));   // 1 + 3
  REQUIRE(pd({true, true, false}) == Approx(6.0));    // 1 + 2 + 3
  REQUIRE(pd({true, true, true}) == Approx(10.0));    // all edges
  REQUIRE(pd({false, false, true}) == Approx(4.0));
  REQUIRE(pd({false, false, false}) == Approx(0.0));
}


TEST_CASE("PhyloTree::evo_heritage", "[phylo]") {
  // Same tree as above, with the inner node numbered before the root to check
  // that the traversal does not rely on ape's node numbering.
  // ((t1:1, t2:2):3, t3:4); tips 1..3, inner 4, root 5
  std::vector<int> parent = {5, 4, 4, 5};
  std::vector<int> child = {4, 1, 2, 3};
  std::vector<double> len = {3, 1, 2, 4};
  PhyloTree tree(4,
    [&](int e) { return parent[e]; },
    [&](int e) { return child[e]; },
    [&](int e) { return len[e]; },
    3);

  auto evo = [&](std::vector<bool> present, double rho) {
    return tree.evo_heritage([&](int tip) { return present[tip]; }, rho);
  };
  auto pd = [&](std::vector<bool> present) {
    return tree.faith_pd([&](int tip) { return present[tip]; });
  };

  SECTION("rho = 0 is Faith's PD") {
    for (auto set : std::vector<std::vector<bool>>{
           {true, false, false}, {true, true, false}, {true, true, true},
           {false, true, true}, {false, false, false}}) {
      REQUIRE(evo(set, 0.0) == Approx(pd(set)));
    }
    REQUIRE(evo({true, true, false}, 1e-9) == Approx(6.0).epsilon(1e-6));
  }

  SECTION("large rho tends to the number of present tips") {
    REQUIRE(evo({true, false, false}, 50.0) == Approx(1.0));
    REQUIRE(evo({true, true, false}, 50.0) == Approx(2.0));
    REQUIRE(evo({true, true, true}, 50.0) == Approx(3.0));
  }

  SECTION("hand-computed value at rho = 1, tips t1 and t2") {
    // p(t1) = p(t2) = 1; p(inner) = 1 - (1 - e^-1)(1 - e^-2)
    // phi = [(1 - e^-1) + (1 - e^-2) + (1 - e^-3) p(inner)] / (1 - e^-1)
    const double e1 = std::exp(-1.0), e2 = std::exp(-2.0), e3 = std::exp(-3.0);
    const double p_inner = 1.0 - (1.0 - e1) * (1.0 - e2);
    const double expected = ((1.0 - e1) + (1.0 - e2) + (1.0 - e3) * p_inner) / (1.0 - e1);
    REQUIRE(evo({true, true, false}, 1.0) == Approx(expected));
  }

  SECTION("empty set gives zero") {
    REQUIRE(evo({false, false, false}, 0.5) == Approx(0.0));
  }
}


TEST_CASE("functional_dispersion", "[func]") {
  Mat traits{{{0, 0}, {2, 0}, {1, 1}}};

  SECTION("fewer than two species gives zero") {
    REQUIRE(functional_dispersion(traits, 2, {0}, one) == Approx(0.0));
  }

  SECTION("equal weights: mean distance to the centroid") {
    // centroid of (0,0),(2,0) is (1,0); both at distance 1
    REQUIRE(functional_dispersion(traits, 2, {0, 1}, one) == Approx(1.0));
  }

  SECTION("weights move the centroid") {
    std::vector<double> w = {3.0, 1.0, 0.0};
    auto weight = [&](int i) { return w[i]; };
    // centroid (0.5, 0); distances 0.5 and 1.5 weighted 3/4 and 1/4
    REQUIRE(functional_dispersion(traits, 2, {0, 1}, weight) == Approx(0.75));
  }
}


TEST_CASE("functional_range_volume", "[func]") {
  Mat traits{{{0, 0}, {2, 0}, {1, 3}}};
  REQUIRE(functional_range_volume(traits, 2, {0, 1}) == Approx(0.0));      // too few species
  REQUIRE(functional_range_volume(traits, 2, {0, 1, 2}) == Approx(6.0));   // 2 * 3
  Mat flat{{{0, 0}, {2, 0}, {1, 0}}};
  REQUIRE(functional_range_volume(flat, 2, {0, 1, 2}) == Approx(0.0));     // zero range
}


TEST_CASE("rao_trait_entropy", "[func]") {
  Mat traits{{{0, 0}, {3, 4}}};
  std::vector<double> ab = {1.0, 1.0};
  // 2 * 0.25 * 5 = 2.5
  REQUIRE(rao_trait_entropy(traits, 2, 2, [&](int i) { return ab[i]; }) == Approx(2.5));
}
