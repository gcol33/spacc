// test_beta.cpp
// Catch2 tests for spacc beta diversity calculations

#include "catch.hpp"
#include "../core/beta_core.h"

using namespace spacc;

TEST_CASE("calc_beta_sorensen computes correctly", "[beta]") {
  SECTION("identical communities (a=10, b=0, c=0)") {
    auto result = calc_beta_sorensen(10, 0, 0);
    REQUIRE(result.total == Approx(0.0));
    REQUIRE(result.turnover == Approx(0.0));
    REQUIRE(result.nestedness == Approx(0.0));
  }

  SECTION("completely different communities (a=0, b=5, c=5)") {
    auto result = calc_beta_sorensen(0, 5, 5);
    REQUIRE(result.total == Approx(1.0));
    REQUIRE(result.turnover == Approx(1.0));
    REQUIRE(result.nestedness == Approx(0.0));
  }

  SECTION("partial overlap (a=5, b=3, c=2)") {
    auto result = calc_beta_sorensen(5, 3, 2);
    // total = (3+2)/(2*5+3+2) = 5/15 = 0.333
    REQUIRE(result.total == Approx(5.0 / 15.0));
    // turnover = min(3,2)/(5+min(3,2)) = 2/7
    REQUIRE(result.turnover == Approx(2.0 / 7.0));
    // nestedness = total - turnover
    REQUIRE(result.nestedness == Approx(result.total - result.turnover));
  }

  SECTION("nested community (a=5, b=5, c=0)") {
    auto result = calc_beta_sorensen(5, 5, 0);
    // total = 5/15 = 0.333
    REQUIRE(result.total == Approx(5.0 / 15.0));
    // turnover = 0/5 = 0 (min_bc = 0)
    REQUIRE(result.turnover == Approx(0.0));
    // All beta is due to nestedness
    REQUIRE(result.nestedness == Approx(result.total));
  }

  SECTION("empty communities") {
    auto result = calc_beta_sorensen(0, 0, 0);
    REQUIRE(result.total == Approx(0.0));
    REQUIRE(result.turnover == Approx(0.0));
    REQUIRE(result.nestedness == Approx(0.0));
  }

  SECTION("turnover + nestedness equals total") {
    auto result = calc_beta_sorensen(10, 7, 3);
    REQUIRE(result.total == Approx(result.turnover + result.nestedness));
  }
}

TEST_CASE("calc_beta_jaccard computes correctly", "[beta]") {
  SECTION("identical communities") {
    auto result = calc_beta_jaccard(10, 0, 0);
    REQUIRE(result.total == Approx(0.0));
    REQUIRE(result.turnover == Approx(0.0));
    REQUIRE(result.nestedness == Approx(0.0));
  }

  SECTION("completely different communities") {
    auto result = calc_beta_jaccard(0, 5, 5);
    REQUIRE(result.total == Approx(1.0));
    REQUIRE(result.turnover == Approx(1.0));
    REQUIRE(result.nestedness == Approx(0.0));
  }

  SECTION("partial overlap (a=5, b=3, c=2)") {
    auto result = calc_beta_jaccard(5, 3, 2);
    // total = (3+2)/(5+3+2) = 5/10 = 0.5
    REQUIRE(result.total == Approx(0.5));
    // turnover = 2*min(3,2)/(5+2*min(3,2)) = 4/9
    REQUIRE(result.turnover == Approx(4.0 / 9.0));
    REQUIRE(result.nestedness == Approx(result.total - result.turnover));
  }

  SECTION("empty communities") {
    auto result = calc_beta_jaccard(0, 0, 0);
    REQUIRE(result.total == Approx(0.0));
  }

  SECTION("turnover + nestedness equals total") {
    auto result = calc_beta_jaccard(8, 4, 6);
    REQUIRE(result.total == Approx(result.turnover + result.nestedness));
  }
}

TEST_CASE("Sorensen vs Jaccard relationship", "[beta]") {
  SECTION("Jaccard total >= Sorensen total") {
    for (int a = 0; a <= 10; a++) {
      for (int b = 0; b <= 10; b++) {
        for (int c = 0; c <= 10; c++) {
          if (a + b + c > 0) {
            auto sor = calc_beta_sorensen(a, b, c);
            auto jac = calc_beta_jaccard(a, b, c);
            REQUIRE(jac.total >= sor.total - 1e-10);
          }
        }
      }
    }
  }
}

TEST_CASE("count_abc works correctly", "[beta]") {
  SECTION("identical sets") {
    SpeciesSet s1 = {1, 2, 3, 4, 5};
    SpeciesSet s2 = {1, 2, 3, 4, 5};
    auto [a, b, c] = count_abc(s1, s2);
    REQUIRE(a == 5);
    REQUIRE(b == 0);
    REQUIRE(c == 0);
  }

  SECTION("disjoint sets") {
    SpeciesSet s1 = {1, 2, 3};
    SpeciesSet s2 = {4, 5, 6};
    auto [a, b, c] = count_abc(s1, s2);
    REQUIRE(a == 0);
    REQUIRE(b == 3);
    REQUIRE(c == 3);
  }

  SECTION("partial overlap") {
    SpeciesSet s1 = {1, 2, 3, 4};
    SpeciesSet s2 = {3, 4, 5, 6};
    auto [a, b, c] = count_abc(s1, s2);
    REQUIRE(a == 2);  // {3, 4}
    REQUIRE(b == 2);  // {1, 2}
    REQUIRE(c == 2);  // {5, 6}
  }

  SECTION("empty sets") {
    SpeciesSet s1 = {};
    SpeciesSet s2 = {};
    auto [a, b, c] = count_abc(s1, s2);
    REQUIRE(a == 0);
    REQUIRE(b == 0);
    REQUIRE(c == 0);
  }

  SECTION("one empty set") {
    SpeciesSet s1 = {1, 2, 3};
    SpeciesSet s2 = {};
    auto [a, b, c] = count_abc(s1, s2);
    REQUIRE(a == 0);
    REQUIRE(b == 3);
    REQUIRE(c == 0);
  }

  SECTION("subset relationship") {
    SpeciesSet s1 = {1, 2};
    SpeciesSet s2 = {1, 2, 3, 4, 5};
    auto [a, b, c] = count_abc(s1, s2);
    REQUIRE(a == 2);
    REQUIRE(b == 0);
    REQUIRE(c == 3);
  }
}

TEST_CASE("count_abc_vectors works correctly", "[beta]") {
  SECTION("with duplicates in input") {
    std::vector<int> v1 = {1, 2, 2, 3, 3, 3};  // Unique: {1, 2, 3}
    std::vector<int> v2 = {2, 3, 4};            // Unique: {2, 3, 4}
    auto [a, b, c] = count_abc_vectors(v1, v2);
    REQUIRE(a == 2);  // {2, 3}
    REQUIRE(b == 1);  // {1}
    REQUIRE(c == 1);  // {4}
  }
}
