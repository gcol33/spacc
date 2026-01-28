// test_hill.cpp
// Catch2 tests for spacc Hill number calculations

#include "catch.hpp"
#include "../core/hill_core.h"
#include <cmath>

using namespace spacc;

TEST_CASE("calc_hill_number q=0 (species richness)", "[hill]") {
  SECTION("empty community") {
    std::vector<double> abundances = {};
    REQUIRE(calc_hill_number(abundances, 0.0) == Approx(0.0));
  }

  SECTION("all zeros") {
    std::vector<double> abundances = {0, 0, 0, 0};
    REQUIRE(calc_hill_number(abundances, 0.0) == Approx(0.0));
  }

  SECTION("single species") {
    std::vector<double> abundances = {10};
    REQUIRE(calc_hill_number(abundances, 0.0) == Approx(1.0));
  }

  SECTION("multiple species") {
    std::vector<double> abundances = {5, 10, 3, 0, 8};
    REQUIRE(calc_hill_number(abundances, 0.0) == Approx(4.0));  // 4 species present
  }

  SECTION("presence-absence data") {
    std::vector<int> pa = {1, 1, 1, 0, 1, 1};
    REQUIRE(calc_hill_number(pa, 0.0) == Approx(5.0));
  }
}

TEST_CASE("calc_hill_number q=1 (Shannon exponential)", "[hill]") {
  SECTION("single species") {
    std::vector<double> abundances = {100};
    REQUIRE(calc_hill_number(abundances, 1.0) == Approx(1.0));
  }

  SECTION("two equal species") {
    std::vector<double> abundances = {50, 50};
    // H = -2*(0.5*ln(0.5)) = ln(2)
    // exp(H) = 2
    REQUIRE(calc_hill_number(abundances, 1.0) == Approx(2.0));
  }

  SECTION("three equal species") {
    std::vector<double> abundances = {10, 10, 10};
    REQUIRE(calc_hill_number(abundances, 1.0) == Approx(3.0));
  }

  SECTION("unequal abundances") {
    std::vector<double> abundances = {90, 10};
    // More dominant = lower diversity
    double hill = calc_hill_number(abundances, 1.0);
    REQUIRE(hill < 2.0);
    REQUIRE(hill > 1.0);
  }

  SECTION("highly skewed community") {
    std::vector<double> abundances = {1000, 1, 1, 1, 1};
    double hill = calc_hill_number(abundances, 1.0);
    // Should be close to 1 (dominated by single species)
    REQUIRE(hill < 2.0);
    REQUIRE(hill >= 1.0);
  }
}

TEST_CASE("calc_hill_number q=2 (inverse Simpson)", "[hill]") {
  SECTION("single species") {
    std::vector<double> abundances = {100};
    REQUIRE(calc_hill_number(abundances, 2.0) == Approx(1.0));
  }

  SECTION("two equal species") {
    std::vector<double> abundances = {50, 50};
    // D = 2*(0.5)^2 = 0.5
    // 1/D = 2
    REQUIRE(calc_hill_number(abundances, 2.0) == Approx(2.0));
  }

  SECTION("four equal species") {
    std::vector<double> abundances = {25, 25, 25, 25};
    REQUIRE(calc_hill_number(abundances, 2.0) == Approx(4.0));
  }

  SECTION("q=2 <= q=1 <= q=0 for uneven communities") {
    std::vector<double> abundances = {100, 50, 25, 10, 5, 1};
    double q0 = calc_hill_number(abundances, 0.0);
    double q1 = calc_hill_number(abundances, 1.0);
    double q2 = calc_hill_number(abundances, 2.0);

    REQUIRE(q2 <= q1);
    REQUIRE(q1 <= q0);
  }
}

TEST_CASE("calc_hill_number general q", "[hill]") {
  SECTION("q=0.5") {
    std::vector<double> abundances = {50, 30, 20};
    double q05 = calc_hill_number(abundances, 0.5);
    double q0 = calc_hill_number(abundances, 0.0);
    double q1 = calc_hill_number(abundances, 1.0);

    REQUIRE(q05 < q0);
    REQUIRE(q05 > q1);
  }

  SECTION("q=3 (emphasizes dominance more)") {
    std::vector<double> abundances = {100, 10, 1};
    double q2 = calc_hill_number(abundances, 2.0);
    double q3 = calc_hill_number(abundances, 3.0);

    REQUIRE(q3 <= q2);
  }

  SECTION("negative q not typical but should work") {
    std::vector<double> abundances = {50, 30, 20};
    double q_neg = calc_hill_number(abundances, -1.0);
    double q0 = calc_hill_number(abundances, 0.0);

    // Negative q emphasizes rare species
    REQUIRE(q_neg >= q0);
  }
}

TEST_CASE("calc_shannon_entropy", "[hill]") {
  SECTION("single species") {
    std::vector<double> abundances = {100};
    REQUIRE(calc_shannon_entropy(abundances) == Approx(0.0));
  }

  SECTION("two equal species") {
    std::vector<double> abundances = {50, 50};
    REQUIRE(calc_shannon_entropy(abundances) == Approx(std::log(2.0)));
  }

  SECTION("empty community") {
    std::vector<double> abundances = {};
    REQUIRE(calc_shannon_entropy(abundances) == Approx(0.0));
  }

  SECTION("relation to Hill q=1") {
    std::vector<double> abundances = {40, 30, 20, 10};
    double H = calc_shannon_entropy(abundances);
    double hill_q1 = calc_hill_number(abundances, 1.0);
    REQUIRE(std::exp(H) == Approx(hill_q1));
  }
}

TEST_CASE("calc_simpson_index", "[hill]") {
  SECTION("single species") {
    std::vector<double> abundances = {100};
    REQUIRE(calc_simpson_index(abundances) == Approx(1.0));
  }

  SECTION("two equal species") {
    std::vector<double> abundances = {50, 50};
    REQUIRE(calc_simpson_index(abundances) == Approx(0.5));
  }

  SECTION("four equal species") {
    std::vector<double> abundances = {25, 25, 25, 25};
    REQUIRE(calc_simpson_index(abundances) == Approx(0.25));
  }

  SECTION("empty community") {
    std::vector<double> abundances = {};
    REQUIRE(calc_simpson_index(abundances) == Approx(0.0));
  }
}

TEST_CASE("calc_inverse_simpson", "[hill]") {
  SECTION("relation to Hill q=2") {
    std::vector<double> abundances = {40, 30, 20, 10};
    double inv_simp = calc_inverse_simpson(abundances);
    double hill_q2 = calc_hill_number(abundances, 2.0);
    REQUIRE(inv_simp == Approx(hill_q2));
  }

  SECTION("empty community") {
    std::vector<double> abundances = {};
    REQUIRE(calc_inverse_simpson(abundances) == Approx(0.0));
  }
}

TEST_CASE("count_species_richness", "[hill]") {
  SECTION("with zeros") {
    std::vector<int> abundances = {5, 0, 3, 0, 0, 8, 1};
    REQUIRE(count_species_richness(abundances) == 4);
  }

  SECTION("all present") {
    std::vector<int> abundances = {1, 2, 3, 4, 5};
    REQUIRE(count_species_richness(abundances) == 5);
  }

  SECTION("all absent") {
    std::vector<int> abundances = {0, 0, 0};
    REQUIRE(count_species_richness(abundances) == 0);
  }

  SECTION("empty") {
    std::vector<int> abundances = {};
    REQUIRE(count_species_richness(abundances) == 0);
  }
}

TEST_CASE("Hill numbers with integer vector", "[hill]") {
  SECTION("works with int type") {
    std::vector<int> abundances = {10, 20, 30, 40};
    REQUIRE(calc_hill_number(abundances, 0.0) == Approx(4.0));
    REQUIRE(calc_hill_number(abundances, 1.0) > 0);
    REQUIRE(calc_hill_number(abundances, 2.0) > 0);
  }
}
