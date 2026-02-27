// test_coverage.cpp
// Catch2 tests for spacc coverage calculations

#include "catch.hpp"
#include "../core/coverage_core.h"

using namespace spacc;

TEST_CASE("count_singletons", "[coverage]") {
  SECTION("no singletons") {
    std::vector<int> abundances = {5, 3, 2, 10};
    REQUIRE(count_singletons(abundances) == 0);
  }

  SECTION("all singletons") {
    std::vector<int> abundances = {1, 1, 1, 1};
    REQUIRE(count_singletons(abundances) == 4);
  }

  SECTION("mixed") {
    std::vector<int> abundances = {1, 5, 1, 3, 1, 2};
    REQUIRE(count_singletons(abundances) == 3);
  }

  SECTION("empty") {
    std::vector<int> abundances = {};
    REQUIRE(count_singletons(abundances) == 0);
  }

  SECTION("with zeros") {
    std::vector<int> abundances = {0, 1, 0, 1, 0};
    REQUIRE(count_singletons(abundances) == 2);
  }
}

TEST_CASE("count_doubletons", "[coverage]") {
  SECTION("no doubletons") {
    std::vector<int> abundances = {5, 3, 1, 10};
    REQUIRE(count_doubletons(abundances) == 0);
  }

  SECTION("all doubletons") {
    std::vector<int> abundances = {2, 2, 2};
    REQUIRE(count_doubletons(abundances) == 3);
  }

  SECTION("mixed") {
    std::vector<int> abundances = {1, 2, 3, 2, 1, 2};
    REQUIRE(count_doubletons(abundances) == 3);
  }
}

TEST_CASE("count_total_individuals", "[coverage]") {
  SECTION("normal case") {
    std::vector<int> abundances = {5, 3, 2, 10};
    REQUIRE(count_total_individuals(abundances) == 20);
  }

  SECTION("empty") {
    std::vector<int> abundances = {};
    REQUIRE(count_total_individuals(abundances) == 0);
  }

  SECTION("with zeros") {
    std::vector<int> abundances = {0, 5, 0, 3, 0};
    REQUIRE(count_total_individuals(abundances) == 8);
  }
}

TEST_CASE("calc_chao_coverage", "[coverage]") {
  SECTION("empty community") {
    std::vector<int> abundances = {};
    REQUIRE(calc_chao_coverage(abundances) == Approx(0.0));
  }

  SECTION("single individual") {
    std::vector<int> abundances = {1};
    REQUIRE(calc_chao_coverage(abundances) == Approx(0.0));
  }

  SECTION("no singletons or doubletons") {
    std::vector<int> abundances = {10, 20, 30};
    REQUIRE(calc_chao_coverage(abundances) == Approx(1.0));
  }

  SECTION("only singletons") {
    std::vector<int> abundances = {1, 1, 1, 1, 1};
    double C = calc_chao_coverage(abundances);
    // n=5, f1=5, f2=0
    // C = 1 - (5/5) * (4/5) = 1 - 0.8 = 0.2
    REQUIRE(C == Approx(0.2));
  }

  SECTION("with doubletons") {
    // n=10, f1=2, f2=2
    // Formula: C = 1 - (f1/n) * ((n-1)*f1 / ((n-1)*f1 + 2*f2))
    // C = 1 - (2/10) * (9*2 / (9*2 + 4))
    // C = 1 - 0.2 * (18/22) = 1 - 0.2 * 0.818 = 1 - 0.1636 = 0.836
    std::vector<int> abundances = {1, 1, 2, 2, 4};  // n=10, f1=2, f2=2
    double C = calc_chao_coverage(abundances);
    REQUIRE(C == Approx(0.8364).epsilon(0.01));
  }

  SECTION("high coverage community") {
    std::vector<int> abundances = {100, 50, 30, 20, 10};
    double C = calc_chao_coverage(abundances);
    REQUIRE(C > 0.9);
    REQUIRE(C <= 1.0);
  }

  SECTION("coverage is between 0 and 1") {
    std::vector<int> abundances = {1, 1, 1, 2, 2, 5, 10};
    double C = calc_chao_coverage(abundances);
    REQUIRE(C >= 0.0);
    REQUIRE(C <= 1.0);
  }

  SECTION("more individuals = higher coverage") {
    std::vector<int> small = {1, 1, 1, 2, 3};
    std::vector<int> large = {10, 10, 10, 20, 30};

    double C_small = calc_chao_coverage(small);
    double C_large = calc_chao_coverage(large);

    REQUIRE(C_large >= C_small);
  }
}

TEST_CASE("calc_coverage_deficit", "[coverage]") {
  SECTION("complement of coverage") {
    std::vector<int> abundances = {1, 1, 2, 5, 10};
    double C = calc_chao_coverage(abundances);
    double D = calc_coverage_deficit(abundances);
    REQUIRE(C + D == Approx(1.0));
  }

  SECTION("full coverage has zero deficit") {
    std::vector<int> abundances = {10, 20, 30};
    REQUIRE(calc_coverage_deficit(abundances) == Approx(0.0));
  }

  SECTION("empty has full deficit... wait no") {
    std::vector<int> abundances = {};
    // Empty returns 0 coverage, so deficit = 1
    REQUIRE(calc_coverage_deficit(abundances) == Approx(1.0));
  }
}

TEST_CASE("interpolate_richness_at_coverage", "[coverage]") {
  SECTION("target at exact point") {
    std::vector<double> richness = {1, 2, 3, 4, 5};
    std::vector<double> coverage = {0.2, 0.4, 0.6, 0.8, 1.0};

    REQUIRE(interpolate_richness_at_coverage(richness, coverage, 0.4) == Approx(2.0));
    REQUIRE(interpolate_richness_at_coverage(richness, coverage, 0.8) == Approx(4.0));
  }

  SECTION("linear interpolation") {
    std::vector<double> richness = {0, 10};
    std::vector<double> coverage = {0, 1.0};

    REQUIRE(interpolate_richness_at_coverage(richness, coverage, 0.5) == Approx(5.0));
    REQUIRE(interpolate_richness_at_coverage(richness, coverage, 0.25) == Approx(2.5));
  }

  SECTION("target below minimum") {
    std::vector<double> richness = {5, 10, 15};
    std::vector<double> coverage = {0.3, 0.6, 0.9};

    REQUIRE(interpolate_richness_at_coverage(richness, coverage, 0.1) == Approx(5.0));
  }

  SECTION("target above maximum returns -1") {
    std::vector<double> richness = {5, 10, 15};
    std::vector<double> coverage = {0.3, 0.6, 0.9};

    REQUIRE(interpolate_richness_at_coverage(richness, coverage, 0.95) == Approx(-1.0));
  }

  SECTION("empty vectors") {
    std::vector<double> richness = {};
    std::vector<double> coverage = {};

    REQUIRE(interpolate_richness_at_coverage(richness, coverage, 0.5) == Approx(0.0));
  }

  SECTION("non-linear accumulation") {
    std::vector<double> richness = {5, 15, 20, 22, 23};
    std::vector<double> coverage = {0.2, 0.5, 0.7, 0.9, 1.0};

    // Interpolate between 15 and 20 at coverage 0.6
    // c0=0.5, c1=0.7, r0=15, r1=20
    // r = 15 + (0.6-0.5)*(20-15)/(0.7-0.5) = 15 + 0.1*5/0.2 = 15 + 2.5 = 17.5
    REQUIRE(interpolate_richness_at_coverage(richness, coverage, 0.6) == Approx(17.5));
  }
}

TEST_CASE("coverage with double abundances", "[coverage]") {
  SECTION("works with double type") {
    std::vector<double> abundances = {1.0, 1.0, 2.0, 5.0};
    double C = calc_chao_coverage(abundances);
    REQUIRE(C >= 0.0);
    REQUIRE(C <= 1.0);
  }
}

TEST_CASE("calc_chiu_coverage", "[coverage]") {
  SECTION("empty community") {
    std::vector<int> abundances = {};
    std::vector<int> incidences = {};
    REQUIRE(calc_chiu_coverage(abundances, incidences, 0) == Approx(0.0));
  }

  SECTION("single site returns 0") {
    std::vector<int> abundances = {5, 3, 1};
    std::vector<int> incidences = {1, 1, 1};
    REQUIRE(calc_chiu_coverage(abundances, incidences, 1) == Approx(0.0));
  }

  SECTION("no uniques gives full coverage") {
    // All species in >= 2 sites -> Q1 = 0 -> coverage = 1
    std::vector<int> abundances = {10, 20, 30};
    std::vector<int> incidences = {3, 5, 4};
    REQUIRE(calc_chiu_coverage(abundances, incidences, 5) == Approx(1.0));
  }

  SECTION("known values") {
    // T=10, n+=100
    // sp A: abundance=5, in 1 site  -> Q1 species, contributes to G1
    // sp B: abundance=3, in 1 site  -> Q1 species, contributes to G1
    // sp C: abundance=20, in 2 sites -> Q2 species
    // sp D: abundance=30, in 5 sites
    // sp E: abundance=42, in 8 sites
    // Q1=2, Q2=1, G1=5+3=8, n+=100, T=10
    // a0_hat = (9/10) * 8 * 2 / (2*1) = 0.9 * 16 / 2 = 7.2
    // C = 100 / (100 + 7.2) = 100 / 107.2 ~ 0.9328
    std::vector<int> abundances = {5, 3, 20, 30, 42};
    std::vector<int> incidences = {1, 1, 2, 5, 8};
    double C = calc_chiu_coverage(abundances, incidences, 10);
    REQUIRE(C == Approx(0.9328).epsilon(0.001));
  }

  SECTION("all species unique (Q2 = 0 uses Q2_safe = 1)") {
    // All species in exactly 1 site
    // Q1=4, Q2=0 (capped to 1), G1=1+1+1+1=4, n+=4, T=4
    // a0_hat = (3/4) * 4 * 4 / (2*1) = 0.75 * 16 / 2 = 6
    // C = 4 / (4 + 6) = 0.4
    std::vector<int> abundances = {1, 1, 1, 1};
    std::vector<int> incidences = {1, 1, 1, 1};
    double C = calc_chiu_coverage(abundances, incidences, 4);
    REQUIRE(C == Approx(0.4).epsilon(0.001));
  }

  SECTION("coverage is between 0 and 1") {
    std::vector<int> abundances = {1, 2, 5, 10, 3, 1};
    std::vector<int> incidences = {1, 2, 3, 5, 1, 2};
    double C = calc_chiu_coverage(abundances, incidences, 5);
    REQUIRE(C >= 0.0);
    REQUIRE(C <= 1.0);
  }

  SECTION("more sites generally increases coverage") {
    // With 3 sites
    std::vector<int> abund3 = {2, 1, 1, 3};
    std::vector<int> incid3 = {1, 1, 2, 2};
    double C3 = calc_chiu_coverage(abund3, incid3, 3);

    // Same species, but accumulated over 10 sites (more sites, higher incidences)
    std::vector<int> abund10 = {20, 10, 15, 30};
    std::vector<int> incid10 = {5, 3, 7, 9};
    double C10 = calc_chiu_coverage(abund10, incid10, 10);

    REQUIRE(C10 >= C3);
  }

  SECTION("species with zero abundance are ignored") {
    std::vector<int> abundances = {0, 5, 0, 10};
    std::vector<int> incidences = {0, 3, 0, 5};
    double C = calc_chiu_coverage(abundances, incidences, 5);
    // Q1=0 -> coverage = 1.0
    REQUIRE(C == Approx(1.0));
  }
}
