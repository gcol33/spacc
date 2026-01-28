// test_accumulation.cpp
// Catch2 tests for spacc accumulation algorithms

#include "catch.hpp"
#include "../core/accumulation_core.h"

using namespace spacc;

TEST_CASE("accumulate_species", "[accumulation]") {
  SECTION("adds new species") {
    std::vector<int> species_pa = {1, 0, 1, 1, 0};
    SpeciesSet seen;

    int new_count = accumulate_species(species_pa, seen);

    REQUIRE(new_count == 3);
    REQUIRE(seen.size() == 3);
    REQUIRE(seen.count(0) == 1);
    REQUIRE(seen.count(2) == 1);
    REQUIRE(seen.count(3) == 1);
  }

  SECTION("no duplicates") {
    std::vector<int> species_pa = {1, 1, 0};
    SpeciesSet seen = {0};  // Species 0 already seen

    int new_count = accumulate_species(species_pa, seen);

    REQUIRE(new_count == 1);  // Only species 1 is new
    REQUIRE(seen.size() == 2);
  }

  SECTION("empty site") {
    std::vector<int> species_pa = {0, 0, 0};
    SpeciesSet seen = {0, 1};

    int new_count = accumulate_species(species_pa, seen);

    REQUIRE(new_count == 0);
    REQUIRE(seen.size() == 2);
  }

  SECTION("all already seen") {
    std::vector<int> species_pa = {1, 1, 1};
    SpeciesSet seen = {0, 1, 2};

    int new_count = accumulate_species(species_pa, seen);

    REQUIRE(new_count == 0);
    REQUIRE(seen.size() == 3);
  }
}

TEST_CASE("accumulate_abundances", "[accumulation]") {
  SECTION("adds abundances") {
    std::vector<int> site = {5, 3, 0, 2};
    std::vector<int> cumulative = {10, 0, 5, 0};

    accumulate_abundances(site, cumulative);

    REQUIRE(cumulative[0] == 15);
    REQUIRE(cumulative[1] == 3);
    REQUIRE(cumulative[2] == 5);
    REQUIRE(cumulative[3] == 2);
  }
}

TEST_CASE("knn_accumulate_single", "[accumulation]") {
  SECTION("simple 3-site case") {
    // Sites arranged in a line: 0 -- 1 -- 2
    // Site 0: species A
    // Site 1: species A, B
    // Site 2: species B, C
    SiteSpeciesMatrix species_pa = {
      {1, 0, 0},  // Site 0: species 0
      {1, 1, 0},  // Site 1: species 0, 1
      {0, 1, 1}   // Site 2: species 1, 2
    };

    DistanceMatrix dist_mat = {
      {0.0, 1.0, 2.0},
      {1.0, 0.0, 1.0},
      {2.0, 1.0, 0.0}
    };

    // Start from site 0
    Curve curve = knn_accumulate_single(species_pa, dist_mat, 0);

    REQUIRE(curve.size() == 3);
    REQUIRE(curve[0] == 1);  // Site 0: 1 species
    REQUIRE(curve[1] == 2);  // Site 0+1: 2 species (0 and 1)
    REQUIRE(curve[2] == 3);  // All: 3 species
  }

  SECTION("monotonic increase") {
    SiteSpeciesMatrix species_pa = {
      {1, 0, 0, 0},
      {0, 1, 0, 0},
      {0, 0, 1, 0},
      {0, 0, 0, 1}
    };

    DistanceMatrix dist_mat = {
      {0.0, 1.0, 2.0, 3.0},
      {1.0, 0.0, 1.0, 2.0},
      {2.0, 1.0, 0.0, 1.0},
      {3.0, 2.0, 1.0, 0.0}
    };

    Curve curve = knn_accumulate_single(species_pa, dist_mat, 0);

    for (size_t i = 1; i < curve.size(); i++) {
      REQUIRE(curve[i] >= curve[i-1]);
    }
  }

  SECTION("same species everywhere") {
    SiteSpeciesMatrix species_pa = {
      {1, 1},
      {1, 1},
      {1, 1}
    };

    DistanceMatrix dist_mat = {
      {0.0, 1.0, 2.0},
      {1.0, 0.0, 1.0},
      {2.0, 1.0, 0.0}
    };

    Curve curve = knn_accumulate_single(species_pa, dist_mat, 0);

    // All sites have same 2 species
    REQUIRE(curve[0] == 2);
    REQUIRE(curve[1] == 2);
    REQUIRE(curve[2] == 2);
  }

  SECTION("single site") {
    SiteSpeciesMatrix species_pa = {{1, 1, 0, 1}};
    DistanceMatrix dist_mat = {{0.0}};

    Curve curve = knn_accumulate_single(species_pa, dist_mat, 0);

    REQUIRE(curve.size() == 1);
    REQUIRE(curve[0] == 3);
  }
}

TEST_CASE("kncn_accumulate_single", "[accumulation]") {
  SECTION("simple case") {
    SiteSpeciesMatrix species_pa = {
      {1, 0, 0},
      {0, 1, 0},
      {0, 0, 1}
    };

    // Sites at (0,0), (2,0), (1,2)
    std::vector<double> x = {0, 2, 1};
    std::vector<double> y = {0, 0, 2};

    Curve curve = kncn_accumulate_single(species_pa, x, y, 0);

    REQUIRE(curve.size() == 3);
    REQUIRE(curve[0] == 1);
    // After site 0, centroid is (0,0)
    // Nearest to (0,0): site 1 at (2,0) dist=2, site 2 at (1,2) dist=sqrt(5)
    // So site 1 is next
    REQUIRE(curve[1] == 2);
    REQUIRE(curve[2] == 3);
  }

  SECTION("monotonic increase") {
    SiteSpeciesMatrix species_pa = {
      {1, 0, 0},
      {0, 1, 0},
      {0, 0, 1}
    };

    std::vector<double> x = {0, 1, 2};
    std::vector<double> y = {0, 1, 0};

    Curve curve = kncn_accumulate_single(species_pa, x, y, 0);

    for (size_t i = 1; i < curve.size(); i++) {
      REQUIRE(curve[i] >= curve[i-1]);
    }
  }
}

TEST_CASE("knn_accumulate_multi", "[accumulation]") {
  SECTION("multiple seeds") {
    SiteSpeciesMatrix species_pa = {
      {1, 0},
      {0, 1},
      {1, 1}
    };

    DistanceMatrix dist_mat = {
      {0.0, 1.0, 1.5},
      {1.0, 0.0, 0.5},
      {1.5, 0.5, 0.0}
    };

    std::vector<int> seeds = {0, 1, 2};
    CurveMatrix curves = knn_accumulate_multi(species_pa, dist_mat, seeds);

    REQUIRE(curves.size() == 3);
    REQUIRE(curves[0].size() == 3);
    REQUIRE(curves[1].size() == 3);
    REQUIRE(curves[2].size() == 3);
  }
}

TEST_CASE("mean_curve", "[accumulation]") {
  SECTION("simple average") {
    CurveMatrix curves = {
      {2, 4, 6},
      {4, 6, 8}
    };

    CurveDouble mean = mean_curve(curves);

    REQUIRE(mean.size() == 3);
    REQUIRE(mean[0] == Approx(3.0));
    REQUIRE(mean[1] == Approx(5.0));
    REQUIRE(mean[2] == Approx(7.0));
  }

  SECTION("single curve") {
    CurveMatrix curves = {{1, 2, 3, 4, 5}};

    CurveDouble mean = mean_curve(curves);

    REQUIRE(mean.size() == 5);
    for (size_t i = 0; i < 5; i++) {
      REQUIRE(mean[i] == Approx(static_cast<double>(i + 1)));
    }
  }

  SECTION("empty") {
    CurveMatrix curves = {};
    CurveDouble mean = mean_curve(curves);
    REQUIRE(mean.empty());
  }
}

TEST_CASE("curve_quantile", "[accumulation]") {
  SECTION("median of odd count") {
    CurveMatrix curves = {
      {1, 10, 100},
      {2, 20, 200},
      {3, 30, 300}
    };

    CurveDouble median = curve_quantile(curves, 0.5);

    REQUIRE(median.size() == 3);
    REQUIRE(median[0] == Approx(2.0));
    REQUIRE(median[1] == Approx(20.0));
    REQUIRE(median[2] == Approx(200.0));
  }

  SECTION("extremes") {
    CurveMatrix curves = {
      {1, 1},
      {5, 5},
      {10, 10}
    };

    CurveDouble q0 = curve_quantile(curves, 0.0);
    CurveDouble q1 = curve_quantile(curves, 1.0);

    REQUIRE(q0[0] == Approx(1.0));
    REQUIRE(q1[0] == Approx(10.0));
  }

  SECTION("interpolation at 0.25") {
    CurveMatrix curves = {
      {0},
      {4},
      {8},
      {12},
      {16}
    };

    CurveDouble q25 = curve_quantile(curves, 0.25);
    // Index = 0.25 * 4 = 1
    REQUIRE(q25[0] == Approx(4.0));
  }

  SECTION("empty") {
    CurveMatrix curves = {};
    CurveDouble result = curve_quantile(curves, 0.5);
    REQUIRE(result.empty());
  }
}
