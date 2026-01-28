// test_distance.cpp
// Catch2 tests for spacc distance calculations

#include "catch.hpp"
#include "../core/distance_core.h"
#include <cmath>

using namespace spacc;

TEST_CASE("to_radians converts correctly", "[distance]") {
  REQUIRE(to_radians(0.0) == Approx(0.0));
  REQUIRE(to_radians(180.0) == Approx(PI));
  REQUIRE(to_radians(90.0) == Approx(PI / 2));
  REQUIRE(to_radians(360.0) == Approx(2 * PI));
  REQUIRE(to_radians(-90.0) == Approx(-PI / 2));
}

TEST_CASE("euclidean_distance computes correctly", "[distance]") {
  SECTION("same point") {
    REQUIRE(euclidean_distance(0, 0, 0, 0) == Approx(0.0));
    REQUIRE(euclidean_distance(5, 5, 5, 5) == Approx(0.0));
  }

  SECTION("unit distance") {
    REQUIRE(euclidean_distance(0, 0, 1, 0) == Approx(1.0));
    REQUIRE(euclidean_distance(0, 0, 0, 1) == Approx(1.0));
  }

  SECTION("3-4-5 triangle") {
    REQUIRE(euclidean_distance(0, 0, 3, 4) == Approx(5.0));
  }

  SECTION("diagonal") {
    REQUIRE(euclidean_distance(0, 0, 1, 1) == Approx(std::sqrt(2.0)));
  }

  SECTION("negative coordinates") {
    REQUIRE(euclidean_distance(-1, -1, 1, 1) == Approx(std::sqrt(8.0)));
  }

  SECTION("symmetry") {
    REQUIRE(euclidean_distance(1, 2, 3, 4) == Approx(euclidean_distance(3, 4, 1, 2)));
  }
}

TEST_CASE("haversine_distance computes correctly", "[distance]") {
  SECTION("same point") {
    REQUIRE(haversine_distance(0, 0, 0, 0) == Approx(0.0));
    REQUIRE(haversine_distance(48.8566, 2.3522, 48.8566, 2.3522) == Approx(0.0));
  }

  SECTION("equator points") {
    // 1 degree longitude at equator ~ 111 km
    double d = haversine_distance(0, 0, 0, 1);
    REQUIRE(d == Approx(111.19).epsilon(0.01));
  }

  SECTION("London to Paris") {
    // London: 51.5074 N, 0.1278 W
    // Paris: 48.8566 N, 2.3522 E
    double d = haversine_distance(51.5074, -0.1278, 48.8566, 2.3522);
    // Known distance: ~344 km
    REQUIRE(d == Approx(344.0).epsilon(0.05));
  }

  SECTION("antipodal points") {
    // Opposite sides of Earth ~ 20015 km (half circumference)
    double d = haversine_distance(0, 0, 0, 180);
    REQUIRE(d == Approx(EARTH_RADIUS_KM * PI).epsilon(0.01));
  }

  SECTION("symmetry") {
    double d1 = haversine_distance(40, -74, 51, 0);
    double d2 = haversine_distance(51, 0, 40, -74);
    REQUIRE(d1 == Approx(d2));
  }

  SECTION("pole to equator") {
    // North pole to equator = 1/4 circumference ~ 10008 km
    double d = haversine_distance(90, 0, 0, 0);
    REQUIRE(d == Approx(EARTH_RADIUS_KM * PI / 2).epsilon(0.01));
  }
}

TEST_CASE("compute_distance_matrix works correctly", "[distance]") {
  SECTION("single point") {
    std::vector<double> x = {0.0};
    std::vector<double> y = {0.0};
    auto dm = compute_distance_matrix(x, y);
    REQUIRE(dm.size() == 1);
    REQUIRE(dm[0][0] == Approx(0.0));
  }

  SECTION("two points euclidean") {
    std::vector<double> x = {0.0, 3.0};
    std::vector<double> y = {0.0, 4.0};
    auto dm = compute_distance_matrix(x, y, "euclidean");

    REQUIRE(dm.size() == 2);
    REQUIRE(dm[0][0] == Approx(0.0));
    REQUIRE(dm[1][1] == Approx(0.0));
    REQUIRE(dm[0][1] == Approx(5.0));
    REQUIRE(dm[1][0] == Approx(5.0));  // Symmetry
  }

  SECTION("three points") {
    std::vector<double> x = {0.0, 1.0, 0.0};
    std::vector<double> y = {0.0, 0.0, 1.0};
    auto dm = compute_distance_matrix(x, y);

    REQUIRE(dm.size() == 3);
    // Diagonal is zero
    REQUIRE(dm[0][0] == Approx(0.0));
    REQUIRE(dm[1][1] == Approx(0.0));
    REQUIRE(dm[2][2] == Approx(0.0));
    // (0,0) to (1,0)
    REQUIRE(dm[0][1] == Approx(1.0));
    // (0,0) to (0,1)
    REQUIRE(dm[0][2] == Approx(1.0));
    // (1,0) to (0,1)
    REQUIRE(dm[1][2] == Approx(std::sqrt(2.0)));
  }

  SECTION("haversine method") {
    std::vector<double> lon = {0.0, 1.0};
    std::vector<double> lat = {0.0, 0.0};
    auto dm = compute_distance_matrix(lon, lat, "haversine");

    REQUIRE(dm.size() == 2);
    REQUIRE(dm[0][0] == Approx(0.0));
    REQUIRE(dm[0][1] == Approx(111.19).epsilon(0.01));
  }
}

TEST_CASE("find_nearest_unvisited works correctly", "[distance]") {
  SECTION("simple case") {
    DistanceMatrix dm = {
      {0.0, 1.0, 2.0},
      {1.0, 0.0, 1.5},
      {2.0, 1.5, 0.0}
    };
    VisitedFlags visited = {true, false, false};

    int nearest = find_nearest_unvisited(dm, 0, visited);
    REQUIRE(nearest == 1);  // Distance 1.0 < 2.0
  }

  SECTION("all visited except one") {
    DistanceMatrix dm = {
      {0.0, 1.0, 2.0},
      {1.0, 0.0, 1.5},
      {2.0, 1.5, 0.0}
    };
    VisitedFlags visited = {true, true, false};

    int nearest = find_nearest_unvisited(dm, 0, visited);
    REQUIRE(nearest == 2);
  }

  SECTION("returns -1 when all visited") {
    DistanceMatrix dm = {
      {0.0, 1.0},
      {1.0, 0.0}
    };
    VisitedFlags visited = {true, true};

    int nearest = find_nearest_unvisited(dm, 0, visited);
    REQUIRE(nearest == -1);
  }

  SECTION("ties broken by index order") {
    DistanceMatrix dm = {
      {0.0, 1.0, 1.0},
      {1.0, 0.0, 1.0},
      {1.0, 1.0, 0.0}
    };
    VisitedFlags visited = {true, false, false};

    int nearest = find_nearest_unvisited(dm, 0, visited);
    REQUIRE(nearest == 1);  // First with min distance
  }
}
