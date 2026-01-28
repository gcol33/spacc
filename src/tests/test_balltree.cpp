// test_balltree.cpp
// Catch2 tests for BallTree (haversine nearest-neighbor)

#include "catch.hpp"
#include "../balltree.h"

TEST_CASE("haversine distance", "[balltree]") {
  SECTION("same point gives zero") {
    double d = haversine(0.0, 0.0, 0.0, 0.0);
    REQUIRE(d == Approx(0.0).margin(1e-10));
  }

  SECTION("known distance Paris to London") {
    // Paris: 2.3522 E, 48.8566 N
    // London: -0.1276 W, 51.5074 N
    double d = haversine(2.3522, 48.8566, -0.1276, 51.5074);
    // Actual distance is about 344 km
    REQUIRE(d == Approx(344.0).margin(10.0));
  }

  SECTION("known distance New York to Los Angeles") {
    // NYC: -74.006, 40.7128
    // LA: -118.2437, 34.0522
    double d = haversine(-74.006, 40.7128, -118.2437, 34.0522);
    // Actual distance is about 3944 km
    REQUIRE(d == Approx(3944.0).margin(50.0));
  }

  SECTION("antipodal points") {
    // 0,0 to 180,0
    double d = haversine(0.0, 0.0, 180.0, 0.0);
    // Half Earth circumference ~ 20015 km
    REQUIRE(d == Approx(20015.0).margin(100.0));
  }

  SECTION("symmetric") {
    double d1 = haversine(10.0, 20.0, 30.0, 40.0);
    double d2 = haversine(30.0, 40.0, 10.0, 20.0);
    REQUIRE(d1 == Approx(d2));
  }

  SECTION("equatorial distance") {
    // 1 degree of longitude at equator ~ 111 km
    double d = haversine(0.0, 0.0, 1.0, 0.0);
    REQUIRE(d == Approx(111.0).margin(5.0));
  }
}

TEST_CASE("BallTree construction", "[balltree]") {
  SECTION("small tree") {
    std::vector<double> x = {0.0, 1.0, 2.0};
    std::vector<double> y = {0.0, 0.0, 0.0};

    BallTree tree(x, y);

    REQUIRE(tree.n_points == 3);
    REQUIRE(!tree.nodes.empty());
  }

  SECTION("larger tree") {
    std::vector<double> x, y;
    for (int i = 0; i < 100; i++) {
      x.push_back(static_cast<double>(i % 10));
      y.push_back(static_cast<double>(i / 10));
    }

    BallTree tree(x, y);

    REQUIRE(tree.n_points == 100);
    REQUIRE(tree.indices.size() == 100);
  }

  SECTION("single point") {
    std::vector<double> x = {5.0};
    std::vector<double> y = {10.0};

    BallTree tree(x, y);

    REQUIRE(tree.n_points == 1);
  }
}

TEST_CASE("BallTree find_nearest_unvisited", "[balltree]") {
  SECTION("simple linear case") {
    // Points at (0,0), (1,0), (2,0) degrees
    std::vector<double> x = {0.0, 1.0, 2.0};
    std::vector<double> y = {0.0, 0.0, 0.0};

    BallTree tree(x, y);
    std::vector<bool> visited(3, false);

    // Query from (0,0) - nearest should be itself (index 0)
    int nearest = tree.find_nearest_unvisited(0.0, 0.0, visited);
    REQUIRE(nearest == 0);

    // Mark index 0 as visited
    visited[0] = true;
    nearest = tree.find_nearest_unvisited(0.0, 0.0, visited);
    REQUIRE(nearest == 1);  // Next nearest is (1,0)

    // Mark index 1 as visited
    visited[1] = true;
    nearest = tree.find_nearest_unvisited(0.0, 0.0, visited);
    REQUIRE(nearest == 2);  // Only (2,0) left

    // Mark all visited
    visited[2] = true;
    nearest = tree.find_nearest_unvisited(0.0, 0.0, visited);
    REQUIRE(nearest == -1);  // None left
  }

  SECTION("off-grid query") {
    std::vector<double> x = {0.0, 10.0, 20.0};
    std::vector<double> y = {0.0, 0.0, 0.0};

    BallTree tree(x, y);
    std::vector<bool> visited(3, false);

    // Query from (11,0) - nearest should be (10,0) which is index 1
    int nearest = tree.find_nearest_unvisited(11.0, 0.0, visited);
    REQUIRE(nearest == 1);
  }

  SECTION("all visited returns -1") {
    std::vector<double> x = {0.0, 1.0};
    std::vector<double> y = {0.0, 0.0};

    BallTree tree(x, y);
    std::vector<bool> visited = {true, true};

    int nearest = tree.find_nearest_unvisited(0.5, 0.0, visited);
    REQUIRE(nearest == -1);
  }

  SECTION("larger dataset") {
    // Create a grid of 100 points
    std::vector<double> x, y;
    for (int i = 0; i < 10; i++) {
      for (int j = 0; j < 10; j++) {
        x.push_back(static_cast<double>(i));
        y.push_back(static_cast<double>(j));
      }
    }

    BallTree tree(x, y);
    std::vector<bool> visited(100, false);

    // Query from (4.5, 4.5) - should find one of the nearby points
    int nearest = tree.find_nearest_unvisited(4.5, 4.5, visited);
    REQUIRE(nearest >= 0);
    REQUIRE(nearest < 100);

    // The nearest point should be within the 4x4, 4x5, 5x4, or 5x5 region
    int row = nearest / 10;
    int col = nearest % 10;
    REQUIRE(row >= 4);
    REQUIRE(row <= 5);
    REQUIRE(col >= 4);
    REQUIRE(col <= 5);
  }

  SECTION("progressive search covers all points") {
    std::vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<double> y = {0.0, 0.0, 0.0, 0.0, 0.0};

    BallTree tree(x, y);
    std::vector<bool> visited(5, false);
    std::vector<int> order;

    // Starting from (0,0), progressively find nearest
    double qx = 0.0, qy = 0.0;
    for (int i = 0; i < 5; i++) {
      int nearest = tree.find_nearest_unvisited(qx, qy, visited);
      REQUIRE(nearest != -1);
      order.push_back(nearest);
      visited[nearest] = true;
      qx = x[nearest];
      qy = y[nearest];
    }

    // Should visit all 5 points
    REQUIRE(order.size() == 5);
    std::sort(order.begin(), order.end());
    for (int i = 0; i < 5; i++) {
      REQUIRE(order[i] == i);
    }
  }
}

TEST_CASE("BallTree with geographic coordinates", "[balltree]") {
  SECTION("European cities") {
    // Some European cities (lon, lat)
    std::vector<double> lon = {
      2.3522,    // Paris
      -0.1276,   // London
      13.4050,   // Berlin
      -3.7038,   // Madrid
      12.4964,   // Rome
      4.9041,    // Brussels
      4.8952,    // Amsterdam
      16.3738    // Vienna
    };
    std::vector<double> lat = {
      48.8566,   // Paris
      51.5074,   // London
      52.5200,   // Berlin
      40.4168,   // Madrid
      41.9028,   // Rome
      50.8503,   // Brussels
      52.3702,   // Amsterdam
      48.2082    // Vienna
    };

    BallTree tree(lon, lat);
    std::vector<bool> visited(8, false);

    // From Paris, nearest should be Brussels (index 5) at ~265 km
    int nearest = tree.find_nearest_unvisited(2.3522, 48.8566, visited);
    REQUIRE(nearest == 0);  // Paris itself

    visited[0] = true;  // Mark Paris visited
    nearest = tree.find_nearest_unvisited(2.3522, 48.8566, visited);
    // Brussels (5) or London (1) should be closest
    REQUIRE((nearest == 5 || nearest == 1));
  }
}

TEST_CASE("BallTree leaf size behavior", "[balltree]") {
  SECTION("exactly leaf size points") {
    std::vector<double> x, y;
    for (int i = 0; i < BallTree::LEAF_SIZE; i++) {
      x.push_back(static_cast<double>(i));
      y.push_back(0.0);
    }

    BallTree tree(x, y);
    REQUIRE(tree.n_points == BallTree::LEAF_SIZE);
    // Should create a single leaf node
    REQUIRE(tree.nodes.size() >= 1);
  }

  SECTION("just above leaf size") {
    std::vector<double> x, y;
    for (int i = 0; i < BallTree::LEAF_SIZE + 1; i++) {
      x.push_back(static_cast<double>(i));
      y.push_back(0.0);
    }

    BallTree tree(x, y);
    // Should split into multiple nodes
    REQUIRE(tree.nodes.size() > 1);
  }
}
