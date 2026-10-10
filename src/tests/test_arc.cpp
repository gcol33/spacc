// test_arc.cpp
// Catch2 tests for the accumulation rate curve (ARC) estimator

#include "catch.hpp"
#include "../core/arc_core.h"
#include <random>

using namespace spacc::arc;

// Simpson integral of f(n) over [n0, inf) on u = log n.
template<typename F>
static double tail_numeric(F f, double n0, double span) {
  const int steps = 200000;
  const double a = std::log(n0), h = span / steps;
  double total = 0.0;
  for (int i = 0; i <= steps; i++) {
    const double u = a + i * h;
    const double wgt = (i == 0 || i == steps) ? 1.0 : (i % 2 ? 4.0 : 2.0);
    total += wgt * f(std::exp(u)) * std::exp(u);
  }
  return total * h / 3.0;
}

TEST_CASE("special functions", "[arc]") {
  REQUIRE(norm_cdf(0.0) == Approx(0.5));
  REQUIRE(norm_quantile(0.975) == Approx(1.959964).epsilon(1e-6));
  // I_0.3(2, 3) = sum_{j=2}^{4} C(4, j) 0.3^j 0.7^(4-j)
  REQUIRE(incomplete_beta(0.3, 2.0, 3.0) == Approx(0.3483).epsilon(1e-10));
  REQUIRE(incomplete_beta(0.9, 2.0, 3.0) == Approx(1.0 - incomplete_beta(0.1, 3.0, 2.0)));
  // Gamma(1, x) = exp(-x), Gamma(2, x) = (x + 1) exp(-x); both branches
  for (double x : {0.5, 3.0, 20.0}) {
    REQUIRE(log_upper_gamma(1.0, x) == Approx(-x).epsilon(1e-12));
    REQUIRE(log_upper_gamma(2.0, x) == Approx(std::log(x + 1.0) - x).epsilon(1e-12));
  }
}

TEST_CASE("curve tails match numerical integrals", "[arc]") {
  const double n0 = 500.0;
  SECTION("log-normal") {
    const double th[2] = {5.0, std::log(1.5)};
    auto f = [&](double n) { return curve_eval(LOGNORMAL, std::log(n), th, nullptr); };
    REQUIRE(curve_tail(LOGNORMAL, th, n0) == Approx(tail_numeric(f, n0, 40.0)).epsilon(1e-6));
  }
  SECTION("log-logistic") {
    const double th[2] = {std::log(1e-4), std::log(2.5)};
    auto f = [&](double n) { return curve_eval(LOGLOGISTIC, std::log(n), th, nullptr); };
    REQUIRE(curve_tail(LOGLOGISTIC, th, n0) == Approx(tail_numeric(f, n0, 60.0)).epsilon(1e-6));
    const double flat[2] = {std::log(1e-4), std::log(0.9)};
    REQUIRE(std::isinf(curve_tail(LOGLOGISTIC, flat, n0)));
  }
  SECTION("Weibull") {
    const double th[2] = {std::log(0.05), std::log(0.6)};
    auto f = [&](double n) { return curve_eval(WEIBULL, std::log(n), th, nullptr); };
    REQUIRE(curve_tail(WEIBULL, th, n0) == Approx(tail_numeric(f, n0, 40.0)).epsilon(1e-6));
  }
}

TEST_CASE("curve gradients match finite differences", "[arc]") {
  const double th0[3][2] = {{5.0, 0.3}, {-8.0, 0.7}, {-3.0, -0.5}};
  for (int c = 0; c < 3; c++) {
    const double ln = std::log(300.0);
    double grad[2];
    curve_eval(c, ln, th0[c], grad);
    for (int k = 0; k < 2; k++) {
      double up[2] = {th0[c][0], th0[c][1]}, dn[2] = {th0[c][0], th0[c][1]};
      up[k] += 1e-6; dn[k] -= 1e-6;
      const double fd = (curve_eval(c, ln, up, nullptr) - curve_eval(c, ln, dn, nullptr)) / 2e-6;
      REQUIRE(grad[k] == Approx(fd).margin(1e-8));
    }
  }
}

TEST_CASE("rate curve differences the expected accumulation", "[arc]") {
  // y_1 = 4, y_2 = 2, y_3 = 1: N = 11 tags, unthinned points n = 2..11
  const std::vector<double> yk = {4.0, 2.0, 1.0};
  const RateCurve rc = rate_curve(yk, 1000);
  REQUIRE(rc.n.size() == 10);
  REQUIRE(rc.n.front() == 2.0);
  REQUIRE(rc.n.back() == 11.0);
  auto S = [&](double n) {
    double t = 0.0;
    for (int k = 0; k < 3; k++) t += yk[k] * (1.0 - std::pow(1.0 - n / 11.0, k + 1));
    return t;
  };
  double sum = 0.0;
  for (double v : rc.s) sum += v;
  REQUIRE(sum == Approx(S(11.0) - S(1.0)));
  REQUIRE(S(11.0) == Approx(7.0));

  const RateCurve thin = rate_curve(std::vector<double>{3000.0, 500.0}, 100);
  REQUIRE(thin.n.size() == 100);
  REQUIRE(thin.n.back() == 4000.0);
}

TEST_CASE("weights are the minimum over the capped box", "[arc]") {
  std::mt19937 rng(7);
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  const int m = 40;
  std::vector<double> s(m), omega(m, 1.0);
  std::vector<std::vector<double>> f(3, std::vector<double>(m));
  for (int i = 0; i < m; i++) {
    const double n = 10.0 + 5.0 * i;
    f[0][i] = std::exp(-0.02 * n);
    f[1][i] = 1.0 / (1.0 + 0.01 * n);
    f[2][i] = std::exp(-0.05 * n);
    s[i] = 0.6 * f[0][i] + 0.2 * f[1][i] + 0.1 * f[2][i] + 0.01 * (unif(rng) - 0.5);
  }
  const double cap[3] = {0.8, 0.3, 0.4};
  double w[3];
  optimise_weights(f, s.data(), omega.data(), m, cap, w);

  auto loss = [&](const double* v) {
    double total = 0.0;
    for (int i = 0; i < m; i++) {
      const double r = s[i] - (v[0] * f[0][i] + v[1] * f[1][i] + v[2] * f[2][i]);
      total += omega[i] * (r * r + v[1] * v[2] / s[i]);
    }
    return total;
  };
  for (int a = 0; a < 3; a++) {
    REQUIRE(w[a] >= 0.0);
    REQUIRE(w[a] <= cap[a]);
  }
  const double best = loss(w);
  for (int trial = 0; trial < 20000; trial++) {
    const double v[3] = {cap[0] * unif(rng), cap[1] * unif(rng), cap[2] * unif(rng)};
    REQUIRE(loss(v) >= best - 1e-12);
  }
}

TEST_CASE("ridge NNLS satisfies the KKT conditions", "[arc]") {
  std::mt19937 rng(11);
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  const int R = 12, J = 8;
  std::vector<double> A(R * J), y(R);
  for (double& a : A) a = unif(rng);
  for (double& v : y) v = 5.0 * unif(rng) - 1.0;
  std::vector<double> H(J * J, 0.0), g(J, 0.0);
  for (int j = 0; j < J; j++) {
    for (int r = 0; r < R; r++) g[j] += A[r * J + j] * y[r];
    for (int l = 0; l < J; l++)
      for (int r = 0; r < R; r++) H[j * J + l] += A[r * J + j] * A[r * J + l];
  }
  for (double u : {0.0, 0.5}) {
    std::vector<double> v;
    ridge_nnls(H, g, J, u, v);
    for (int j = 0; j < J; j++) {
      double grad = g[j] - u * v[j];
      for (int l = 0; l < J; l++) grad -= H[j * J + l] * v[l];
      REQUIRE(v[j] >= 0.0);
      if (v[j] > 0.0) REQUIRE(grad == Approx(0.0).margin(1e-9));
      else REQUIRE(grad <= 1e-9);
    }
  }
}

TEST_CASE("smoothing keeps the number of tags", "[arc]") {
  const std::vector<double> y = {50.0, 20.0, 9.0, 5.0, 0.0, 2.0, 7.0};
  std::vector<double> out;
  REQUIRE(smooth_counts(y, 0.5, out));
  REQUIRE(count_tags(out) == Approx(count_tags(y)));
  REQUIRE(out[4] == 0.0);
  std::vector<double> single;
  REQUIRE_FALSE(smooth_counts(std::vector<double>{3.0, 0.0, 1.0}, 0.5, single));
}

TEST_CASE("estimate returns S_obs plus a positive number of undetected species", "[arc]") {
  // frequency counts of a sparse Poisson-Gamma sample (y_1..y_20, y_21+)
  const std::vector<double> y = {264, 133, 53, 27, 18, 10, 6, 3, 0, 1, 0, 1,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0};
  const Estimate e = estimate(y, 516.0, false, 0.5, 1000);
  REQUIRE(e.n_portions == 7);
  REQUIRE(e.estimate > 516.0);
  REQUIRE(e.S[WEIBULL] < e.S[LOGNORMAL]);
  REQUIRE(e.S[LOGNORMAL] < e.S[LOGLOGISTIC]);
}
