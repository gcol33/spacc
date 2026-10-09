// core/inext_core.h
// Size- and coverage-based rarefaction and extrapolation of Hill numbers for
// one abundance sample (Chao & Jost 2012; Chao et al. 2014), following the
// iNEXT.3D estimators: rarefied Hill numbers at integer sizes interpolated
// linearly in between, the observed sample at its own size, and the
// asymptote-anchored extrapolation beyond it. Pure C++, no Rcpp.

#ifndef SPACC_CORE_INEXT_H
#define SPACC_CORE_INEXT_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace spacc {

// Digamma by upward recurrence to x >= 10 and the asymptotic series there.
inline double digamma(double x) {
  double result = 0.0;
  while (x < 10.0) {
    result -= 1.0 / x;
    x += 1.0;
  }
  const double f = 1.0 / (x * x);
  return result + std::log(x) - 0.5 / x -
         f * (1.0 / 12 - f * (1.0 / 120 - f * (1.0 / 252 - f * (1.0 / 240 - f / 132))));
}

inline double lchoose(double n, double k) {
  return std::lgamma(n + 1.0) - std::lgamma(k + 1.0) - std::lgamma(n - k + 1.0);
}

// An abundance sample reduced to its positive counts.
struct AbundanceSample {
  std::vector<double> x;   // positive abundances
  double n = 0.0;          // sample size, sum of x
  int f1 = 0, f2 = 0;      // singletons, doubletons

  AbundanceSample() = default;

  template<class T>
  explicit AbundanceSample(const std::vector<T>& abundances) {
    for (const T& a : abundances) add(static_cast<double>(a));
  }

  void add(double a) {
    if (a <= 0) return;
    x.push_back(a);
    n += a;
    if (a == 1.0) f1++;
    if (a == 2.0) f2++;
  }

  // A = n f0 / (n f0 + f1) with the Chao1 undetected count f0: the factor by
  // which the coverage deficit shrinks per added individual beyond n. Equals
  // (n-1) f1 / ((n-1) f1 + 2 f2) when f2 > 0; 1 without singletons.
  double deficit_ratio() const {
    if (f1 == 0) return 1.0;
    const double f0 = f2 > 0 ? (n - 1) / n * f1 * f1 / (2.0 * f2)
                             : (n - 1) / n * f1 * (f1 - 1) / 2.0;
    return n * f0 / (n * f0 + f1);
  }
};

// Hill number of order q of the sample's relative abundances.
inline double observed_hill(const AbundanceSample& s, double q) {
  if (s.x.empty()) return 0.0;
  if (q == 0.0) return static_cast<double>(s.x.size());
  if (q == 1.0) {
    double h = 0.0;
    for (double a : s.x) { const double p = a / s.n; h -= p * std::log(p); }
    return std::exp(h);
  }
  double sum = 0.0;
  for (double a : s.x) sum += std::pow(a / s.n, q);
  return std::pow(sum, 1.0 / (1.0 - q));
}

// Expected Hill number of a random subsample of m < n individuals drawn
// without replacement: the species counts Y_i are hypergeometric, and
//   qD(m) = (sum_i E[(Y_i / m)^q 1{Y_i > 0}])^(1/(1-q)),
// exp(-sum_i E[(Y_i/m) log(Y_i/m)]) at q = 1. q = 0 and q = 2 use the closed
// forms of P(Y_i > 0) and E[Y_i^2]; other orders sum the pmf over the k where
// it is not negligible (log pmf within 40 of its mode).
inline double rarefied_hill(const AbundanceSample& s, double m, double q) {
  const double n = s.n;
  if (m <= 0) return 0.0;
  if (m >= n) return observed_hill(s, q);

  if (q == 0.0) {
    double S = 0.0;
    for (double a : s.x) {
      S += (n - a < m) ? 1.0 : 1.0 - std::exp(lchoose(n - a, m) - lchoose(n, m));
    }
    return S;
  }
  if (q == 2.0) {
    double sum = 0.0;
    for (double a : s.x) {
      const double mean = m * a / n;
      const double var = mean * (1.0 - a / n) * (n - m) / (n - 1.0);
      sum += (var + mean * mean) / (m * m);
    }
    return 1.0 / sum;
  }

  const double log_total = lchoose(n, m);
  double sum = 0.0;
  for (double a : s.x) {
    const double k_lo = std::max(1.0, m - (n - a));
    const double k_hi = std::min(a, m);
    if (k_lo > k_hi) continue;
    auto log_pmf = [&](double k) { return lchoose(a, k) + lchoose(n - a, m - k) - log_total; };
    auto term = [&](double k) {
      const double p = k / m;
      return q == 1.0 ? -p * std::log(p) : std::pow(p, q);
    };
    const double mode = std::min(k_hi, std::max(k_lo, std::floor((m + 1) * (a + 1) / (n + 2))));
    const double peak = log_pmf(mode);
    sum += term(mode) * std::exp(peak);
    for (double k = mode + 1; k <= k_hi; k++) {
      const double lp = log_pmf(k);
      if (lp < peak - 40.0) break;
      sum += term(k) * std::exp(lp);
    }
    for (double k = mode - 1; k >= k_lo; k--) {
      const double lp = log_pmf(k);
      if (lp < peak - 40.0) break;
      sum += term(k) * std::exp(lp);
    }
  }
  return q == 1.0 ? std::exp(sum) : std::pow(sum, 1.0 / (1.0 - q));
}

// Second-order term of the Chao et al. (2013) Shannon entropy estimator.
inline double shannon_second_order(double n, double f1, double f2) {
  double A;
  if (f2 > 0) A = 2 * f2 / ((n - 1) * f1 + 2 * f2);
  else if (f1 != 0) A = 2 / ((n - 1) * (f1 - 1) + 2);
  else A = 1;
  if (A == 1 || f1 == 0) return 0.0;
  const double pow_term = std::pow(1 - A, -n + 1);
  if (!std::isfinite(pow_term)) return 0.0;
  double partial = 0.0;
  for (int r = 1; r < n; r++) partial += std::pow(1 - A, r) / r;
  return f1 / n * pow_term * (-std::log(A) - partial);
}

// Asymptotic Hill number estimators for integer q: Chao1 (q = 0), the Chao
// et al. (2013) Shannon estimator (q = 1) and the unbiased estimator of
// sum p^q for integer q >= 2. Never below the observed value. NaN for
// non-integer q.
inline double asymptotic_hill(const AbundanceSample& s, double q) {
  const double n = s.n;
  const double obs = observed_hill(s, q);
  if (s.x.empty()) return 0.0;
  double est;
  if (q == 0.0) {
    est = s.x.size() + (n - 1) / n *
          (s.f2 > 0 ? static_cast<double>(s.f1) * s.f1 / (2.0 * s.f2)
                    : s.f1 * (s.f1 - 1) / 2.0);
  } else if (q == 1.0) {
    double A = 0.0;
    const double dg_n = digamma(n);
    for (double a : s.x) A += a / n * (dg_n - digamma(a));
    est = std::exp(A + shannon_second_order(n, s.f1, s.f2));
  } else if (q == std::round(q)) {
    if (s.f1 == static_cast<int>(s.x.size()) && q == 2.0) return obs;
    double A = 0.0;
    for (double a : s.x) if (a >= q) A += std::exp(lchoose(a, q) - lchoose(n, q));
    est = std::pow(A, 1.0 / (1.0 - q));
  } else {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return std::max(est, obs);
}

// Hill number at sample size m: rarefied below n (linear between the integer
// sizes around m), the observed value at n, and extrapolated above n,
//   qD(n + m*) = obs + (asy - obs) (1 - (1 - b)^m*),
// with b fixed by the rarefied value at n - 1; q = 2 uses its closed form.
inline double hill_at_size(const AbundanceSample& s, double m, double q) {
  const double n = s.n;
  if (m < n) {
    const double lo = std::floor(m), hi = std::ceil(m);
    if (lo == hi) return rarefied_hill(s, m, q);
    return (hi - m) * rarefied_hill(s, lo, q) + (m - lo) * rarefied_hill(s, hi, q);
  }
  const double obs = observed_hill(s, q);
  if (m == n) return obs;

  const double asy = asymptotic_hill(s, q);
  if (std::isnan(asy)) return asy;
  const double r_prev = rarefied_hill(s, n - 1, q);
  const double b = (asy != obs) ? (obs - r_prev) / (asy - r_prev) : 0.0;
  if (q == 2.0) {
    if (b == 0.0) return asy;
    double lambda = 0.0;
    for (double a : s.x) lambda += a * (a - 1) / (n * (n - 1));
    return 1.0 / (1.0 / m + (1.0 - 1.0 / m) * lambda);
  }
  return obs + (asy - obs) * (1.0 - std::pow(1.0 - b, m - n));
}

// Expected sample coverage at size m (Chao & Jost 2012): rarefied below n,
// 1 - f1/n A at n, and 1 - f1/n A^(m - n + 1) above it.
inline double coverage_at_size(const AbundanceSample& s, double m) {
  const double n = s.n;
  if (n <= 0) return 0.0;
  const double A = s.deficit_ratio();
  auto at_integer = [&](double k) {
    if (k >= n) return 1.0 - s.f1 / n * A;
    double deficit = 0.0;
    for (double a : s.x) {
      if (n - a < k) continue;
      deficit += a / n * std::exp(std::lgamma(n - a + 1) - std::lgamma(n - a - k + 1) -
                                  std::lgamma(n) + std::lgamma(n - k));
    }
    return 1.0 - deficit;
  };
  if (m < n) {
    const double lo = std::floor(m), hi = std::ceil(m);
    if (lo == hi) return at_integer(m);
    return (hi - m) * at_integer(lo) + (m - lo) * at_integer(hi);
  }
  if (m == n) return 1.0 - s.f1 / n * A;
  return 1.0 - s.f1 / n * std::pow(A, m - n + 1);
}

// Sample size at which the expected coverage equals C: bisection on the
// increasing rarefied coverage below the observed coverage, the closed-form
// inverse of the extrapolated coverage above it (infinite at C = 1). A sample
// whose coverage cannot grow (A = 0) stays at n.
inline double size_at_coverage(const AbundanceSample& s, double C) {
  const double n = s.n;
  const double ref = coverage_at_size(s, n);
  if (C <= 0.0) return 0.0;
  if (C == ref) return n;
  if (C < ref) {
    double lo = 0.0, hi = n;
    for (int it = 0; it < 200 && hi - lo > 1e-10 * n; it++) {
      const double mid = 0.5 * (lo + hi);
      if (coverage_at_size(s, mid) < C) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
  }
  const double A = s.deficit_ratio();
  if (A <= 0.0 || A >= 1.0) return n;
  if (C >= 1.0) return std::numeric_limits<double>::infinity();
  return n + (std::log(n / s.f1) + std::log(1 - C)) / std::log(A) - 1;
}

inline double hill_at_coverage(const AbundanceSample& s, double C, double q) {
  return hill_at_size(s, size_at_coverage(s, C), q);
}

} // namespace spacc

#endif // SPACC_CORE_INEXT_H
