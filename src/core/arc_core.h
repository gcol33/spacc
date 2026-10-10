// core/arc_core.h
// Accumulation rate curve (ARC) estimator of total richness (Shestopaloff
// 2025, Environmental and Ecological Statistics 32: 311-345).
// Pure C++ - testable with Catch2.
//
// Data are frequency counts y = (y_1, ..., y_K, y_K+): y_k species seen k
// times, y_K+ species seen more than K times. Only species seen at most K
// times enter the rate curve.

#ifndef SPACC_CORE_ARC_H
#define SPACC_CORE_ARC_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace spacc {
namespace arc {

constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kPi = 3.14159265358979323846;

// ---- special functions --------------------------------------------------

inline double norm_cdf(double z) { return 0.5 * std::erfc(-z / std::sqrt(2.0)); }

inline double norm_pdf(double z) { return std::exp(-0.5 * z * z) / std::sqrt(2.0 * kPi); }

// Inverse of norm_cdf by bisection; only used for starting values.
inline double norm_quantile(double p) {
  double lo = -40.0, hi = 40.0;
  for (int i = 0; i < 200; i++) {
    const double mid = 0.5 * (lo + hi);
    if (norm_cdf(mid) < p) lo = mid; else hi = mid;
  }
  return 0.5 * (lo + hi);
}

// Modified Lentz evaluation of the continued fraction of the regularised
// incomplete beta function (Press et al., Numerical Recipes, 6.4).
inline double beta_cf(double a, double b, double x) {
  const double tiny = 1e-300, eps = 1e-15;
  double c = 1.0, d = 1.0 - (a + b) * x / (a + 1.0);
  if (std::fabs(d) < tiny) d = tiny;
  d = 1.0 / d;
  double h = d;
  for (int m = 1; m <= 10000; m++) {
    const int m2 = 2 * m;
    double aa = m * (b - m) * x / ((a + m2 - 1.0) * (a + m2));
    d = 1.0 + aa * d; if (std::fabs(d) < tiny) d = tiny;
    c = 1.0 + aa / c; if (std::fabs(c) < tiny) c = tiny;
    d = 1.0 / d; h *= d * c;
    aa = -(a + m) * (a + b + m) * x / ((a + m2) * (a + m2 + 1.0));
    d = 1.0 + aa * d; if (std::fabs(d) < tiny) d = tiny;
    c = 1.0 + aa / c; if (std::fabs(c) < tiny) c = tiny;
    d = 1.0 / d;
    const double del = d * c;
    h *= del;
    if (std::fabs(del - 1.0) < eps) break;
  }
  return h;
}

// Regularised incomplete beta function I_x(a, b).
inline double incomplete_beta(double x, double a, double b) {
  if (x <= 0.0) return 0.0;
  if (x >= 1.0) return 1.0;
  const double log_front = std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b) +
                           a * std::log(x) + b * std::log1p(-x);
  if (x < (a + 1.0) / (a + b + 2.0)) return std::exp(log_front) * beta_cf(a, b, x) / a;
  return 1.0 - std::exp(log_front) * beta_cf(b, a, 1.0 - x) / b;
}

// log of the upper incomplete gamma function Gamma(a, x) = int_x^inf t^(a-1) e^-t dt.
inline double log_upper_gamma(double a, double x) {
  if (x <= 0.0) return std::lgamma(a);
  if (x < a + 1.0) {
    // series for the lower function, Gamma(a, x) = Gamma(a) (1 - P(a, x))
    double ap = a, sum = 1.0 / a, del = sum;
    for (int n = 0; n < 100000; n++) {
      ap += 1.0; del *= x / ap; sum += del;
      if (std::fabs(del) < std::fabs(sum) * 1e-16) break;
    }
    const double P = std::exp(-x + a * std::log(x) - std::lgamma(a)) * sum;
    return std::lgamma(a) + std::log1p(-std::min(P, 1.0));
  }
  // continued fraction (Lentz) for Gamma(a, x) e^x x^-a
  const double tiny = 1e-300;
  double b = x + 1.0 - a, c = 1.0 / tiny, d = 1.0 / b, h = d;
  for (int i = 1; i < 100000; i++) {
    const double an = -i * (i - a);
    b += 2.0;
    d = an * d + b; if (std::fabs(d) < tiny) d = tiny;
    c = b + an / c; if (std::fabs(c) < tiny) c = tiny;
    d = 1.0 / d;
    const double del = d * c;
    h *= del;
    if (std::fabs(del - 1.0) < 1e-16) break;
  }
  return -x + a * std::log(x) + std::log(h);
}

// ---- the three decay curves (Eq. 4) --------------------------------------
// Parameters are kept on an unconstrained scale:
//   log-normal    f = 1 - Phi((log n - mu) / sigma),  th = (mu, log sigma)
//   log-logistic  f = 1 / (1 + alpha n^gamma),        th = (log alpha, log gamma)
//   Weibull       f = exp(-lambda n^gamma),           th = (log lambda, log gamma)

enum Curve { LOGNORMAL = 0, LOGLOGISTIC = 1, WEIBULL = 2 };

// Value of curve c at n = exp(ln), and its gradient in th when grad != nullptr.
inline double curve_eval(int c, double ln, const double* th, double* grad) {
  if (c == LOGNORMAL) {
    const double sigma = std::exp(th[1]);
    const double z = (ln - th[0]) / sigma;
    if (grad) {
      const double dens = norm_pdf(z);
      grad[0] = dens / sigma;
      grad[1] = dens * z;
    }
    return norm_cdf(-z);
  }
  const double gamma = std::exp(th[1]);
  const double log_g = std::min(th[0] + gamma * ln, 700.0);
  const double g = std::exp(log_g);
  if (c == LOGLOGISTIC) {
    const double f = 1.0 / (1.0 + g);
    if (grad) {
      const double dfdlogg = -g * f * f;
      grad[0] = dfdlogg;
      grad[1] = dfdlogg * gamma * ln;
    }
    return f;
  }
  const double f = std::exp(-g);
  if (grad) {
    grad[0] = -g * f;
    grad[1] = -g * f * gamma * ln;
  }
  return f;
}

// Integral of curve c from n0 to infinity: the undetected species it implies.
// The log-logistic integral is finite only for gamma > 1.
inline double curve_tail(int c, const double* th, double n0) {
  const double ln0 = std::log(n0);
  if (c == LOGNORMAL) {
    // E[(X - n0)^+] for X log-normal(mu, sigma)
    const double mu = th[0], sigma = std::exp(th[1]);
    return std::exp(mu + 0.5 * sigma * sigma) * norm_cdf((mu + sigma * sigma - ln0) / sigma) -
           n0 * norm_cdf((mu - ln0) / sigma);
  }
  const double gamma = std::exp(th[1]), a = 1.0 / gamma;
  if (c == WEIBULL) {
    // (1/gamma) lambda^(-1/gamma) Gamma(1/gamma, lambda n0^gamma)
    const double x = std::exp(th[0] + gamma * ln0);
    return std::exp(-std::log(gamma) - a * th[0] + log_upper_gamma(a, x));
  }
  if (gamma <= 1.0) return kInf;
  // alpha^(-1/gamma) / gamma * B(1 - a, a) * I_t0(1 - a, a), t0 = 1 / (1 + alpha n0^gamma)
  const double log_g = th[0] + gamma * ln0;
  const double t0 = 1.0 / (1.0 + std::exp(std::min(log_g, 700.0)));
  return std::exp(-a * th[0] - std::log(gamma)) * kPi / std::sin(kPi * a) *
         incomplete_beta(t0, 1.0 - a, a);
}

// ---- rate curve (Eq. 2-3) ------------------------------------------------

struct RateCurve {
  std::vector<double> n;  // number of tags drawn
  std::vector<double> s;  // expected new species per tag, S(n) - S(n - 1)
};

// Expected species after n of the N tags of species seen at most K times,
// S(n) = sum_k y_k (1 - (1 - n / N)^k), differenced at n = 2, ..., N and
// thinned to at most max_points evenly spaced points. yk holds y_1..y_K.
inline RateCurve rate_curve(const std::vector<double>& yk, int max_points) {
  double tags = 0.0;
  for (std::size_t k = 0; k < yk.size(); k++) tags += (k + 1) * yk[k];
  const double N = std::round(tags);
  RateCurve rc;
  if (N < 3.0) return rc;

  auto S = [&](double n) {
    double total = 0.0;
    const double base = 1.0 - n / N;
    for (std::size_t k = 0; k < yk.size(); k++) {
      if (yk[k] != 0.0) total += yk[k] * (1.0 - std::pow(base, static_cast<double>(k + 1)));
    }
    return total;
  };

  const long L = static_cast<long>(N) - 1;  // points n = 2..N
  const long m = std::min<long>(L, max_points);
  long last = -1;
  for (long j = 0; j < m; j++) {
    const long i = m == 1 ? 0 : static_cast<long>(std::nearbyint(
      static_cast<double>(L - 1) * j / (m - 1)));
    if (i == last) continue;
    last = i;
    const double n = 2.0 + i;
    rc.n.push_back(n);
    rc.s.push_back(S(n) - S(n - 1.0));
  }
  return rc;
}

// ---- weighted least-squares fit of one curve (Eq. 5) ----------------------

struct CurveFit {
  double th[2];
  double sse;
  bool ok;
};

// The Weibull is fitted on the log scale, log s against -lambda n^gamma,
// where it is linear in lambda; the other two on the scale of s.
inline double fit_scale(int c, double ln, const double* th, double* grad) {
  if (c != WEIBULL) return curve_eval(c, ln, th, grad);
  const double gamma = std::exp(th[1]);
  const double g = std::exp(std::min(th[0] + gamma * ln, 700.0));
  if (grad) {
    grad[0] = -g;
    grad[1] = -g * gamma * ln;
  }
  return -g;
}

// Levenberg-Marquardt on sum_i w_i (t_i - fit_scale(n_i))^2 from start th0,
// with t = s, or log s for the Weibull; ln holds log n_i.
inline CurveFit lm_fit(int c, const double* ln, const double* t, const double* w, int m,
                       const double* th0) {
  CurveFit fit;
  fit.th[0] = th0[0]; fit.th[1] = th0[1];
  fit.ok = false;

  auto sse_at = [&](const double* th) {
    double total = 0.0;
    for (int i = 0; i < m; i++) {
      const double r = t[i] - fit_scale(c, ln[i], th, nullptr);
      total += w[i] * r * r;
    }
    return std::isfinite(total) ? total : kInf;
  };

  double sse = sse_at(fit.th);
  if (!std::isfinite(sse)) { fit.sse = kInf; return fit; }
  double damping = 1e-3;
  for (int iter = 0; iter < 300; iter++) {
    double JtJ[3] = {0.0, 0.0, 0.0}, Jtr[2] = {0.0, 0.0}, grad[2];
    for (int i = 0; i < m; i++) {
      const double r = t[i] - fit_scale(c, ln[i], fit.th, grad);
      JtJ[0] += w[i] * grad[0] * grad[0];
      JtJ[1] += w[i] * grad[0] * grad[1];
      JtJ[2] += w[i] * grad[1] * grad[1];
      Jtr[0] += w[i] * grad[0] * r;
      Jtr[1] += w[i] * grad[1] * r;
    }
    bool accepted = false;
    while (damping < 1e12) {
      const double a = JtJ[0] * (1.0 + damping), b = JtJ[1], d = JtJ[2] * (1.0 + damping);
      const double det = a * d - b * b;
      if (!(det > 0.0)) { damping *= 10.0; continue; }
      const double trial[2] = {fit.th[0] + (d * Jtr[0] - b * Jtr[1]) / det,
                               fit.th[1] + (a * Jtr[1] - b * Jtr[0]) / det};
      const double trial_sse = sse_at(trial);
      if (trial_sse < sse) {
        const double drop = (sse - trial_sse) / sse;
        fit.th[0] = trial[0]; fit.th[1] = trial[1];
        sse = trial_sse;
        damping = std::max(damping / 10.0, 1e-12);
        accepted = true;
        if (drop < 1e-10) iter = 300;
        break;
      }
      damping *= 10.0;
    }
    if (!accepted) break;
  }
  fit.sse = sse;
  fit.ok = std::isfinite(sse);
  return fit;
}

// Best of several starts, each placing the curve through the point
// (n_mid, s_mid) of the fitted range; ln holds log n.
inline CurveFit fit_curve(int c, const double* ln, const double* s, const double* w, int m) {
  std::vector<double> log_s;
  const double* t = s;
  if (c == WEIBULL) {
    log_s.resize(m);
    for (int i = 0; i < m; i++) log_s[i] = std::log(s[i]);
    t = log_s.data();
  }
  const int mid = m / 2;
  const double s_mid = std::min(std::max(s[mid], 1e-12), 1.0 - 1e-12);
  const double ln_mid = ln[mid];

  std::vector<double> shapes;
  if (c == LOGNORMAL) shapes = {0.25, 0.5, 1.0, 2.0, 4.0, 8.0};
  else if (c == LOGLOGISTIC) shapes = {0.5, 1.0, 1.5, 2.0, 3.0, 5.0};
  else shapes = {0.1, 0.25, 0.5, 1.0, 2.0};

  // starting values, ranked by their loss; the two best are refined
  std::vector<CurveFit> starts;
  for (double shape : shapes) {
    double th0[2];
    if (c == LOGNORMAL) {
      th0[0] = ln_mid - shape * norm_quantile(1.0 - s_mid);
      th0[1] = std::log(shape);
    } else if (c == LOGLOGISTIC) {
      th0[0] = std::log(1.0 / s_mid - 1.0) - shape * ln_mid;
      th0[1] = std::log(shape);
    } else {
      th0[0] = std::log(-std::log(s_mid)) - shape * ln_mid;
      th0[1] = std::log(shape);
    }
    CurveFit start;
    start.th[0] = th0[0]; start.th[1] = th0[1];
    double sse = 0.0;
    for (int i = 0; i < m; i++) {
      const double r = t[i] - fit_scale(c, ln[i], th0, nullptr);
      sse += w[i] * r * r;
    }
    start.sse = std::isfinite(sse) ? sse : kInf;
    starts.push_back(start);
  }
  std::sort(starts.begin(), starts.end(),
            [](const CurveFit& a, const CurveFit& b) { return a.sse < b.sse; });

  CurveFit best;
  best.ok = false; best.sse = kInf; best.th[0] = best.th[1] = kNaN;
  for (std::size_t i = 0; i < std::min<std::size_t>(2, starts.size()); i++) {
    const CurveFit f = lm_fit(c, ln, t, w, m, starts[i].th);
    if (f.ok && f.sse < best.sse) best = f;
  }
  return best;
}

// ---- weighted estimator (Eq. 6-8) ----------------------------------------

struct Combination {
  double unobs;  // weighted undetected species
  double S[3];   // undetected species per curve
  double w[3];   // weights
  bool ok;
};

// Upper bounds on the weights from the three curve estimates and S_obs (Eq. 7).
inline void weight_bounds(const double* S, double s_obs, double* cap) {
  const double S_ln = S[LOGNORMAL], S_ll = S[LOGLOGISTIC], S_wb = S[WEIBULL];
  auto ratios = [&](double s_new, double& r_wb, double& r_ll) {
    if (std::isinf(s_new)) { r_wb = 1.0; r_ll = 0.0; return; }
    r_wb = s_new / (s_obs + s_new);
    r_ll = s_obs / (s_obs + s_new);
  };
  double r_wb, r_ll;
  ratios(S_ln, r_wb, r_ll);
  const bool have_ll = std::isfinite(S_ll) || std::isinf(S_ll);
  const bool have_wb = std::isfinite(S_wb);

  // midpoint with the log-logistic when sparse, with the Weibull otherwise
  double S_star = S_ln;
  if (r_wb > r_ll) { if (have_ll) S_star = 0.5 * (S_ln + S_ll); }
  else if (have_wb) S_star = 0.5 * (S_ln + S_wb);
  double r_wb_star, r_ll_star;
  ratios(S_star, r_wb_star, r_ll_star);

  const double c_wb = have_wb ? S_wb / S_ln * r_wb : 0.0;
  const double c_ll = (have_ll && std::isfinite(S_ll)) ? std::min(S_ln / S_ll * r_ll, c_wb) : 0.0;
  double c_ln = std::max(0.5 * (r_wb + r_wb_star), 0.5 * (r_ll + r_ll_star));
  c_ln = std::max(c_ln, 1.0 - c_ll - c_wb + 1e-6);

  cap[LOGNORMAL] = c_ln; cap[LOGLOGISTIC] = c_ll; cap[WEIBULL] = c_wb;
  for (int i = 0; i < 3; i++) cap[i] = std::min(std::max(cap[i], 0.0), 1.0);
}

// The weight loss sum_n w_n [(s_n - sum_c w_c f_c(n))^2 + w_ll w_wb / s_n] as
// the quadratic 1/2 v'Hv - 2 b'v + const, H = 2G + P (e_ll e_wb' + e_wb e_ll').
struct WeightLoss {
  double H[3][3];
  double b[3];
};

inline WeightLoss weight_loss(const std::vector<std::vector<double>>& f, const double* s,
                              const double* omega, int m) {
  WeightLoss q = {{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {0, 0, 0}};
  double P = 0.0;
  for (int i = 0; i < m; i++) {
    for (int a = 0; a < 3; a++) {
      q.b[a] += omega[i] * f[a][i] * s[i];
      for (int c = 0; c < 3; c++) q.H[a][c] += 2.0 * omega[i] * f[a][i] * f[c][i];
    }
    P += omega[i] / s[i];
  }
  q.H[LOGLOGISTIC][WEIBULL] += P;
  q.H[WEIBULL][LOGLOGISTIC] += P;
  return q;
}

inline void weight_gradient(const WeightLoss& q, const double* v, double* g) {
  for (int a = 0; a < 3; a++) {
    g[a] = -2.0 * q.b[a];
    for (int c = 0; c < 3; c++) g[a] += q.H[a][c] * v[c];
  }
}

// Eigendecomposition A = U diag(lam) U' of the leading n x n block (n <= 3)
// of a symmetric matrix by cyclic Jacobi rotations. A is overwritten.
inline void symmetric_eigen(double A[3][3], int n, double* lam, double U[3][3]) {
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) U[i][j] = i == j ? 1.0 : 0.0;
  for (int sweep = 0; sweep < 50; sweep++) {
    double off = 0.0, total = 0.0;
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) {
        total += A[i][j] * A[i][j];
        if (i != j) off += A[i][j] * A[i][j];
      }
    if (off <= 1e-32 * total) break;
    for (int p = 0; p < n; p++) {
      for (int r = p + 1; r < n; r++) {
        if (A[p][r] == 0.0) continue;
        const double theta = (A[r][r] - A[p][p]) / (2.0 * A[p][r]);
        const double t = (theta >= 0.0 ? 1.0 : -1.0) / (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
        const double c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
        for (int k = 0; k < n; k++) {
          const double akp = A[k][p], akr = A[k][r];
          A[k][p] = c * akp - s * akr;
          A[k][r] = s * akp + c * akr;
        }
        for (int k = 0; k < n; k++) {
          const double apk = A[p][k], ark = A[r][k];
          A[p][k] = c * apk - s * ark;
          A[r][k] = s * apk + c * ark;
        }
        for (int k = 0; k < n; k++) {
          const double ukp = U[k][p], ukr = U[k][r];
          U[k][p] = c * ukp - s * ukr;
          U[k][r] = s * ukp + c * ukr;
        }
      }
    }
  }
  for (int i = 0; i < n; i++) lam[i] = A[i][i];
}

enum Bound { AT_ZERO, AT_CAP, FREE };

// Weights within [0, cap] (Eq. 7-8). The product penalty makes the loss
// indefinite, with one local minimum favouring the log-logistic and one the
// Weibull; the weights are the end point of the projected gradient flow
// dv/dt = -grad L(v) started at v = 0. On a face of the box (each weight at
// 0, at its cap, or free) the flow is linear, v_F' = -(H_FF v_F - r), and is
// solved in closed form in the eigenbasis of H_FF. A free weight reaching a
// bound is fixed there; a fixed weight whose gradient turns into the box is
// freed. The flow ends where no event follows.
inline void optimise_weights(const std::vector<std::vector<double>>& f, const double* s,
                             const double* omega, int m, const double* cap, double* w) {
  const WeightLoss q = weight_loss(f, s, omega, m);
  double b_scale = 0.0;
  for (int a = 0; a < 3; a++) b_scale = std::max(b_scale, std::fabs(q.b[a]));
  w[0] = w[1] = w[2] = 0.0;
  if (!(b_scale > 0.0)) return;

  double v[3] = {0.0, 0.0, 0.0};
  int state[3];
  for (int a = 0; a < 3; a++) state[a] = (cap[a] > 0.0 && q.b[a] > 0.0) ? FREE : AT_ZERO;
  const double tol = 1e-12;

  for (int face = 0; face < 100; face++) {
    int F[3], n = 0;
    for (int a = 0; a < 3; a++) if (state[a] == FREE) F[n++] = a;
    if (n == 0) break;

    double A[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, U[3][3], lam[3], rho[3], z0[3];
    double r[3];
    for (int i = 0; i < n; i++) {
      r[i] = 2.0 * q.b[F[i]];
      for (int c = 0; c < 3; c++) if (state[c] != FREE) r[i] -= q.H[F[i]][c] * v[c];
      for (int j = 0; j < n; j++) A[i][j] = q.H[F[i]][F[j]];
    }
    symmetric_eigen(A, n, lam, U);
    double lam_max = 0.0;
    for (int k = 0; k < n; k++) {
      rho[k] = 0.0; z0[k] = 0.0;
      for (int i = 0; i < n; i++) { rho[k] += U[i][k] * r[i]; z0[k] += U[i][k] * v[F[i]]; }
      lam_max = std::max(lam_max, std::fabs(lam[k]));
    }
    if (!(lam_max > 0.0)) break;

    // position at time t on this face
    auto flow = [&](double t, double* out) {
      for (int a = 0; a < 3; a++) out[a] = v[a];
      double z[3];
      for (int k = 0; k < n; k++) {
        const double x = -lam[k] * t;
        const double decay = std::exp(x);
        const double gain = lam[k] == 0.0 ? t : -std::expm1(x) / lam[k];
        z[k] = z0[k] * decay + rho[k] * gain;
      }
      for (int i = 0; i < n; i++) {
        out[F[i]] = 0.0;
        for (int k = 0; k < n; k++) out[F[i]] += U[i][k] * z[k];
      }
    };
    // smallest constraint margin at time t: free weights stay in [0, cap],
    // fixed ones keep their gradient pointing out of the box
    auto margin = [&](double t, int& which) {
      double x[3], g[3];
      flow(t, x);
      weight_gradient(q, x, g);
      double worst = kInf;
      which = -1;
      for (int a = 0; a < 3; a++) {
        double mg;
        if (state[a] == FREE) mg = std::min(x[a], cap[a] - x[a]) / (cap[a] > 0.0 ? cap[a] : 1.0);
        else if (!(cap[a] > 0.0)) continue;
        else mg = (state[a] == AT_ZERO ? g[a] : -g[a]) / b_scale;
        if (!(mg >= worst)) { worst = mg; which = a; }
      }
      return worst;
    };

    const double t0 = 1e-6 / lam_max;
    const int per_decade = 20, n_grid = per_decade * 46;
    double lo = 0.0, hi = -1.0;
    int which = -1;
    for (int i = 0; i <= n_grid; i++) {
      const double t = t0 * std::pow(10.0, static_cast<double>(i) / per_decade);
      if (!(margin(t, which) >= -tol)) { hi = t; break; }
      lo = t;
    }
    if (hi < 0.0) {
      flow(lo, v);
      break;
    }
    for (int it = 0; it < 200 && hi - lo > 1e-15 * hi; it++) {
      const double mid = 0.5 * (lo + hi);
      if (!(margin(mid, which) >= -tol)) hi = mid; else lo = mid;
    }
    margin(hi, which);
    flow(hi, v);
    for (int a = 0; a < 3; a++) v[a] = std::min(std::max(v[a], 0.0), cap[a]);
    if (which < 0) break;
    if (state[which] == FREE) {
      state[which] = v[which] < 0.5 * cap[which] ? AT_ZERO : AT_CAP;
      v[which] = state[which] == AT_CAP ? cap[which] : 0.0;
    } else {
      state[which] = FREE;
    }
  }
  for (int a = 0; a < 3; a++) w[a] = std::min(std::max(v[a], 0.0), cap[a]);
}

// Fit the three curves to the part of the rate curve at or beyond quantile
// `portion` of n and combine them.
inline Combination fit_portion(const RateCurve& rc, double portion, double s_obs) {
  Combination out;
  out.ok = false; out.unobs = kNaN;
  for (int i = 0; i < 3; i++) { out.S[i] = kNaN; out.w[i] = 0.0; }

  // type-7 quantile of n
  const std::size_t L = rc.n.size();
  const double h = (L - 1) * portion;
  const std::size_t lo = static_cast<std::size_t>(std::floor(h));
  const double cut = lo + 1 < L ? rc.n[lo] + (h - lo) * (rc.n[lo + 1] - rc.n[lo]) : rc.n[lo];
  std::size_t first = 0;
  while (first < L && rc.n[first] < cut) first++;
  const int m = static_cast<int>(L - first);
  if (m < 3) return out;

  const double* s = rc.s.data() + first;
  std::vector<double> ln(m);
  for (int i = 0; i < m; i++) ln[i] = std::log(rc.n[first + i]);
  std::vector<double> omega(m);
  double total = 0.0;
  for (int i = 0; i < m; i++) { omega[i] = static_cast<double>(i + 1) * (i + 1); total += omega[i]; }
  for (int i = 0; i < m; i++) omega[i] *= m / total;

  const double n_end = rc.n.back();
  std::vector<std::vector<double>> f(3, std::vector<double>(m, 0.0));
  bool fitted[3];
  for (int c = 0; c < 3; c++) {
    const CurveFit fit = fit_curve(c, ln.data(), s, omega.data(), m);
    fitted[c] = fit.ok;
    if (!fit.ok) continue;
    out.S[c] = curve_tail(c, fit.th, n_end);
    if (std::isnan(out.S[c])) { fitted[c] = false; continue; }
    for (int i = 0; i < m; i++) f[c][i] = curve_eval(c, ln[i], fit.th, nullptr);
  }
  if (!fitted[LOGNORMAL] || !(out.S[LOGNORMAL] > 0.0) || std::isinf(out.S[LOGNORMAL])) return out;

  double cap[3];
  weight_bounds(out.S, s_obs, cap);
  for (int c = 0; c < 3; c++) if (!fitted[c]) cap[c] = 0.0;
  optimise_weights(f, s, omega.data(), m, cap, out.w);

  out.unobs = 0.0;
  for (int c = 0; c < 3; c++) if (out.w[c] > 0.0) out.unobs += out.w[c] * out.S[c];
  out.ok = std::isfinite(out.unobs);
  return out;
}

// ---- smoothing of the counts (Eq. 10) ------------------------------------

// Number of tags represented by y_full; the y_K+ slot counts K + 1 each.
inline double count_tags(const std::vector<double>& y) {
  double t = 0.0;
  for (std::size_t k = 0; k < y.size(); k++) t += (k + 1) * y[k];
  return t;
}

// Mix the first contiguous positive counts y_1..y_k1 (k1 <= K) with a
// log-linear fit log y_k = b0 + b1 k (weights y_k) in proportion p, keeping
// the number of tags. Returns false when fewer than two counts are positive.
inline bool smooth_counts(const std::vector<double>& y_full, double p, std::vector<double>& out) {
  const std::size_t K = y_full.size() - 1;
  std::size_t k1 = 0;
  while (k1 < K && y_full[k1] > 0.0) k1++;
  out = y_full;
  if (k1 < 2) return false;

  double sw = 0, sx = 0, sy = 0, sxx = 0, sxy = 0;
  for (std::size_t i = 0; i < k1; i++) {
    const double x = i + 1.0, yv = std::log(y_full[i]), wt = y_full[i];
    sw += wt; sx += wt * x; sy += wt * yv; sxx += wt * x * x; sxy += wt * x * yv;
  }
  const double b1 = (sw * sxy - sx * sy) / (sw * sxx - sx * sx);
  const double b0 = (sy - b1 * sx) / sw;

  std::vector<double> fitted(k1), head(y_full.begin(), y_full.begin() + k1);
  for (std::size_t i = 0; i < k1; i++) fitted[i] = std::exp(b0 + b1 * (i + 1.0));
  const double scale = p * count_tags(head) / count_tags(fitted);
  for (std::size_t i = 0; i < k1; i++) out[i] = y_full[i] * (1.0 - p) + fitted[i] * scale;
  const double rescale = count_tags(y_full) / count_tags(out);
  for (double& v : out) v *= rescale;
  return true;
}

// ---- point estimate ------------------------------------------------------

struct Estimate {
  double estimate;  // S_obs + undetected species
  double S[3];      // undetected species per curve, averaged over portions
  double w[3];      // weights, averaged over portions
  int n_portions;   // portions that produced an estimate
};

// Undetected species averaged over the curve portions; y_full = (y_1..y_K, y_K+).
inline Estimate estimate_counts(const std::vector<double>& y_full, double s_obs, int max_points) {
  Estimate e;
  e.estimate = kNaN; e.n_portions = 0;
  for (int c = 0; c < 3; c++) { e.S[c] = 0.0; e.w[c] = 0.0; }

  const std::vector<double> yk(y_full.begin(), y_full.end() - 1);
  const RateCurve rc = rate_curve(yk, max_points);
  if (rc.n.size() < 3) return e;

  static const double portions[] = {0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30};
  double unobs = 0.0;
  for (double q : portions) {
    const Combination comb = fit_portion(rc, q, s_obs);
    if (!comb.ok) continue;
    unobs += comb.unobs;
    for (int c = 0; c < 3; c++) { e.S[c] += comb.S[c]; e.w[c] += comb.w[c]; }
    e.n_portions++;
  }
  if (e.n_portions == 0) return e;
  for (int c = 0; c < 3; c++) { e.S[c] /= e.n_portions; e.w[c] /= e.n_portions; }
  e.estimate = s_obs + unobs / e.n_portions;
  return e;
}

// Point estimate, optionally with the smoothing adjustment of Eq. 11,
// S_adj = (S_p - S_1) / (1 - p) + S_1, from the counts smoothed in
// proportion p and fully.
inline Estimate estimate(const std::vector<double>& y_full, double s_obs, bool adjust,
                         double p, int max_points) {
  if (!adjust) return estimate_counts(y_full, s_obs, max_points);
  std::vector<double> y_p, y_1;
  if (!smooth_counts(y_full, p, y_p) || !smooth_counts(y_full, 1.0, y_1)) {
    return estimate_counts(y_full, s_obs, max_points);
  }
  Estimate e_p = estimate_counts(y_p, s_obs, max_points);
  const Estimate e_1 = estimate_counts(y_1, s_obs, max_points);
  e_p.estimate = (e_p.estimate - e_1.estimate) / (1.0 - p) + e_1.estimate;
  return e_p;
}

// ---- Poisson mixture for the count covariances (Eq. 12-15) ----------------

// Solve M z = r for the symmetric system on the index set idx (Gaussian
// elimination with partial pivoting). Returns false when singular.
inline bool solve_subsystem(const std::vector<double>& H, const std::vector<double>& g, int J,
                            double u, const std::vector<int>& idx, std::vector<double>& z) {
  const int p = static_cast<int>(idx.size());
  std::vector<double> M(static_cast<std::size_t>(p) * (p + 1));
  for (int i = 0; i < p; i++) {
    for (int j = 0; j < p; j++) M[i * (p + 1) + j] = H[idx[i] * J + idx[j]] + (i == j ? u : 0.0);
    M[i * (p + 1) + p] = g[idx[i]];
  }
  for (int col = 0; col < p; col++) {
    int piv = col;
    for (int r = col + 1; r < p; r++)
      if (std::fabs(M[r * (p + 1) + col]) > std::fabs(M[piv * (p + 1) + col])) piv = r;
    if (std::fabs(M[piv * (p + 1) + col]) < 1e-300) return false;
    if (piv != col)
      for (int j = 0; j <= p; j++) std::swap(M[col * (p + 1) + j], M[piv * (p + 1) + j]);
    for (int r = col + 1; r < p; r++) {
      const double factor = M[r * (p + 1) + col] / M[col * (p + 1) + col];
      for (int j = col; j <= p; j++) M[r * (p + 1) + j] -= factor * M[col * (p + 1) + j];
    }
  }
  z.assign(p, 0.0);
  for (int i = p - 1; i >= 0; i--) {
    double acc = M[i * (p + 1) + p];
    for (int j = i + 1; j < p; j++) acc -= M[i * (p + 1) + j] * z[j];
    z[i] = acc / M[i * (p + 1) + i];
  }
  return true;
}

// v >= 0 minimising ||y - A v||^2 + u ||v||^2, from H = A'A and g = A'y, by
// the active-set method of Lawson and Hanson in the normal-equation form of
// Bro and de Jong (1997).
inline void ridge_nnls(const std::vector<double>& H, const std::vector<double>& g, int J,
                       double u, std::vector<double>& v) {
  v.assign(J, 0.0);
  std::vector<char> active(J, 0);
  std::vector<double> grad(g), z, s(J);
  std::vector<int> idx;
  double tol = 0.0;
  for (int j = 0; j < J; j++) tol = std::max(tol, std::fabs(g[j]));
  tol = 1e-12 * std::max(tol, 1.0);

  for (int outer = 0; outer < 3 * J; outer++) {
    int enter = -1;
    double best = tol;
    for (int j = 0; j < J; j++) if (!active[j] && grad[j] > best) { best = grad[j]; enter = j; }
    if (enter < 0) break;
    active[enter] = 1;

    for (int inner = 0; inner < 3 * J; inner++) {
      idx.clear();
      for (int j = 0; j < J; j++) if (active[j]) idx.push_back(j);
      if (!solve_subsystem(H, g, J, u, idx, z)) { active[enter] = 0; break; }
      std::fill(s.begin(), s.end(), 0.0);
      bool positive = true;
      for (std::size_t i = 0; i < idx.size(); i++) {
        s[idx[i]] = z[i];
        if (z[i] <= 0.0) positive = false;
      }
      if (positive) { v = s; break; }
      // step towards s until the first active coefficient reaches zero
      double alpha = 1.0;
      for (int j : idx) if (s[j] <= 0.0) alpha = std::min(alpha, v[j] / (v[j] - s[j]));
      for (int j = 0; j < J; j++) {
        v[j] += alpha * (s[j] - v[j]);
        if (active[j] && v[j] <= 1e-14) { v[j] = 0.0; active[j] = 0; }
      }
    }
    for (int j = 0; j < J; j++) {
      double hv = u * v[j];
      for (int l = 0; l < J; l++) hv += H[j * J + l] * v[l];
      grad[j] = g[j] - hv;
    }
  }
}

}  // namespace arc
}  // namespace spacc

#endif  // SPACC_CORE_ARC_H
