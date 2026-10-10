// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include "core/arc_core.h"
#include "parallel_utils.h"
using namespace Rcpp;
using namespace RcppParallel;


// ARC point estimate from frequency counts y_full = (y_1..y_K, y_K+).
// [[Rcpp::export]]
List cpp_arc_point(NumericVector y_full, double s_obs, bool adjust, double p, int max_points) {
  const std::vector<double> y(y_full.begin(), y_full.end());
  const spacc::arc::Estimate e = spacc::arc::estimate(y, s_obs, adjust, p, max_points);
  NumericVector S = NumericVector::create(e.S[0], e.S[1], e.S[2]);
  NumericVector w = NumericVector::create(e.w[0], e.w[1], e.w[2]);
  const CharacterVector curves = CharacterVector::create("lognormal", "loglogistic", "weibull");
  S.attr("names") = curves;
  w.attr("names") = curves;
  return List::create(Named("estimate") = e.estimate, Named("undetected") = S,
                      Named("weights") = w, Named("n_portions") = e.n_portions);
}


// ARC point estimate of every row of a matrix of frequency counts; each row's
// S_obs is its sum.
struct ArcDrawsWorker : public Worker {
  const RMatrix<double> draws;
  RVector<double> out;
  const bool adjust;
  const double p;
  const int max_points;

  ArcDrawsWorker(const NumericMatrix& draws_, NumericVector& out_, bool adjust_, double p_,
                 int max_points_)
    : draws(draws_), out(out_), adjust(adjust_), p(p_), max_points(max_points_) {}

  void operator()(std::size_t begin, std::size_t end) {
    std::vector<double> y(draws.ncol());
    for (std::size_t b = begin; b < end; b++) {
      double s_obs = 0.0;
      for (std::size_t k = 0; k < y.size(); k++) { y[k] = draws(b, k); s_obs += y[k]; }
      out[b] = spacc::arc::estimate(y, s_obs, adjust, p, max_points).estimate;
    }
  }
};

// [[Rcpp::export]]
NumericVector cpp_arc_draws(NumericMatrix draws, bool adjust, double p, int max_points,
                            int n_cores = 1) {
  NumericVector out(draws.nrow());
  ArcDrawsWorker worker(draws, out, adjust, p, max_points);
  run_workers(static_cast<std::size_t>(draws.nrow()), n_cores, worker);
  return out;
}


// Poisson-mixture weights for each penalty u (Eq. 13): fitted to every
// in-sample row with the out-of-sample error, and to the observed counts.
struct ArcPenaltyWorker : public Worker {
  const RMatrix<double> A, in, out;
  const RVector<double> penalties, y_obs;
  RMatrix<double> err, v_obs;
  std::vector<double> H;

  ArcPenaltyWorker(const NumericMatrix& A_, const NumericMatrix& in_, const NumericMatrix& out_,
                   const NumericVector& penalties_, const NumericVector& y_obs_,
                   NumericMatrix& err_, NumericMatrix& v_obs_)
    : A(A_), in(in_), out(out_), penalties(penalties_), y_obs(y_obs_), err(err_), v_obs(v_obs_) {
    const int R = A.nrow(), J = A.ncol();
    H.assign(static_cast<std::size_t>(J) * J, 0.0);
    for (int j = 0; j < J; j++)
      for (int l = 0; l < J; l++)
        for (int r = 0; r < R; r++) H[j * J + l] += A(r, j) * A(r, l);
  }

  template<typename Row>
  void fit(const Row& y, double u, std::vector<double>& g, std::vector<double>& v) const {
    const int R = A.nrow(), J = A.ncol();
    for (int j = 0; j < J; j++) {
      g[j] = 0.0;
      for (int r = 0; r < R; r++) g[j] += A(r, j) * y[r];
    }
    spacc::arc::ridge_nnls(H, g, J, u, v);
  }

  void operator()(std::size_t begin, std::size_t end) {
    const int R = A.nrow(), J = A.ncol(), B = in.nrow();
    std::vector<double> g(J), v(J), y(R);
    for (std::size_t ui = begin; ui < end; ui++) {
      const double u = penalties[ui];
      for (int b = 0; b < B; b++) {
        for (int r = 0; r < R; r++) y[r] = in(b, r);
        fit(y, u, g, v);
        double sse = 0.0;
        for (int r = 0; r < R; r++) {
          double fitted = 0.0;
          for (int j = 0; j < J; j++) fitted += A(r, j) * v[j];
          const double d = fitted - out(b, r);
          sse += d * d;
        }
        err(b, ui) = sse / R;
      }
      for (int r = 0; r < R; r++) y[r] = y_obs[r];
      fit(y, u, g, v);
      for (int j = 0; j < J; j++) v_obs(j, ui) = v[j];
    }
  }
};

// [[Rcpp::export]]
List cpp_arc_penalty(NumericMatrix A, NumericMatrix in_sample, NumericMatrix out_sample,
                     NumericVector penalties, NumericVector y_obs, int n_cores = 1) {
  NumericMatrix err(in_sample.nrow(), penalties.size());
  NumericMatrix v_obs(A.ncol(), penalties.size());
  ArcPenaltyWorker worker(A, in_sample, out_sample, penalties, y_obs, err, v_obs);
  run_workers(static_cast<std::size_t>(penalties.size()), n_cores, worker);
  return List::create(Named("err") = err, Named("v_obs") = v_obs);
}
