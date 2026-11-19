/*
 *  BART: Bayesian Additive Regression Trees - Stateful Sampler
 *  Modified for the semi-parametric probit BART model (spBART)
 *
 *  This file implements a stateful BART sampler that can be used within
 *  Gibbs samplers for semi-parametric models. Unlike the standard wbart()
 *  which runs a complete MCMC chain in one call, this implementation:
 *
 *  1. Creates a sampler object ONCE (wbart_create)
 *  2. Allows updating the response without reinitializing trees (wbart_update_response)
 *  3. Runs ONE iteration at a time (wbart_run_iteration)
 *  4. Switches on the Dirichlet (DART) splitting-probability update when
 *     requested (wbart_start_dart), mirroring pbart(), which calls
 *     startdart() at iteration burn/2
 *  5. Cleans up when done (wbart_destroy)
 *
 *  Version 2.9.9.1 changes (2026-10-11):
 *   - Error s.d. fixed at 1 (probit latent scale); previously sd(initial response).
 *   - R's RNG state is synchronised around every iteration (GetRNGstate /
 *     PutRNGstate); previously C++ draws did not advance .Random.seed, so the
 *     next R-level draws reused the same random numbers.
 *   - Cutpoints from bartModelMatrix() are passed to the model via setxinfo(),
 *     exactly as cpbart.cpp does; previously they were computed but ignored.
 *   - DART is switched on by an explicit call (cwbart_start_dart) instead of
 *     at creation, so the driver can use the same schedule as pbart().
 */

#include "tree.h"
#include "treefuns.h"
#include "info.h"
#include "bartfuns.h"
#include "bd.h"
#include "bart.h"
#include "heterbart.h"
#include "rn.h"

#ifndef NoRcpp

// Stateful sampler state structure
struct wbart_stateful_state {
  // Core BART objects
  heterbart *bm;        // BART model (tree ensemble)

  // Data pointers
  double *ix;           // Predictors (training) - FIXED
  double *ixp;          // Predictors (test) - FIXED or NULL
  double *iy;           // Response (training) - MUTABLE
  double *iw;           // Weights - FIXED

  // Dimensions
  size_t n;             // Number of training observations
  size_t p;             // Number of predictors
  size_t np;            // Number of test observations
  size_t m;             // Number of trees

  // MCMC parameters
  double sigma;         // Error s.d. on the latent scale (fixed at 1 for probit)
  double *allfit;       // Current fit (sum of all trees)
  double *r;            // Residuals
  double *ftemp;        // Temporary fit storage

  // Test set predictions
  double *allfit_test;  // Test set predictions (if test data provided)
  double *ftemp_test;   // Temporary test fit storage

  // Prior/model parameters
  double alpha;         // BART tree prior (base)
  double mybeta;        // BART tree prior (power)
  double tau;           // Prior s.d. of leaf parameters
  bool dart;            // Whether Dirichlet sparsity is requested
  bool dart_on;         // Whether the Dirichlet update has been switched on
  double a;             // Dirichlet parameter
  double b;             // Dirichlet parameter
  double rho;           // Dirichlet concentration

  // RNG
  arn gen;              // R-based random number generator

  // Cutpoint information
  xinfo xi;             // Cutpoint information

  // Variable counts for PIP computation
  std::vector<int> varcnt;
};

//=============================================================================
// Create stateful BART sampler
//=============================================================================
RcppExport SEXP cwbart_create(
  SEXP _in,         // number of observations in training data
  SEXP _ip,         // dimension of x
  SEXP _inp,        // number of observations in test data
  SEXP _ix,         // x, train, pxn (transposed)
  SEXP _iy,         // y, train, nx1
  SEXP _ixp,        // x, test, pxnp (transposed)
  SEXP _im,         // number of trees
  SEXP _inc,        // number of cut points
  SEXP _ipower,     // power parameter for tree prior
  SEXP _ibase,      // base parameter for tree prior
  SEXP _itau,       // prior s.d. of leaf parameters
  SEXP _idart,      // whether to use DART
  SEXP _itheta,     // theta for DART
  SEXP _iomega,     // omega for DART
  SEXP _ia,         // a parameter for DART
  SEXP _ib,         // b parameter for DART
  SEXP _irho,       // rho parameter for DART
  SEXP _iaug,       // whether to use data augmentation
  SEXP _Xinfo       // cutpoint information (p x max(numcut) matrix)
)
{
  //--------------------------------------------------------------------------
  // Process arguments
  //--------------------------------------------------------------------------
  size_t n = Rcpp::as<int>(_in);
  size_t p = Rcpp::as<int>(_ip);
  size_t np = Rcpp::as<int>(_inp);

  Rcpp::NumericVector xv(_ix);
  double *ix = &xv[0];

  Rcpp::NumericVector yv(_iy);
  double *iy = &yv[0];

  double *ixp = nullptr;
  if(np > 0) {
    Rcpp::NumericMatrix xpv(_ixp);
    ixp = &xpv[0];
  }

  size_t m = Rcpp::as<int>(_im);

  Rcpp::IntegerVector _nc(_inc);
  int *numcut = &_nc[0];

  double mybeta = Rcpp::as<double>(_ipower);
  double alpha = Rcpp::as<double>(_ibase);
  double tau = Rcpp::as<double>(_itau);

  bool dart = (Rcpp::as<int>(_idart) == 1);

  double theta = Rcpp::as<double>(_itheta);
  double omega = Rcpp::as<double>(_iomega);
  double a = Rcpp::as<double>(_ia);
  double b = Rcpp::as<double>(_ib);
  double rho = Rcpp::as<double>(_irho);

  bool aug = (Rcpp::as<int>(_iaug) == 1);

  Rcpp::NumericMatrix Xinfo(_Xinfo);

  //--------------------------------------------------------------------------
  // Allocate state structure
  //--------------------------------------------------------------------------
  wbart_stateful_state *state = new wbart_stateful_state;

  state->n = n;
  state->p = p;
  state->np = np;
  state->m = m;

  // Copy data (y is copied because it is updated every iteration)
  state->ix = new double[p * n];
  std::copy(ix, ix + p * n, state->ix);

  state->iy = new double[n];
  std::copy(iy, iy + n, state->iy);

  if(np > 0) {
    state->ixp = new double[p * np];
    std::copy(ixp, ixp + p * np, state->ixp);
  } else {
    state->ixp = nullptr;
  }

  state->iw = new double[n];
  for(size_t i = 0; i < n; i++) state->iw[i] = 1.0;

  state->alpha = alpha;
  state->mybeta = mybeta;
  state->tau = tau;
  state->dart = dart;
  state->dart_on = false;
  state->a = a;
  state->b = b;
  state->rho = rho;

  //--------------------------------------------------------------------------
  // Error s.d.: fixed at 1 on the probit latent scale (Albert & Chib 1993)
  //--------------------------------------------------------------------------
  state->sigma = 1.0;

  //--------------------------------------------------------------------------
  // Allocate working arrays
  //--------------------------------------------------------------------------
  state->allfit = new double[n];
  state->r = new double[n];
  state->ftemp = new double[n];

  for(size_t i = 0; i < n; i++) {
    state->allfit[i] = 0.0;
    state->r[i] = state->iy[i];
    state->ftemp[i] = 0.0;
  }

  if(np > 0) {
    state->allfit_test = new double[np];
    state->ftemp_test = new double[np];
    for(size_t i = 0; i < np; i++) {
      state->allfit_test[i] = 0.0;
      state->ftemp_test[i] = 0.0;
    }
  } else {
    state->allfit_test = nullptr;
    state->ftemp_test = nullptr;
  }

  state->varcnt.resize(p);
  for(size_t i = 0; i < p; i++) state->varcnt[i] = 0;

  //--------------------------------------------------------------------------
  // Create BART model
  //--------------------------------------------------------------------------
  state->bm = new heterbart(m);
  state->bm->setprior(alpha, mybeta, tau);

  // Cutpoints: use the bartModelMatrix() grid, as cpbart.cpp does.
  // setxinfo() must be called before setdata(); otherwise setdata() builds
  // its own grid with makexinfo().
  if(Xinfo.size() > 0) {
    state->xi.resize(p);
    for(size_t i = 0; i < p; i++) {
      state->xi[i].resize(numcut[i]);
      for(size_t j = 0; j < (size_t)numcut[i]; j++) {
        state->xi[i][j] = Xinfo(i, j);
      }
    }
    state->bm->setxinfo(state->xi);
  }

  state->bm->setdata(p, n, state->ix, state->iy, numcut);

  // DART parameters (the update itself is switched on by cwbart_start_dart)
  if(dart) {
    state->bm->setdart(a, b, rho, aug, true, theta, omega);
  }

  Rcpp::XPtr<wbart_stateful_state> ptr(state, true);
  return ptr;
}

//=============================================================================
// Switch on the Dirichlet splitting-probability update (once)
//=============================================================================
RcppExport SEXP cwbart_start_dart(SEXP _state)
{
  Rcpp::XPtr<wbart_stateful_state> ptr(_state);
  wbart_stateful_state *state = ptr;

  if(state->dart && !state->dart_on) {
    state->bm->startdart();   // bart::startdart() toggles, so call it only once
    state->dart_on = true;
  }

  return Rcpp::wrap(state->dart_on);
}

//=============================================================================
// Update response (y) without reinitializing trees
//=============================================================================
RcppExport SEXP cwbart_update_response(SEXP _state, SEXP _iy_new)
{
  Rcpp::XPtr<wbart_stateful_state> ptr(_state);
  wbart_stateful_state *state = ptr;

  Rcpp::NumericVector yv(_iy_new);
  double *iy_new = &yv[0];

  // The bart object holds a pointer to state->iy, so copying in place updates it
  std::copy(iy_new, iy_new + state->n, state->iy);

  for(size_t i = 0; i < state->n; i++) {
    state->r[i] = state->iy[i] - state->allfit[i];
  }

  return R_NilValue;
}

//=============================================================================
// Run ONE iteration from current tree state
//=============================================================================
RcppExport SEXP cwbart_run_iteration(SEXP _state)
{
  Rcpp::XPtr<wbart_stateful_state> ptr(_state);
  wbart_stateful_state *state = ptr;

  for(size_t j = 0; j < state->p; j++) state->varcnt[j] = 0;

  // Read .Random.seed before drawing and write it back afterwards, so the
  // C++ draws advance R's RNG stream like any other R random-number call.
  GetRNGstate();

  // One backfitting sweep over all trees (birth/death moves + leaf draws),
  // followed by the Dirichlet update when DART is switched on.
  std::vector<double> svec(state->n, state->sigma);
  state->bm->draw(svec.data(), state->gen);

  PutRNGstate();

  for(size_t i = 0; i < state->n; i++) {
    state->allfit[i] = state->bm->f(i);
  }

  // Number of splits on each variable in the current ensemble
  std::vector<size_t>& nv = state->bm->getnv();
  for(size_t j = 0; j < state->p && j < nv.size(); j++) {
    state->varcnt[j] = nv[j];
  }

  if(state->np > 0) {
    state->bm->predict(state->p, state->np, state->ixp, state->allfit_test);
  }

  Rcpp::List result;

  Rcpp::NumericVector yhat_train(state->n);
  std::copy(state->allfit, state->allfit + state->n, yhat_train.begin());
  result["yhat.train"] = yhat_train;

  if(state->np > 0) {
    Rcpp::NumericVector yhat_test(state->np);
    std::copy(state->allfit_test, state->allfit_test + state->np, yhat_test.begin());
    result["yhat.test"] = yhat_test;
  }

  Rcpp::IntegerVector varcount(state->p);
  std::copy(state->varcnt.begin(), state->varcnt.end(), varcount.begin());
  result["varcount"] = varcount;

  result["sigma"] = state->sigma;
  result["dart_on"] = state->dart_on;

  return result;
}

//=============================================================================
// Destroy sampler and free memory
//=============================================================================
RcppExport SEXP cwbart_destroy(SEXP _state)
{
  Rcpp::XPtr<wbart_stateful_state> ptr(_state);
  wbart_stateful_state *state = ptr;
  if(state->bm == nullptr) return R_NilValue;   // already destroyed

  delete state->bm;
  state->bm = nullptr;

  delete[] state->ix;
  delete[] state->iy;
  delete[] state->iw;
  if(state->ixp) delete[] state->ixp;
  delete[] state->allfit;
  delete[] state->r;
  delete[] state->ftemp;
  if(state->allfit_test) delete[] state->allfit_test;
  if(state->ftemp_test) delete[] state->ftemp_test;
  state->ixp = nullptr;
  state->allfit_test = nullptr;
  state->ftemp_test = nullptr;

  // The state struct itself is freed by the XPtr finalizer
  return R_NilValue;
}

#endif
