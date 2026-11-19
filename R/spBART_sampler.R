################################################################################
# spBART_sampler.R
#
# Shared spBART Gibbs sampler, sourced by every simulation and real-data script
# so that all analyses run exactly the same code.
#
# Model (Section 2.2 of the manuscript):
#   Y_i = I(U_i > 0),  U_i = c + f(X_i) + Z_i' beta + eps_i,  eps_i ~ N(0, 1)
#   f      : sum of n_trees trees with the DART prior (Linero 2018)
#   beta   : N(0, sigma_beta^2 I)
#   c      : fixed offset Phi^{-1}(mean(Y)), the same centring pbart() uses
#            (binaryOffset); it is a constant, not a parameter. c = 0 when
#            use_offset = FALSE.
#
# Gibbs sampler (Appendix A):
#   1. U    | Y, f, beta  truncated normals (Albert and Chib 1993)
#   2. beta | U, f        conjugate Gaussian
#   3. f    | U, beta     one BART backfitting sweep on U - c - Z'beta
#                         (stateful sampler in the modified BART package)
#
# Settings shared with the BART comparator pbart(), so the two are compared on
# equal footing:
#   leaf prior s.d. 3/(k*sqrt(n_trees)) with k = 2; numcut = 100;
#   usequants = FALSE; power = 2, base = 0.95; a = 0.5, b = 1;
#   DART switched on at iteration floor(n_burn/2) + 1 (pbart: i == burn/2,
#   0-based); latent error s.d. fixed at 1.
#
# Two analysis-specific settings (see spBART_analysis_settings()):
#   "simulation"  : offset c = Phi^{-1}(mean(Y)) and leaf prior s.d.
#                   3/(k*sqrt(n_trees)) = 1.5/sqrt(T), as in pbart(); these are
#                   the settings of the simulation studies (Section 5).
#   "application" : no offset (c = 0) and leaf prior s.d.
#                   0.5*sd(U0 - Z beta0)/sqrt(T), about 0.5/sqrt(T), computed
#                   from the initial latent variables; these are the settings of
#                   the code that produced the MM application (Section 4).
# fit_spBART() defaults to "simulation"; fit_spBART_analysis() selects either.
#
# Requires the modified BART package, version 2.9.9.1
# (BART_package/BART_2.9.9.1.tar.gz).
################################################################################

SPBART_BART_VERSION <- "2.9.9.1"

#' Personal R library used for the modified BART package (same rule on the
#' laptop and on the cluster): R_LIBS_USER if set, else ~/R/library/<major.minor>
spbart_user_lib <- function() {
  lib <- Sys.getenv("R_LIBS_USER")
  if (lib == "") {
    lib <- file.path(Sys.getenv("HOME"), "R", "library",
                     paste0(R.version$major, ".", strsplit(R.version$minor, "\\.")[[1]][1]))
  }
  path.expand(lib)
}

#' Load the modified BART package and stop if the wrong build is installed
load_spBART_BART <- function() {
  lib <- spbart_user_lib()
  if (dir.exists(lib)) .libPaths(unique(c(lib, .libPaths())))
  extra <- Sys.getenv("SPBART_EXTRA_LIBS")   # optional, local testing only (":"-separated)
  if (nzchar(extra)) .libPaths(unique(c(strsplit(extra, ":")[[1]], .libPaths())))
  suppressPackageStartupMessages(library(BART))
  installed <- as.character(utils::packageVersion("BART"))
  needed <- c("wbart_create", "wbart_update_response", "wbart_run_iteration",
              "wbart_start_dart", "wbart_destroy")
  missing <- needed[!vapply(needed, exists, logical(1))]
  if (installed != SPBART_BART_VERSION || length(missing) > 0) {
    stop(sprintf(paste0(
      "Modified BART %s is required but BART %s is loaded%s.\n",
      "Install it once with:  Rscript BART_package/00_install_BART.R\n",
      "(or submit slurm/00_install_BART.sbatch)."),
      SPBART_BART_VERSION, installed,
      if (length(missing)) paste0(" (missing: ", paste(missing, collapse = ", "), ")") else ""))
  }
  invisible(installed)
}

#' Fit spBART
#'
#' @param X        n x D matrix for the tree ensemble (genes, plus study indicator)
#' @param Z        n x J matrix for the linear part (no intercept column)
#' @param Y        binary outcome (0/1)
#' @param X_test,Z_test optional test matrices; posterior predictive
#'                 probabilities are returned for them
#' @param n_burn,n_iter,n_thin burn-in, number of saved draws, thinning
#' @param sigma_beta prior s.d. of beta (manuscript: sigma_beta^2 = 10)
#' @param beta_init starting value of beta (default 0); used for dispersed
#'   starts in multi-chain runs
#' @param use_offset if TRUE, centre the latent scale at the fixed offset
#'   qnorm(mean(Y)) (pbart's binaryOffset); if FALSE, no offset
#' @param leaf_prior "pbart": leaf prior s.d. 3/(k*sqrt(n_trees));
#'   "application": 0.5*sd(U0 - Z beta0)/sqrt(n_trees), from the initial
#'   latent variables U0 and coefficients beta0 (k is then not used)
#' @param seed if not NULL, set.seed(seed) before any random draw
#' @return list with beta_draws, f_train_draws, varcount (n_iter x D split
#'   counts, for PIPs), loglik (per saved draw), offset, f_hat, beta_hat,
#'   mu_hat (posterior mean of the latent index), prob_train_mean, and, when
#'   test data are given, f_test_draws, prob.test (n_iter x n_test) and
#'   prob.test.mean
fit_spBART <- function(X, Z, Y,
                       X_test = NULL, Z_test = NULL,
                       n_burn = 2000, n_iter = 1000, n_thin = 5,
                       n_trees = 200, sparse = TRUE,
                       theta = 0, omega = 1,
                       a = 0.5, b = 1, rho = NULL,
                       k = 2, numcut = 100, usequants = FALSE,
                       sigma_beta = sqrt(10),
                       use_offset = TRUE,
                       leaf_prior = c("pbart", "application"),
                       beta_init = NULL, seed = NULL,
                       verbose = TRUE) {
  leaf_prior <- match.arg(leaf_prior)
  if (!is.null(seed)) set.seed(seed)

  X <- as.matrix(X)
  Z <- as.matrix(Z)
  Y <- as.integer(Y)
  N <- length(Y)
  D <- ncol(X)
  J <- ncol(Z)
  if (nrow(X) != N || nrow(Z) != N) stop("X, Z and Y must have the same number of rows")
  if (is.null(rho)) rho <- D

  has_test <- !is.null(X_test) && !is.null(Z_test)
  if (has_test) {
    X_test <- as.matrix(X_test)
    Z_test <- as.matrix(Z_test)
    N_test <- nrow(X_test)
  }

  # Fixed latent-scale offset, as pbart(): binaryOffset = qnorm(mean(y))
  offset <- if (use_offset) qnorm(mean(Y)) else 0

  # Truncation limits for U | Y
  lower <- ifelse(Y == 1, 0, -Inf)
  upper <- ifelse(Y == 1, Inf, 0)

  # beta | U, f has a fixed posterior covariance (Z does not change)
  V_beta <- solve(crossprod(Z) + diag(J) / sigma_beta^2)
  R_beta <- chol(V_beta)            # V_beta = t(R_beta) %*% R_beta
  VZt <- V_beta %*% t(Z)

  # Initial values
  beta <- if (is.null(beta_init)) rep(0, J) else as.numeric(beta_init)
  U <- rnorm(N)
  f_train <- rep(0, N)
  Zbeta <- as.numeric(Z %*% beta)

  # Storage
  beta_draws <- matrix(NA_real_, n_iter, J, dimnames = list(NULL, colnames(Z)))
  f_draws <- matrix(NA_real_, n_iter, N)
  varcount <- matrix(0L, n_iter, D, dimnames = list(NULL, colnames(X)))
  loglik <- rep(NA_real_, n_iter)
  if (has_test) {
    f_test_draws <- matrix(NA_real_, n_iter, N_test)
    prob_test_draws <- matrix(NA_real_, n_iter, N_test)
  }

  # Leaf prior s.d.
  tau <- if (leaf_prior == "pbart") {
    3 / (k * sqrt(n_trees))
  } else {
    0.5 * sd(U - offset - Zbeta) / sqrt(n_trees)
  }

  # Stateful BART sampler (created once; trees persist across iterations)
  create_args <- list(
    x.train = X, y.train = U - offset - Zbeta,
    sparse = sparse, theta = theta, omega = omega,
    a = a, b = b, rho = rho, augment = FALSE,
    ntree = n_trees, numcut = numcut, usequants = usequants,
    k = k, tau = tau
  )
  if (has_test) create_args$x.test <- X_test
  bart_sampler <- do.call(wbart_create, create_args)
  on.exit(wbart_destroy(bart_sampler), add = TRUE)
  kept <- attr(bart_sampler, "rm.const")   # columns of X used by BART

  dart_start <- floor(n_burn / 2) + 1
  n_total <- n_burn + n_iter * n_thin

  for (iter in seq_len(n_total)) {
    if (verbose && iter %% 500 == 0) message("  Iteration ", iter, "/", n_total)

    # 1. U | Y, f, beta
    U <- truncnorm::rtruncnorm(N, a = lower, b = upper,
                               mean = offset + f_train + Zbeta, sd = 1)

    # 2. beta | U, f
    beta <- as.numeric(VZt %*% (U - offset - f_train) + t(R_beta) %*% rnorm(J))
    Zbeta <- as.numeric(Z %*% beta)

    # 3. f | U, beta (DART switched on at the same point as pbart)
    if (sparse && iter == dart_start) wbart_start_dart(bart_sampler)
    wbart_update_response(bart_sampler, U - offset - Zbeta)
    bart_result <- wbart_run_iteration(bart_sampler)
    f_train <- bart_result$yhat.train

    # Save after burn-in, every n_thin-th iteration
    if (iter > n_burn && (iter - n_burn) %% n_thin == 0) {
      q <- (iter - n_burn) / n_thin
      beta_draws[q, ] <- beta
      f_draws[q, ] <- f_train
      varcount[q, kept] <- bart_result$varcount
      eta <- offset + f_train + Zbeta
      loglik[q] <- sum(pnorm(eta[Y == 1], log.p = TRUE)) +
                   sum(pnorm(-eta[Y == 0], log.p = TRUE))
      if (has_test) {
        f_test_draws[q, ] <- bart_result$yhat.test
        prob_test_draws[q, ] <- pnorm(offset + bart_result$yhat.test +
                                        as.numeric(Z_test %*% beta))
      }
    }
  }

  f_hat <- colMeans(f_draws)
  beta_hat <- colMeans(beta_draws)
  mu_draws <- offset + f_draws + beta_draws %*% t(Z)

  result <- list(
    beta_draws = beta_draws,
    f_train_draws = f_draws,
    varcount = varcount,
    loglik = loglik,
    offset = offset,
    f_hat = f_hat,
    beta_hat = beta_hat,
    mu_hat = colMeans(mu_draws),
    prob_train_mean = colMeans(pnorm(mu_draws)),
    n_burn = n_burn, n_iter = n_iter, n_thin = n_thin,
    settings = list(n_trees = n_trees, k = k, tau = tau,
                    use_offset = use_offset, leaf_prior = leaf_prior,
                    numcut = numcut, sparse = sparse, a = a, b = b, rho = rho,
                    sigma_beta = sigma_beta, dart_start = dart_start,
                    BART_version = as.character(utils::packageVersion("BART")))
  )
  if (has_test) {
    result$f_test_draws <- f_test_draws
    result$prob.test <- prob_test_draws
    result$prob.test.mean <- colMeans(prob_test_draws)
  }
  result
}

#' Offset and leaf-prior settings of each analysis in the manuscript
#'
#' @param analysis "simulation" (Section 5: offset qnorm(mean(Y)), leaf prior
#'   s.d. 3/(k*sqrt(T)), as in pbart) or "application" (Section 4: no offset,
#'   leaf prior s.d. 0.5*sd(U0 - Z beta0)/sqrt(T))
#' @return list(use_offset, leaf_prior), arguments of fit_spBART()
spBART_analysis_settings <- function(analysis = c("simulation", "application")) {
  analysis <- match.arg(analysis)
  switch(analysis,
         simulation  = list(use_offset = TRUE,  leaf_prior = "pbart"),
         application = list(use_offset = FALSE, leaf_prior = "application"))
}

#' Fit spBART with the offset and leaf prior of a given analysis
#'
#' @param ... arguments of fit_spBART() other than use_offset and leaf_prior
#' @param analysis "simulation" or "application"; see spBART_analysis_settings()
fit_spBART_analysis <- function(..., analysis = c("simulation", "application")) {
  args <- list(...)
  if (any(c("use_offset", "leaf_prior") %in% names(args))) {
    stop("use_offset and leaf_prior are set by 'analysis'; call fit_spBART() to set them directly")
  }
  do.call(fit_spBART, c(args, spBART_analysis_settings(analysis)))
}

#' Fit several independent chains and pool their draws
#'
#' Each chain is a call of fit_spBART_analysis() with its own seed, and hence
#' its own initial latent variables. The pooled result has the same fields as
#' fit_spBART() (draws stacked chain by chain), plus `chains` (each chain's
#' fit) and `diagnostics` (spBART_chain_diagnostics()).
#'
#' @param ... arguments of fit_spBART_analysis(), including X, Z, Y and analysis
#' @param n_chains number of chains
#' @param seeds one seed per chain
#' @param latent_ids training participants whose latent means are monitored
#'   (default: the first, middle and last)
#' @param cores number of chains run in parallel (parallel::mclapply); the
#'   result does not depend on it, because each chain sets its own seed
fit_spBART_chains <- function(..., n_chains = 4, seeds, latent_ids = NULL, cores = 1) {
  stopifnot(length(seeds) == n_chains)
  args <- list(...)
  if ("seed" %in% names(args)) stop("seeds are set per chain through 'seeds'")
  run_chain <- function(ch) do.call(fit_spBART_analysis, c(args, list(seed = seeds[ch])))
  chains <- if (cores > 1) parallel::mclapply(seq_len(n_chains), run_chain, mc.cores = cores) else
    lapply(seq_len(n_chains), run_chain)
  failed <- vapply(chains, inherits, logical(1), what = "try-error")
  if (any(failed)) stop("chain(s) ", paste(which(failed), collapse = ", "), " failed: ", chains[[which(failed)[1]]])
  stack <- function(nm) if (is.null(chains[[1]][[nm]])) NULL else do.call(rbind, lapply(chains, `[[`, nm))
  Z <- as.matrix(args$Z)
  pooled <- chains[[1]]
  for (nm in c("beta_draws", "f_train_draws", "varcount", "f_test_draws", "prob.test")) pooled[nm] <- list(stack(nm))
  pooled$loglik <- unlist(lapply(chains, `[[`, "loglik"))
  mu_draws <- pooled$offset + pooled$f_train_draws + pooled$beta_draws %*% t(Z)
  pooled$f_hat <- colMeans(pooled$f_train_draws)
  pooled$beta_hat <- colMeans(pooled$beta_draws)
  pooled$mu_hat <- colMeans(mu_draws)
  pooled$prob_train_mean <- colMeans(pnorm(mu_draws))
  if (!is.null(pooled$prob.test)) pooled$prob.test.mean <- colMeans(pooled$prob.test)
  pooled$n_iter <- nrow(pooled$beta_draws)
  pooled$n_chains <- n_chains
  pooled$seeds <- seeds
  pooled$chains <- chains
  if (is.null(latent_ids)) latent_ids <- unique(c(1, ceiling(nrow(Z) / 2), nrow(Z)))
  pooled$diagnostics <- spBART_chain_diagnostics(chains, Z, latent_ids)
  pooled
}

#' Convergence diagnostics across chains
#'
#' Rank-normalized split-Rhat and bulk effective sample size (Vehtari et al.
#' 2021; R package posterior) for the log-likelihood, the total number of tree
#' splits, each coefficient and the latent means of the participants latent_ids.
#'
#' @param chains list of fit_spBART() results with the same number of draws
#' @param Z the linear-part design matrix used in the fits
#' @param latent_ids indices of the monitored participants
#' @param rhat_max,ess_min convergence criterion
#' @return data frame with one row per quantity: rhat, ess_bulk and pass
spBART_chain_diagnostics <- function(chains, Z, latent_ids, rhat_max = 1.02, ess_min = 390) {
  per_chain <- function(get) sapply(chains, get)                # draws x chains
  q <- list(loglik = per_chain(function(ch) ch$loglik),
            n_splits = per_chain(function(ch) rowSums(ch$varcount)))
  Z <- as.matrix(Z)
  for (j in seq_len(ncol(Z))) {
    nm <- if (is.null(colnames(Z))) paste0("beta_", j) else paste0("beta_", colnames(Z)[j])
    q[[nm]] <- per_chain(function(ch) ch$beta_draws[, j])
  }
  for (i in latent_ids) {
    q[[paste0("latent_mean_", i)]] <- per_chain(function(ch)
      ch$offset + ch$f_train_draws[, i] + as.numeric(ch$beta_draws %*% Z[i, ]))
  }
  out <- data.frame(quantity = names(q),
                    rhat = vapply(q, posterior::rhat, numeric(1)),
                    ess_bulk = vapply(q, posterior::ess_bulk, numeric(1)),
                    row.names = NULL)
  out$pass <- out$rhat <= rhat_max & out$ess_bulk >= ess_min
  out
}

#' Bayesian FDR selection from PIPs (Section 3.3; Mueller et al. 2006)
fdr_gene_selection <- function(PIPs, alpha = 0.05) {
  PIP_sorted <- sort(PIPs, decreasing = TRUE)
  D <- length(PIP_sorted)
  if (D == 0) {
    return(list(selected_genes = character(0), R_star = 0, tau = NA, FDR_hat = NA))
  }
  FDR <- cumsum(1 - PIP_sorted) / seq_len(D)
  valid_R <- which(FDR <= alpha)
  if (length(valid_R) == 0) {
    return(list(selected_genes = character(0), R_star = 0, tau = NA, FDR_hat = NA))
  }
  R_star <- max(valid_R)
  list(selected_genes = names(PIP_sorted)[seq_len(R_star)],
       R_star = R_star,
       tau = PIP_sorted[R_star],
       FDR_hat = FDR[R_star])
}
