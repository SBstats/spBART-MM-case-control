################################################################################
# study_functions.R
#
# Shared code for simulation studies 07 to 12 (07_ to 12_ drivers):
#   - data-generating models (DGM 2 baseline, smooth DGM 3, correlated genes,
#     covariates correlated with signal genes, partial misspecification,
#     BMI-effect grid, semi-synthetic data built on the real 5hmC matrix)
#   - method wrappers with one common protocol: a single full-data fit, each
#     method's own selection rule, and prediction on an independent test set
#       spBART      PIP + Bayesian FDR (0.05); beta from the linear part
#       BART        pbart(sparse = TRUE) on [X, Z]; PIP + Bayesian FDR
#       stan4bart   probit stan4bart; PIP (varcount) + Bayesian FDR; beta
#       enet        logistic elastic net (alpha = 0.5), Z unpenalized, lambda.1se
#       knockoff    model-X knockoffs (own implementation of the equicorrelated
#                   construction), lasso statistic with Z unpenalized,
#                   target FDR 0.05 (offset 0, see m_knockoff)
#   - metrics, probit + GMM screening (as in the application), CV-union pipeline
#
# Requires: R/spBART_sampler.R sourced first.
################################################################################

suppressPackageStartupMessages({
  library(pROC)
  library(truncnorm)
})

#-------------------------------------------------------------------------------
# Constants (DGM 2 of the manuscript)
#-------------------------------------------------------------------------------
STUDY_BETA  <- c(age = 0.50, sex = 0.60, race = -0.40, bmi = 0.45)
STUDY_CUT_Q <- c(0.45, 0.55, 0.40, 0.60, 0.50, 0.50, 0.35, 0.65, 0.42, 0.58)
SD_TREES    <- 1.7861                 # population s.d. of the DGM 2 gene trees
SMOOTH_SCALE <- SD_TREES / 1.5375     # smooth f rescaled to the same s.d.
LINEAR_SCALE <- SD_TREES / sqrt(10)   # near-linear f rescaled to the same s.d.
BLOCK_SIZE  <- 50

#-------------------------------------------------------------------------------
# Signal functions of the 10 signal genes S (n x 10)
#-------------------------------------------------------------------------------
f_trees <- function(S, cuts) {           # DGM 2: 10 stumps + 5 pairwise trees
  f <- 0
  for (j in 1:10) f <- f + ifelse(S[, j] <= cuts[j], -0.40, 0.40)
  for (k in 1:5) {
    a <- 2 * k - 1; b <- 2 * k
    f <- f + ifelse(S[, a] > cuts[a] & S[, b] > cuts[b], 0.35,
                    ifelse(S[, a] <= cuts[a] & S[, b] <= cuts[b], -0.35, 0))
  }
  f
}

f_smooth <- function(S) {                # DGM 3: smooth nonlinear
  SMOOTH_SCALE * (sin(pi * S[, 1] * S[, 2] / 2) + 0.5 * (S[, 3]^2 - 1) + 0.6 * S[, 4] +
                  tanh(1.5 * S[, 5]) + 0.5 * S[, 6] * S[, 7] + (exp(-S[, 8]^2) - 1 / sqrt(3)) +
                  0.5 * (abs(S[, 9]) - sqrt(2 / pi)) + 0.4 * S[, 10])
}

f_linear <- function(S) LINEAR_SCALE * rowSums(S)   # near-linear f

#-------------------------------------------------------------------------------
# Scenario definition helper
#   type   : "dgm2" (default), "dgm1", "partial"
#   f      : "trees", "smooth", "linear"
#   rho    : within-block exchangeable correlation of genes (0 = independent)
#   layout : "default" (genes 1-10), "sep" (one signal per block), "same"
#            (5 signals in each of blocks 1 and 2)
#   corrZ  : correlation between covariates and signal genes
#   beta_bmi: overrides the BMI coefficient
#-------------------------------------------------------------------------------
scenario <- function(label, n, p, type = "dgm2", f = "trees", rho = 0,
                     layout = "default", corrZ = 0, beta_bmi = NULL) {
  list(label = label, n = n, p = p, type = type, f = f, rho = rho,
       layout = layout, corrZ = corrZ, beta_bmi = beta_bmi)
}

signal_indices <- function(p, layout) {
  idx <- switch(layout,
         default = 1:10,
         sep     = 1 + BLOCK_SIZE * (0:9),
         same    = c(1:5, BLOCK_SIZE + 1:5))
  if (max(idx) > p) stop("layout '", layout, "' needs p >= ", max(idx))
  idx
}

gen_genes <- function(n, p, rho) {
  if (rho == 0) return(matrix(rnorm(n * p), n, p))
  nb <- ceiling(p / BLOCK_SIZE)
  W <- matrix(rnorm(n * nb), n, nb)
  block <- ceiling(seq_len(p) / BLOCK_SIZE)
  sqrt(rho) * W[, block] + sqrt(1 - rho) * matrix(rnorm(n * p), n, p)
}

sigma_genes <- function(p, rho) {        # covariance of gen_genes(), for knockoffs
  if (rho == 0) return(NULL)             # NULL = identity (handled directly)
  block <- ceiling(seq_len(p) / BLOCK_SIZE)
  S <- outer(block, block, "==") * rho
  diag(S) <- 1
  S
}

gen_covariates <- function(n, S, corrZ) {
  e <- matrix(rnorm(n * 4), n, 4)
  c0 <- sqrt(1 - corrZ^2)
  cbind(age  = corrZ * S[, 1] + c0 * e[, 1],
        sex  = as.numeric(corrZ * S[, 3] + c0 * e[, 2] > 0),
        race = as.numeric(corrZ * S[, 4] + c0 * e[, 3] > qnorm(0.7)),
        bmi  = corrZ * S[, 2] + c0 * e[, 4])
}

#' Simulate one data set (training n rows + independent test set) for scenario sc
gen_study_data <- function(sc, seed, n_test = 1000) {
  one <- function(n) {
    if (sc$type == "dgm1") {             # DGM 1 of the manuscript (independent X)
      d <- generate_simulation_data(n = n, p_genes = sc$p, n_true_genes = 10,
                                    true_model = 1, seed = NULL)
      return(list(X = d$X_genes, Z = d$Z_covariates, Y = d$Y, f_true = d$f_true,
                  f_genes = NULL, beta_true = d$beta_true))
    }
    X <- gen_genes(n, sc$p, sc$rho)
    colnames(X) <- paste0("gene_", seq_len(sc$p))
    S <- X[, signal_indices(sc$p, sc$layout), drop = FALSE]
    Z <- gen_covariates(n, S, sc$corrZ)
    fg <- switch(sc$f,
                 trees  = f_trees(S, qnorm(STUDY_CUT_Q)),
                 smooth = f_smooth(S),
                 linear = f_linear(S))
    beta <- STUDY_BETA
    if (!is.null(sc$beta_bmi)) beta["bmi"] <- sc$beta_bmi
    if (sc$type == "partial") {          # age and BMI nonlinear + interacting
      ca <- qnorm(0.55); cb <- qnorm(0.40); c1 <- qnorm(0.45); c3 <- qnorm(0.40)
      h <- 0.4 * tanh(1.5 * Z[, "age"]) + 0.3 * (Z[, "bmi"]^2 - 1) +
        0.40 * ((S[, 1] > c1 & Z[, "age"] > ca) - (S[, 1] <= c1 & Z[, "age"] <= ca)) +
        0.40 * ((S[, 3] > c3 & Z[, "bmi"] > cb) - (S[, 3] <= c3 & Z[, "bmi"] <= cb))
      f_true <- fg + 0.6 * Z[, "sex"] - 0.4 * Z[, "race"] + h
      beta_true <- c(age = NA, sex = 0.6, race = -0.4, bmi = NA)
    } else {
      f_true <- fg + as.numeric(Z %*% beta)
      beta_true <- beta
    }
    list(X = X, Z = Z, Y = rbinom(n, 1, pnorm(f_true)), f_true = f_true,
         f_genes = fg, beta_true = beta_true)
  }
  set.seed(seed)
  tr <- one(sc$n)
  set.seed(seed + 5e8)
  te <- one(n_test)
  true_idx <- if (sc$type == "dgm1") 1:10 else signal_indices(sc$p, sc$layout)
  list(X = tr$X, Z = tr$Z, Y = tr$Y, f_true = tr$f_true, f_genes = tr$f_genes,
       X_test = te$X, Z_test = te$Z, Y_test = te$Y,
       beta_true = tr$beta_true,
       true_genes = colnames(tr$X)[true_idx],
       blocks = if (sc$rho > 0) setNames(ceiling(seq_len(sc$p) / BLOCK_SIZE), colnames(tr$X)) else NULL,
       Sigma_X = sigma_genes(sc$p, sc$rho),
       X_gaussian = sc$type != "semisyn")
}

#-------------------------------------------------------------------------------
# Semi-synthetic data: real development-set X and Z (n = 500) for
# training and the real validation set (n = 369) for testing; Y simulated.
#   base   : list made by real_data_analysis/00_prepare_analysis_base.R
#   set    : "screened" (115 screened genes) or "random2000"
#   signals: "random" (10 genes at random) or "cluster" (a gene + its 9 most
#            correlated genes)
#   f      : "trees" (DGM 2 with empirical-quantile cutpoints) or "smooth"
#-------------------------------------------------------------------------------
gen_semisyn_data <- function(base, set, signals, f, seed) {
  set.seed(seed)
  Xd <- base$X_dev[[set]]; Xv <- base$X_val[[set]]
  p <- ncol(Xd)
  if (signals == "random") {
    sig <- sample(p, 10)
  } else {
    g <- sample(p, 1)
    r <- abs(cor(Xd[, g], Xd))[1, ]
    r[g] <- Inf
    sig <- order(r, decreasing = TRUE)[1:10]
  }
  cuts <- apply(Xd[, sig], 2, quantile, probs = STUDY_CUT_Q)
  cuts <- diag(cuts)                        # j-th gene uses the j-th quantile level
  fun <- function(S) if (f == "trees") f_trees(S, cuts) else f_smooth(S)
  fd <- fun(Xd[, sig]); fv <- fun(Xv[, sig])
  td <- fd + as.numeric(base$Z_dev %*% STUDY_BETA)
  tv <- fv + as.numeric(base$Z_val %*% STUDY_BETA)
  list(X = Xd, Z = base$Z_dev, Y = rbinom(nrow(Xd), 1, pnorm(td)), f_true = td, f_genes = fd,
       X_test = Xv, Z_test = base$Z_val, Y_test = rbinom(nrow(Xv), 1, pnorm(tv)),
       beta_true = STUDY_BETA, true_genes = colnames(Xd)[sig],
       blocks = NULL, Sigma_X = NULL, X_gaussian = FALSE)
}

#-------------------------------------------------------------------------------
# Metrics
#-------------------------------------------------------------------------------
auc_safe <- function(y, p) {
  if (length(unique(y)) < 2 || all(is.na(p))) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(y, p, quiet = TRUE, levels = c(0, 1), direction = "<")))
}

selection_metrics <- function(sel, dat) {
  tg <- dat$true_genes
  out <- list(n_selected = length(sel),
              sensitivity = mean(tg %in% sel),
              FDR = if (length(sel)) mean(!sel %in% tg) else 0)
  if (!is.null(dat$blocks)) {            # block-level view
    tb <- unique(dat$blocks[tg]); sb <- unique(dat$blocks[sel])
    out$block_sensitivity <- mean(tb %in% sb)
    out$block_FDR <- if (length(sb)) mean(!sb %in% tb) else 0
  }
  if (!is.null(dat$proxy_genes)) {       # null genes correlated with a covariate (Study 14)
    out$n_proxy_selected <- sum(sel %in% dat$proxy_genes)
    out$proxy_share_of_selected <- if (length(sel)) mean(sel %in% dat$proxy_genes) else 0
  }
  out
}

#-------------------------------------------------------------------------------
# Average marginal effects (AMEs) of the covariates on the probability scale,
# a scale every method shares: for a binary covariate, the mean over the rows
# dat$ame_rows of P(Y=1 | z_j = 1) - P(Y=1 | z_j = 0); for a continuous one, of
# P(Y=1 | z_j + 1) - P(Y=1 | z_j) (one s.d.). True values come with the data
# (dat$ame_truth); each method returns posterior draws or a point estimate.
#-------------------------------------------------------------------------------
ame_contrasts <- function(Z, rows) {
  out <- lapply(colnames(Z), function(v) {
    z <- Z[rows, , drop = FALSE]; z1 <- z; z0 <- z
    if (all(Z[, v] %in% c(0, 1))) { z1[, v] <- 1; z0[, v] <- 0 } else { z1[, v] <- z[, v] + 1 }
    list(z1 = z1, z0 = z0)
  })
  setNames(out, colnames(Z))
}
ame_summary <- function(D, truth) {                     # D: draws x covariates
  lo <- apply(D, 2, quantile, 0.025); hi <- apply(D, 2, quantile, 0.975); m <- colMeans(D)
  tr <- truth[colnames(D)]
  data.frame(covariate = colnames(D), truth = tr, mean = m, ci_lower = lo, ci_upper = hi,
             bias = m - tr, covers = (tr >= lo & tr <= hi), excludes_zero = (lo > 0 | hi < 0), row.names = NULL)
}
#' True AMEs for a model whose latent index is f_genes + eta_Z(Z)
ame_truth <- function(f_genes, Z, rows, eta_Z) {
  cs <- ame_contrasts(Z, rows)
  vapply(cs, function(c) mean(pnorm(f_genes[rows] + eta_Z(c$z1)) - pnorm(f_genes[rows] + eta_Z(c$z0))), numeric(1))
}

prediction_metrics <- function(dat, prob) {
  if (is.null(prob)) return(list(AUC = NA_real_, Brier = NA_real_))
  list(AUC = auc_safe(dat$Y_test, prob), Brier = mean((dat$Y_test - prob)^2))
}

#' Per-coefficient summaries from posterior draws (rows = draws)
beta_metrics <- function(draws, truth) {
  truth <- truth[colnames(draws)]
  lo <- apply(draws, 2, quantile, 0.025); hi <- apply(draws, 2, quantile, 0.975)
  m <- colMeans(draws)
  data.frame(covariate = colnames(draws), truth = truth, mean = m, sd = apply(draws, 2, sd),
             ci_lower = lo, ci_upper = hi, ci_width = hi - lo,
             bias = m - truth, covers = (truth >= lo & truth <= hi),
             excludes_zero = (lo > 0 | hi < 0),
             sign_prob = vapply(seq_along(m), function(j)
               if (is.na(truth[j])) NA_real_ else mean(sign(draws[, j]) == sign(truth[j])), 1),
             row.names = NULL)
}

#-------------------------------------------------------------------------------
# Random-number streams. Data sets are generated after set.seed(seed); a method
# that called set.seed(seed) again would replay the data's random numbers
# (e.g. "knockoffs" identical to X). Every method therefore uses its own stream.
#-------------------------------------------------------------------------------
method_seed <- function(seed, method) {
  k <- match(method, c("BART", "enet", "knockoff", "pipeline", "stan4bart"))
  if (is.na(k)) stop("unknown method for method_seed(): ", method)
  (seed %% 1e8) + 1e9 + k * 1e8
}

#-------------------------------------------------------------------------------
# Method wrappers. Each returns a list with: method, selected, selection
# metrics, prediction metrics, (beta), minutes. `mc` holds MCMC settings:
#   n_burn, n_iter, n_thin (spBART / BART), s4b_iter, s4b_warmup, s4b_chains
#-------------------------------------------------------------------------------
default_mc <- list(n_burn = 2000, n_iter = 1000, n_thin = 5,
                   s4b_iter = 1500, s4b_warmup = 1000, s4b_chains = 4)

m_spbart <- function(dat, mc = default_mc, sigma2_beta = 10, partial = FALSE,
                     beta_init = NULL, seed = NULL, keep_fit = FALSE) {
  t0 <- Sys.time()
  X <- dat$X; Z <- dat$Z; Xt <- dat$X_test; Zt <- dat$Z_test
  if (partial) {                         # variant: age and BMI go to the trees
    X  <- cbind(X,  Z[,  c("age", "bmi")]);  Z  <- Z[,  c("sex", "race")]
    Xt <- cbind(Xt, Zt[, c("age", "bmi")]); Zt <- Zt[, c("sex", "race")]
  }
  fit <- fit_spBART(X = X, Z = Z, Y = dat$Y, X_test = Xt, Z_test = Zt,
                    n_burn = mc$n_burn, n_iter = mc$n_iter, n_thin = mc$n_thin,
                    rho = ncol(X), sigma_beta = sqrt(sigma2_beta),
                    beta_init = beta_init, seed = seed, verbose = FALSE)
  pip <- colMeans(fit$varcount > 0)[colnames(dat$X)]
  sel <- fdr_gene_selection(pip, 0.05)$selected_genes
  # Variance shares of the linear and tree components
  Zb <- fit$beta_draws %*% t(Z)
  VZ <- apply(Zb, 1, var); Vf <- apply(fit$f_train_draws, 1, var)
  pZ <- colMeans(pnorm(fit$offset + fit$beta_draws %*% t(Zt)))
  pf <- colMeans(pnorm(fit$offset + fit$f_test_draws))
  out <- c(list(method = if (partial) "spBART_partialZ" else "spBART", selected = sel, pip = pip),
           selection_metrics(sel, dat), prediction_metrics(dat, fit$prob.test.mean),
           list(beta = beta_metrics(fit$beta_draws, dat$beta_true),
                var_share_Z = mean(VZ / (VZ + Vf + 1)), var_share_f = mean(Vf / (VZ + Vf + 1)),
                AUC_Z_only = auc_safe(dat$Y_test, pZ), AUC_f_only = auc_safe(dat$Y_test, pf),
                Brier_Z_only = mean((dat$Y_test - pZ)^2), Brier_f_only = mean((dat$Y_test - pf)^2),
                f_RMSE = if (!is.null(dat$f_genes))
                  sqrt(mean(((fit$f_hat - mean(fit$f_hat)) - (dat$f_genes - mean(dat$f_genes)))^2)) else NA_real_,
                mu_RMSE = sqrt(mean((fit$mu_hat - dat$f_true)^2)),
                minutes = as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  if (!is.null(dat$ame_rows) && !partial) {
    cs <- ame_contrasts(Z, dat$ame_rows)
    eta_f <- fit$offset + fit$f_train_draws[, dat$ame_rows, drop = FALSE]
    D <- sapply(cs, function(c) rowMeans(pnorm(eta_f + fit$beta_draws %*% t(c$z1)) -
                                         pnorm(eta_f + fit$beta_draws %*% t(c$z0))))
    out$ame <- ame_summary(D, dat$ame_truth)
  }
  if (keep_fit) out$fit <- fit
  out
}

m_bart <- function(dat, mc = default_mc, seed = NULL, keep_fit = FALSE) {
  t0 <- Sys.time()
  if (!is.null(seed)) set.seed(method_seed(seed, "BART"))
  XZ <- cbind(dat$X, dat$Z); XZt <- cbind(dat$X_test, dat$Z_test)
  fit <- pbart(x.train = XZ, y.train = dat$Y, x.test = XZt, ntree = 200,
               ndpost = mc$n_iter, nskip = mc$n_burn, keepevery = mc$n_thin,
               sparse = TRUE, a = 0.5, b = 1, rho = ncol(XZ), printevery = 1e6)
  colnames(fit$varcount) <- colnames(XZ)
  pip_all <- colMeans(fit$varcount > 0)
  pip <- pip_all[colnames(dat$X)]
  sel <- fdr_gene_selection(pip, 0.05)$selected_genes
  cov_pip <- pip_all[colnames(dat$Z)]
  cov_sel <- fdr_gene_selection(cov_pip, 0.05)$selected_genes   # Bayesian FDR rule on the covariate PIPs
  out <- c(list(method = "BART", selected = sel, pip = pip), selection_metrics(sel, dat),
           prediction_metrics(dat, fit$prob.test.mean),
           list(covariate_pip = cov_pip,
                covariate_selected = setNames(colnames(dat$Z) %in% cov_sel, colnames(dat$Z)),
                mu_RMSE = sqrt(mean((colMeans(fit$yhat.train) - dat$f_true)^2)),
                minutes = as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  if (!is.null(dat$ame_rows)) {
    cs <- ame_contrasts(dat$Z, dat$ame_rows); Xr <- dat$X[dat$ame_rows, , drop = FALSE]
    pr <- function(z) predict(fit, newdata = cbind(Xr, z), mc.cores = 1)$prob.test   # draws x rows
    D <- sapply(cs, function(c) rowMeans(pr(c$z1) - pr(c$z0)))
    out$ame <- ame_summary(D, dat$ame_truth)
  }
  if (keep_fit) out$fit <- fit
  out
}

m_stan4bart <- function(dat, mc = default_mc, seed = 1) {
  t0 <- Sys.time()
  suppressPackageStartupMessages(library(stan4bart))
  zn <- colnames(dat$Z)
  # stan4bart silently drops predictors with non-syntactic names (e.g. the real
  # gene "RP11-192H23.4"), so genes get safe internal names x1..xp
  safe <- paste0("x", seq_len(ncol(dat$X)))
  Xs <- dat$X; Xst <- dat$X_test
  colnames(Xs) <- colnames(Xst) <- safe
  d  <- data.frame(y = dat$Y, Xs, dat$Z)
  dt <- data.frame(Xst, dat$Z_test)
  # stan4bart reads its formula from the unevaluated call, so the formula is
  # inserted into the call as language (a formula object built with
  # as.formula() fails with "object '.' not found")
  fml <- str2lang(paste("y ~", paste(zn, collapse = " + "),
                        "+ bart(. -", paste(zn, collapse = " - "), ")"))
  fit <- eval(bquote(stan4bart(.(fml), data = d, chains = .(mc$s4b_chains), iter = .(mc$s4b_iter),
                               warmup = .(mc$s4b_warmup), cores = 1, seed = .(method_seed(seed, "stan4bart")), verbose = -1,
                               bart_args = list(n.trees = 200, keepTrees = TRUE))))
  vc <- extract(fit, type = "varcount")            # p x draws, rows named x1..xp
  pip <- setNames(rowMeans(vc > 0)[safe], colnames(dat$X))
  if (anyNA(pip)) stop("stan4bart dropped predictors: ", paste(colnames(dat$X)[is.na(pip)], collapse = ", "))
  sel <- fdr_gene_selection(pip, 0.05)$selected_genes
  fe <- t(extract(fit, type = "fixef"))            # draws x J
  prob <- rowMeans(predict(fit, newdata = dt, type = "ev"))
  out <- c(list(method = "stan4bart", selected = sel, pip = pip), selection_metrics(sel, dat),
           prediction_metrics(dat, prob),
           list(beta = beta_metrics(fe[, zn, drop = FALSE], dat$beta_true),
                minutes = as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  if (!is.null(dat$ame_rows)) {
    cs <- ame_contrasts(dat$Z, dat$ame_rows); Xr <- Xs[dat$ame_rows, , drop = FALSE]
    pr <- function(z) predict(fit, newdata = data.frame(Xr, z), type = "ev")       # rows x draws
    D <- sapply(cs, function(c) colMeans(pr(c$z1) - pr(c$z0)))
    out$ame <- ame_summary(D, dat$ame_truth)
  }
  out
}

m_enet <- function(dat, seed = 1) {
  t0 <- Sys.time()
  set.seed(method_seed(seed, "enet"))
  J <- ncol(dat$Z)
  x <- cbind(dat$Z, dat$X)
  cv <- glmnet::cv.glmnet(x, dat$Y, family = "binomial", alpha = 0.5,
                          penalty.factor = c(rep(0, J), rep(1, ncol(dat$X))), nfolds = 10)
  b <- as.numeric(coef(cv, s = "lambda.1se"))[-1]
  names(b) <- colnames(x)
  sel <- names(b)[-(1:J)][b[-(1:J)] != 0]
  prob <- as.numeric(predict(cv, newx = cbind(dat$Z_test, dat$X_test), s = "lambda.1se", type = "response"))
  out <- c(list(method = "enet", selected = sel), selection_metrics(sel, dat), prediction_metrics(dat, prob),
           list(minutes = as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  if (!is.null(dat$ame_rows)) {
    cs <- ame_contrasts(dat$Z, dat$ame_rows); Xr <- dat$X[dat$ame_rows, , drop = FALSE]
    pr <- function(z) as.numeric(predict(cv, newx = cbind(z, Xr), s = "lambda.1se", type = "response"))
    est <- vapply(cs, function(c) mean(pr(c$z1) - pr(c$z0)), numeric(1))
    out$ame <- data.frame(covariate = names(est), truth = dat$ame_truth[names(est)], mean = est,
                          ci_lower = NA_real_, ci_upper = NA_real_, bias = est - dat$ame_truth[names(est)],
                          covers = NA, excludes_zero = NA, row.names = NULL)
  }
  out
}

#-------------------------------------------------------------------------------
# Model-X knockoffs (Candes et al. 2018), implemented here so the study does not
# depend on the knockoff package (whose Rdsdp dependency does not compile on
# every cluster). Reproduces knockoff::create.gaussian(X, mu, Sigma, method = "equi"):
#   s = min(1, 2 * lambda_min(cor(Sigma))), D = s * diag(Sigma),
#   X_k | X ~ N(X - (X - mu) Sigma^{-1} D, 2D - D Sigma^{-1} D).
# Gaussian X with known covariance: mu = 0 and the true Sigma (Sigma = I gives
# exact independent N(0, 1) knockoffs). Real X: knockoff::create.second_order()
# rule -- sample mean and covariance, or corpcor::cov.shrink() if the sample
# covariance is not positive definite (e.g. p > n).
#-------------------------------------------------------------------------------
equi_knockoffs <- function(X, mu, Sigma) {
  n <- nrow(X); p <- ncol(X)
  d <- sqrt(diag(Sigma))
  Xs <- sweep(sweep(X, 2, mu), 2, d, "/")
  e <- eigen(cov2cor(Sigma), symmetric = TRUE)
  lam <- pmax(e$values, 1e-10)
  s <- min(1, 2 * min(lam)) * (1 - 1e-6)
  U <- e$vectors
  Rinv_s <- U %*% (t(U) * (s / lam))                  # s R^{-1}
  Vhalf  <- U %*% (t(U) * sqrt(pmax(2 * s - s^2 / lam, 0)))
  Xk <- Xs - Xs %*% Rinv_s + matrix(rnorm(n * p), n, p) %*% Vhalf
  sweep(sweep(Xk, 2, d, "*"), 2, mu, "+")
}

second_order_knockoffs <- function(X) {
  Sigma <- cov(X)
  if (min(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values) <= 1e-10) {
    Sigma <- matrix(as.numeric(corpcor::cov.shrink(X, verbose = FALSE)), ncol(X))
  }
  equi_knockoffs(X, colMeans(X), Sigma)
}

#' Knockoff selection threshold (knockoff::knockoff.threshold); offset 0 or 1
knockoff_threshold <- function(W, fdr = 0.05, offset = 0) {
  ts <- sort(c(0, abs(W)))
  ratio <- vapply(ts, function(t) (offset + sum(W <= -t)) / max(1, sum(W >= t)), numeric(1))
  ok <- which(ratio <= fdr)
  if (length(ok)) ts[ok[1]] else Inf
}

#' Model-X knockoffs. Statistic: lasso coefficient difference with Z unpenalized.
#' offset = 0 (the "knockoff" rather than "knockoff+" threshold): with only 10
#' signals, knockoff+ at FDR 0.05 needs at least 20 selections and so can never
#' select anything here.
m_knockoff <- function(dat, fdr = 0.05, seed = 1) {
  t0 <- Sys.time()
  set.seed(method_seed(seed, "knockoff"))
  Z <- dat$Z; J <- ncol(Z); X <- dat$X; p <- ncol(X)
  Xk <- if (dat$X_gaussian && is.null(dat$Sigma_X)) {
    matrix(rnorm(length(X)), nrow(X), p)               # Sigma = I: exact knockoffs
  } else if (dat$X_gaussian) {
    equi_knockoffs(X, rep(0, p), dat$Sigma_X)
  } else {
    second_order_knockoffs(X)
  }
  cv <- glmnet::cv.glmnet(cbind(Z, X, Xk), dat$Y, family = "binomial", alpha = 1,
                          penalty.factor = c(rep(0, J), rep(1, 2 * p)), nfolds = 10)
  b <- as.numeric(coef(cv, s = "lambda.min"))[-1][-(1:J)]
  W <- abs(b[1:p]) - abs(b[(p + 1):(2 * p)])
  thr <- knockoff_threshold(W, fdr, offset = 0)
  sel <- colnames(X)[W >= thr]
  c(list(method = "knockoff", selected = sel), selection_metrics(sel, dat),
    prediction_metrics(dat, NULL),
    list(minutes = as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

#' Run a set of methods on one data set; failures are recorded, not fatal
run_methods <- function(dat, methods, mc = default_mc, seed = 1) {
  out <- list()
  for (m in methods) {
    out[[m]] <- tryCatch(switch(m,
      spBART          = m_spbart(dat, mc),
      spBART_partialZ = m_spbart(dat, mc, partial = TRUE),
      BART            = m_bart(dat, mc, seed = seed),
      stan4bart       = m_stan4bart(dat, mc, seed = seed),
      enet            = m_enet(dat, seed = seed),
      knockoff        = m_knockoff(dat, seed = seed)),
      error = function(e) list(method = m, error = conditionMessage(e)))
  }
  out
}

#-------------------------------------------------------------------------------
# Probit + GMM screening, as in the application (Sections 2.1 and 4.1): marginal probit
# per gene adjusted for age and sex, p < 0.05, then the GMM cluster (2-4
# components, BIC) with the largest mean |coefficient|.
#-------------------------------------------------------------------------------
screen_genes <- function(X, Y, Zadj) {
  D <- cbind(1, Zadj)
  res <- t(vapply(seq_len(ncol(X)), function(j) {
    fit <- tryCatch(suppressWarnings(glm.fit(cbind(D, X[, j]), Y,
                    family = binomial(link = "probit"), control = list(maxit = 100))),
                    error = function(e) NULL)
    if (is.null(fit) || !fit$converged || fit$rank < ncol(D) + 1) return(c(0, 1))
    V <- tryCatch(chol2inv(qr.R(fit$qr)), error = function(e) NULL)
    if (is.null(V)) return(c(0, 1))
    k <- ncol(D) + 1
    b <- fit$coefficients[k]; se <- sqrt(V[k, k])
    c(b, 2 * pnorm(-abs(b / se)))
  }, numeric(2)))
  sig <- which(res[, 2] < 0.05)
  if (length(sig) <= 4) return(colnames(X)[sig])
  ab <- abs(res[sig, 1])
  gm <- tryCatch(mclust::Mclust(ab, G = 2:4, modelNames = "V", verbose = FALSE), error = function(e) NULL)
  if (is.null(gm)) return(colnames(X)[sig])
  means <- tapply(ab, gm$classification, mean)
  colnames(X)[sig][gm$classification == as.integer(names(which.max(means)))]
}

#-------------------------------------------------------------------------------
# 5-fold CV-union pipeline (Section 3.3) for spBART or BART on a given gene set,
# then the final fit on the union genes scored on the test set.
#-------------------------------------------------------------------------------
cv_union_pipeline <- function(method, dat, genes, mc = default_mc, K = 5, seed = 1) {
  t0 <- Sys.time()
  set.seed(method_seed(seed, "pipeline"))
  folds <- caret::createFolds(factor(dat$Y), k = K, list = TRUE, returnTrain = FALSE)
  fold_sel <- vector("list", K)
  for (k in seq_len(K)) {
    tr <- setdiff(seq_along(dat$Y), folds[[k]])
    Xk <- dat$X[tr, genes, drop = FALSE]; Zk <- dat$Z[tr, , drop = FALSE]; Yk <- dat$Y[tr]
    if (method == "spBART") {
      fit <- fit_spBART(Xk, Zk, Yk, n_burn = mc$n_burn, n_iter = mc$n_iter, n_thin = mc$n_thin,
                        rho = ncol(Xk), verbose = FALSE)
      pip <- colMeans(fit$varcount > 0)[genes]
    } else {
      XZ <- cbind(Xk, Zk)
      fit <- pbart(x.train = XZ, y.train = Yk, ntree = 200, ndpost = mc$n_iter, nskip = mc$n_burn,
                   keepevery = mc$n_thin, sparse = TRUE, a = 0.5, b = 1, rho = ncol(XZ),
                   printevery = 1e6)
      pip <- colMeans(fit$varcount[, seq_along(genes), drop = FALSE] > 0)
      names(pip) <- genes
    }
    fold_sel[[k]] <- fdr_gene_selection(pip, 0.05)$selected_genes
  }
  union <- unique(unlist(fold_sel))
  final <- if (length(union) == 0) {
    list(prob = NULL)
  } else {
    sub <- dat
    sub$X <- dat$X[, union, drop = FALSE]; sub$X_test <- dat$X_test[, union, drop = FALSE]
    r <- if (method == "spBART") m_spbart(sub, mc) else m_bart(sub, mc, seed = seed + 1)
    r
  }
  out <- c(list(method = paste0(method, "_pipeline"), selected = union, fold_selected = fold_sel),
           selection_metrics(union, dat),
           if (length(union)) list(AUC = final$AUC, Brier = final$Brier) else list(AUC = NA, Brier = NA),
           list(minutes = as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  if (method == "spBART" && length(union)) out$beta <- final$beta
  out
}

#-------------------------------------------------------------------------------
# Array-task bookkeeping shared by the drivers
#-------------------------------------------------------------------------------
task_id <- function() {
  id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  if (is.na(id)) {
    a <- commandArgs(trailingOnly = TRUE)
    id <- if (length(a)) as.integer(a[1]) else 1L
  }
  id
}

study_root <- function(script_dir, study) {
  root <- "results"
  file.path(script_dir, root, study)
}

result_file <- function(dir, label, rep) {
  d <- file.path(dir, label)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  file.path(d, sprintf("rep_%03d.RDS", rep))
}

#-------------------------------------------------------------------------------
# Study 14: realistic adjustment set + null genes that track a strong confounder
#   8 covariates (continuous: age, BMI, education, activity; binary: sex, race,
#   smoking, family history) with linear effects; a strong age effect beta_age;
#   10 signal genes (DGM 2 trees, independent of the covariates) and n_proxy NULL
#   genes correlated with age (corr rho_proxy) that have no effect given age.
#   age_form = "linear" (spBART correctly specified) or "ushape" (control arm:
#   age acts as beta_age * (age^2 - 1)/sqrt(2), which a linear term cannot carry).
#-------------------------------------------------------------------------------
CONF_BETA <- c(age = NA, sex = 0.30, race = -0.30, bmi = 0.25, smoking = 0.20,
               education = -0.15, activity = -0.10, famhx = 0.40)
conf_scenario <- function(label, n, p, beta_age, rho_proxy, n_proxy = 30, age_form = "linear")
  list(label = label, n = n, p = p, beta_age = beta_age, rho_proxy = rho_proxy,
       n_proxy = n_proxy, age_form = age_form, type = "confounder")

gen_confounder_data <- function(sc, seed, n_test = 1000, n_ame = 300) {
  beta <- CONF_BETA; beta["age"] <- sc$beta_age
  proxy_idx <- 10 + seq_len(sc$n_proxy)
  eta_Z <- function(Z) {
    lin <- as.numeric(Z[, names(beta)[-1], drop = FALSE] %*% beta[-1])
    lin + if (sc$age_form == "linear") sc$beta_age * Z[, "age"] else sc$beta_age * (Z[, "age"]^2 - 1) / sqrt(2)
  }
  one <- function(n) {
    Z <- cbind(age = rnorm(n), sex = rbinom(n, 1, 0.5), race = rbinom(n, 1, 0.3), bmi = rnorm(n),
               smoking = rbinom(n, 1, 0.25), education = rnorm(n), activity = rnorm(n), famhx = rbinom(n, 1, 0.1))
    X <- matrix(rnorm(n * sc$p), n, sc$p)
    X[, proxy_idx] <- sc$rho_proxy * Z[, "age"] + sqrt(1 - sc$rho_proxy^2) * X[, proxy_idx]
    colnames(X) <- paste0("gene_", seq_len(sc$p))
    fg <- f_trees(X[, 1:10], qnorm(STUDY_CUT_Q))
    ft <- fg + eta_Z(Z)
    list(X = X, Z = Z, Y = rbinom(n, 1, pnorm(ft)), f_true = ft, f_genes = fg)
  }
  set.seed(seed); tr <- one(sc$n)
  set.seed(seed + 5e8); te <- one(n_test)
  rows <- seq_len(min(n_ame, sc$n))
  btrue <- beta; if (sc$age_form != "linear") btrue["age"] <- NA
  list(X = tr$X, Z = tr$Z, Y = tr$Y, f_true = tr$f_true, f_genes = tr$f_genes,
       X_test = te$X, Z_test = te$Z, Y_test = te$Y, beta_true = btrue,
       true_genes = paste0("gene_", 1:10), proxy_genes = paste0("gene_", proxy_idx),
       blocks = NULL, Sigma_X = NULL, X_gaussian = TRUE,
       ame_rows = rows, ame_truth = ame_truth(tr$f_genes, tr$Z, rows, eta_Z))
}

#' Attach true AMEs to a linear-covariate data set from gen_study_data() (Study 15)
add_ame_truth <- function(dat, n_ame = 300) {
  rows <- seq_len(min(n_ame, nrow(dat$X)))
  b <- dat$beta_true
  dat$ame_rows <- rows
  dat$ame_truth <- ame_truth(dat$f_genes, dat$Z, rows, function(Z) as.numeric(Z[, names(b)] %*% b))
  dat
}
