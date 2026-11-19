################################################################################
# 11_study_multichain.R
#
# One long chain versus many short chains for the full-data fit (screening and
# cross-validation are held out of this study, so that only the MCMC differs).
#   S0     1 chain : 2,000 burn-in + 5,000 iterations thinned by 5 (as in the paper)
#   S10eq 10 chains:   200 burn-in +   500 thinned by 5 (same 7,000 total iterations)
#   S10   10 chains: 1,000 burn-in + 2,000 thinned by 2 (3,000 iterations each)
#   S20   20 chains, S50 50 chains: as S10
# Every chain has its own seed; spBART chains start from beta ~ N(0, 4I).
# Methods: spBART and BART (pbart). Draws are pooled across chains.
# Scenarios (100 replicates): DGM1 and DGM2 at n = 1500, p in {500, 2000};
# covariates correlated with signal genes (c = 0.9, p = 500); application-sized
# semi-synthetic data (115 screened real genes, n = 500).
# Output per strategy: selection, prediction, beta, agreement with S0
# (PIP differences, Jaccard), and across-chain diagnostics (rank-normalized
# split-Rhat, bulk/tail ESS for beta, log-likelihood, number of splits, f at 5
# test points; between-chain PIP s.d.; share of stuck chains).
#
# Array: task = (scenario, replicate) = 6 x 100 = 600 tasks; chains run in
# parallel on --cpus-per-task cores.
################################################################################

script_dir <- tryCatch({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", file_arg[1]), fixed = TRUE)))
}, error = function(e) getwd())
spbart_root <- dirname(script_dir)
source(file.path(spbart_root, "R", "spBART_sampler.R")); load_spBART_BART()
source(file.path(script_dir, "utils", "data_generation.R"))
source(file.path(script_dir, "utils", "study_functions.R"))
suppressPackageStartupMessages({ library(parallel); library(posterior) })

study  <- "11_multichain"
n_reps <- 100
ncores <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
strategies <- list(S0    = list(C = 1,  nb = 2000, ni = 1000, nt = 5),
     S10eq = list(C = 10, nb = 200,  ni = 100,  nt = 5),
     S10   = list(C = 10, nb = 1000, ni = 1000, nt = 2),
     S20   = list(C = 20, nb = 1000, ni = 1000, nt = 2),
     S50   = list(C = 50, nb = 1000, ni = 1000, nt = 2))
nn <- 1500
pp <- c(500, 2000)
scenarios <- list(
  list(label = "DGM1_n1500_p500",  sc = scenario("DGM1_n1500_p500",  nn, pp[1], type = "dgm1")),
  list(label = "DGM1_n1500_p2000", sc = scenario("DGM1_n1500_p2000", nn, pp[2], type = "dgm1")),
  list(label = "DGM2_n1500_p500",  sc = scenario("DGM2_n1500_p500",  nn, pp[1])),
  list(label = "DGM2_n1500_p2000", sc = scenario("DGM2_n1500_p2000", nn, pp[2])),
  list(label = "corrZ_c0.9_n1500_p500", sc = scenario("corrZ_c0.9_n1500_p500", nn, pp[1], corrZ = 0.9)),
  list(label = "semisyn_app_n500_p115", semisyn = TRUE))

#-------------------------------------------------------------------------------
# One chain of each method; returns only what the summaries need
#-------------------------------------------------------------------------------
chain_spbart <- function(dat, s, seed) {
  set.seed(seed)
  b0 <- rnorm(ncol(dat$Z), 0, 2)                       # dispersed start, beta ~ N(0, 4I)
  fit <- fit_spBART(dat$X, dat$Z, dat$Y, X_test = dat$X_test, Z_test = dat$Z_test,
                    n_burn = s$nb, n_iter = s$ni, n_thin = s$nt, rho = ncol(dat$X),
                    beta_init = b0, seed = seed + 1, verbose = FALSE)
  list(pip = colMeans(fit$varcount > 0), beta = fit$beta_draws, loglik = fit$loglik,
       nsplits = rowSums(fit$varcount), f5 = fit$f_test_draws[, 1:5, drop = FALSE],
       prob = fit$prob.test.mean)
}
chain_bart <- function(dat, s, seed) {
  set.seed(seed)
  XZ <- cbind(dat$X, dat$Z)
  fit <- pbart(x.train = XZ, y.train = dat$Y, x.test = cbind(dat$X_test, dat$Z_test),
               ntree = 200, ndpost = s$ni, nskip = s$nb, keepevery = s$nt,
               sparse = TRUE, a = 0.5, b = 1, rho = ncol(XZ), printevery = 1e6)
  colnames(fit$varcount) <- colnames(XZ)
  eta <- fit$yhat.train
  ll <- rowSums(pnorm(eta[, dat$Y == 1, drop = FALSE], log.p = TRUE)) +
        rowSums(pnorm(-eta[, dat$Y == 0, drop = FALSE], log.p = TRUE))
  list(pip = colMeans(fit$varcount > 0)[colnames(dat$X)], beta = NULL, loglik = ll,
       nsplits = rowSums(fit$varcount), f5 = fit$yhat.test[, 1:5, drop = FALSE],
       prob = fit$prob.test.mean)
}

diag_one <- function(M) {                              # M: iterations x chains
  M <- as.matrix(M)
  c(rhat = posterior::rhat(M), ess_bulk = posterior::ess_bulk(M), ess_tail = posterior::ess_tail(M))
}

summarize_chains <- function(chains, dat, ref_pip = NULL) {
  pipM <- sapply(chains, `[[`, "pip")
  pipM <- matrix(pipM, ncol = length(chains), dimnames = list(colnames(dat$X), NULL))
  pip <- rowMeans(pipM)
  sel <- fdr_gene_selection(pip, 0.05)$selected_genes
  prob <- rowMeans(sapply(chains, `[[`, "prob"))
  tg <- dat$true_genes
  out <- c(list(selected = sel), selection_metrics(sel, dat), prediction_metrics(dat, prob))
  C <- length(chains)
  out$pip_between_sd_signal <- if (C > 1) mean(apply(pipM[tg, , drop = FALSE], 1, sd)) else NA
  out$pip_between_sd_null   <- if (C > 1) mean(apply(pipM[!rownames(pipM) %in% tg, , drop = FALSE], 1, sd)) else NA
  out$stuck_share <- mean(pipM[tg, , drop = FALSE] < 0.5)
  qd <- list(loglik = sapply(chains, `[[`, "loglik"), nsplits = sapply(chains, `[[`, "nsplits"))
  for (k in 1:5) qd[[paste0("f_test", k)]] <- sapply(chains, function(ch) ch$f5[, k])
  if (!is.null(chains[[1]]$beta)) {
    for (j in colnames(chains[[1]]$beta)) qd[[paste0("beta_", j)]] <- sapply(chains, function(ch) ch$beta[, j])
    out$beta <- beta_metrics(do.call(rbind, lapply(chains, `[[`, "beta")), dat$beta_true)
  }
  out$diagnostics <- do.call(rbind, lapply(qd, diag_one))
  bt <- grep("^beta_", rownames(out$diagnostics))
  out$max_rhat_beta <- if (length(bt)) max(out$diagnostics[bt, "rhat"]) else NA
  if (!is.null(ref_pip)) {
    out$pip_mean_absdiff_vs_S0 <- mean(abs(pip - ref_pip$pip))
    out$pip_max_absdiff_vs_S0  <- max(abs(pip - ref_pip$pip))
    u <- union(sel, ref_pip$sel)
    out$jaccard_vs_S0 <- if (length(u)) length(intersect(sel, ref_pip$sel)) / length(u) else 1
  }
  list(summary = out, pip = pip, sel = sel)
}

run_strategy <- function(fun, dat, s, seed0) {
  ch <- mclapply(seq_len(s$C), function(c) fun(dat, s, seed0 + 1000 * c),
                 mc.cores = min(ncores, s$C), mc.preschedule = FALSE)
  bad <- vapply(ch, function(x) inherits(x, "try-error") || is.null(x$pip), logical(1))
  if (any(bad)) stop("chain failed: ", paste(unique(sapply(ch[bad], as.character)), collapse = " | "))
  ch
}

#-------------------------------------------------------------------------------
id <- task_id()
todo <- list(c(ceiling(id / n_reps), id - (ceiling(id / n_reps) - 1) * n_reps))
out_dir <- study_root(script_dir, study)

for (sr in todo) {
  s_idx <- sr[1]; rep <- sr[2]; S <- scenarios[[s_idx]]
  f_out <- result_file(out_dir, S$label, rep)
  if (file.exists(f_out)) { cat("Done already:", f_out, "\n"); next }
  seed <- 11e6 + s_idx * 1e4 + rep
  dat <- if (is.null(S$semisyn)) gen_study_data(S$sc, seed, n_test = 1000) else
    gen_semisyn_data(readRDS(file.path(script_dir, "data_semisynthetic", "semisynthetic_base.RDS")),
                     "screened", "random", "trees", seed)
  res <- list(); traces <- list()
  for (m in c("spBART", "BART")) {
    fun <- if (m == "spBART") chain_spbart else chain_bart
    ref <- NULL
    for (st in names(strategies)) {
      t0 <- Sys.time()
      r <- tryCatch({
        ch <- run_strategy(fun, dat, strategies[[st]], seed + 100 * match(st, names(strategies)))
        sm <- summarize_chains(ch, dat, ref)
        if (st == "S0") ref <- list(pip = sm$pip, sel = sm$sel)
        if (rep == 1 && st == "S10") {                 # traceplot data for the Supplement
          traces[[m]] <- lapply(ch[1:min(4, length(ch))], function(x) x[c("beta", "loglik", "nsplits")])
        }
        sm$summary
      }, error = function(e) list(error = conditionMessage(e)))
      r$minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      r$strategy <- strategies[[st]]
      res[[m]][[st]] <- r
      cat(sprintf("%s %s %s: %.1f min\n", S$label, m, st, r$minutes))
    }
  }
  saveRDS(list(study = study, scenario = S$label, rep = rep, seed = seed, ncores = ncores,
               results = res, traces = traces), f_out)
  cat("Saved", f_out, "\n")
}
