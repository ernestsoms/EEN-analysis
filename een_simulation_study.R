# Simulation study for the Enhanced Elastic Net (EEN) estimator.
# Runs 500 Monte Carlo replications in four phases (50, 150, 300, 500)
# under a smooth two-bump DGP with AR(1) correlation, producing all
# figures and tables reported in the manuscript.

library(MASS)
library(glmnet)
library(Matrix)
library(parallel)
library(doParallel)
library(foreach)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(dplyr)
library(tidyr)

# ---- Constants -------------------------------------------------------------
LAMBDA_RIDGE_DEFAULT <- 1
EPSILON_WEIGHT <- 1e-6
EPSILON_ALASSO <- 1e-6
ACTIVE_TOL <- 0.05
BUMP_HALF_WIDTH <- 2.5

BUMP_WIDTH_TARGET <- 6L

METHOD_COLORS <- c(
  "EEN" = "#377eb8",
  "LASSO" = "#984ea3",
  "Elastic Net" = "#4daf4a",
  "Adaptive LASSO" = "#e41a1c"
)

# ---- Data generation --------------------------------------------------------
generate_data <- function(n = 150, p = 500, s = 2, sigma_bump = 1.05,
                          margin = 20L, jitter = 30L, bump_gap = 0,
                          rho = 0.7, sigma = 1, seed = 2026) {
  set.seed(seed)
  idx <- seq_len(p)
  separation <- 2 * BUMP_HALF_WIDTH + bump_gap
  block_mid <- p / 2 + runif(1, -jitter, jitter)
  block_mid <- pmax(margin + BUMP_HALF_WIDTH + 1,
                    pmin(p - margin - BUMP_HALF_WIDTH - 1, block_mid))
  centres <- block_mid + c(-separation / 2, separation / 2)
  centres <- floor(centres) + 0.5
  amp_signs <- c(+1, -1)
  amp_mags <- runif(s, 2.0, 3.0)
  amps <- amp_signs * amp_mags
  beta.true <- rep(0.0, p)
  for (k in seq_len(s)) {
    beta.true <- beta.true + amps[k] * exp(-0.5 * ((idx - centres[k]) / sigma_bump)^2)
  }
  active <- which(abs(beta.true) > ACTIVE_TOL)
  Sigma <- rho^abs(outer(idx, idx, "-"))
  X <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  X <- scale(X)
  y <- X %*% beta.true + rnorm(n, sd = sigma)
  list(X = X, y = y, beta.true = beta.true, active = active,
       bump_centres = centres, bump_widths = rep(sigma_bump, s), bump_amps = amps)
}

# Builds Omega = D^T D directly as a sparse banded matrix (Matrix::dgCMatrix)
# without ever forming a dense p x p matrix. Omega is banded with bandwidth 5
# (Omega[i,j] = 0 whenever |i - j| > 2), which build_omega_neighbors() below
# exploits so that each coordinate update only touches O(1) neighbors.
create_omega <- function(p) {
  if (p >= 3) {
    D <- Matrix::bandSparse(n = p - 2, m = p, k = c(0, 1, 2),
                            diagonals = list(rep(1, p - 2),
                                             rep(-2, p - 2),
                                             rep(1, p - 2)))
    Matrix::crossprod(D)
  } else {
    Matrix::Diagonal(p)
  }
}

# Precomputed once per Omega (not per coordinate, not per sweep): for each j,
# the handful of k != j with Omega[j, k] != 0, plus the corresponding values
# and the diagonal Omega[j, j]. Because Omega is banded, every entry of
# nbr_idx has length <= 4 regardless of p, which is what makes the c_j term
# in the coordinate update O(1) rather than O(p).
build_omega_neighbors <- function(Omega) {
  # Omega = crossprod(D) is returned by Matrix as a symmetric-storage
  # "dsCMatrix", which keeps only one triangle of nonzeros in @i/@p/@x.
  # Coercing straight to "CsparseMatrix" preserves that symmetric storage
  # class (dsCMatrix *is* a CsparseMatrix) and silently drops half of every
  # off-diagonal row/column's entries. Coercing through "generalMatrix"
  # first forces the full, non-symmetric dgCMatrix representation so every
  # neighbor actually gets read.
  Omega <- methods::as(methods::as(Omega, "generalMatrix"), "CsparseMatrix")
  p <- ncol(Omega)
  omega_diag <- Matrix::diag(Omega)
  nbr_idx <- vector("list", p)
  nbr_val <- vector("list", p)
  cptr <- Omega@p; irow <- Omega@i; xval <- Omega@x
  for (j in seq_len(p)) {
    rng <- (cptr[j] + 1):cptr[j + 1]
    if (length(rng) == 0) {
      nbr_idx[[j]] <- integer(0); nbr_val[[j]] <- numeric(0); next
    }
    rows <- irow[rng] + 1L  # CsparseMatrix @i is 0-indexed row within column j
    vals <- xval[rng]
    keep <- rows != j
    nbr_idx[[j]] <- rows[keep]
    nbr_val[[j]] <- vals[keep]
  }
  list(diag = omega_diag, idx = nbr_idx, val = nbr_val)
}

# ---- Ridge initialization ---------------------------------------------------
# Used to build adaptive weights for both EEN and Adaptive LASSO. Falls back
# to a marginal-correlation estimate if the ridge system is singular.
# Purely deterministic: no random draws, so consolidating call sites cannot
# change any simulated value.
ridge_init <- function(X, y, lambda_ridge) {
  tryCatch(
    as.vector(solve(t(X) %*% X + lambda_ridge * diag(ncol(X))) %*% t(X) %*% y),
    error = function(e) as.vector(cor(X, y))
  )
}

# ---- EEN estimator -----------------------------------------------------------
# Coordinate descent with two efficiency properties that the manuscript
# claims but the earlier implementation did not have:
#   1. The residual r = y - X %*% beta is maintained incrementally. X %*% beta
#      is computed once (at beta = 0, so r = y) and thereafter each coordinate
#      update only does O(n) vector arithmetic on r -- it never recomputes a
#      full matrix-vector product inside the j-loop.
#   2. c_j = sum_{k != j} Omega[j, k] * beta[k] is evaluated only over the
#      O(1) banded neighbors of j (via omega_nbrs, built once outside this
#      function and reused across calls), not over all p - 1 coefficients.
# Together this makes a full sweep O(np) as claimed in the manuscript, rather
# than O(np^2). v_j = ||x_j||^2 / n is also precomputed once (matching line 3
# of Algorithm 1), not recomputed inside the j-loop.
een <- function(X, y, lambda1, lambda2, Omega, gamma = 1,
                lambda_ridge = LAMBDA_RIDGE_DEFAULT, maxit = 1000,
                weights = NULL, omega_nbrs = NULL) {
  n <- nrow(X); p <- ncol(X)
  if (is.null(weights)) {
    beta0 <- ridge_init(X, y, lambda_ridge)
    weights <- 1 / (abs(beta0)^gamma + EPSILON_WEIGHT)
  }
  if (is.null(omega_nbrs)) omega_nbrs <- build_omega_neighbors(Omega)
  
  beta <- rep(0, p)
  r <- as.vector(y - X %*% beta)     # one O(np) product; beta = 0 here so r = y
  col_ss <- colSums(X^2) / n         # v_j, precomputed once
  
  for (iter in seq_len(maxit)) {
    beta_old <- beta
    for (j in seq_len(p)) {
      xj <- X[, j]
      beta_j_old <- beta[j]
      
      # r_{-j} recovered from the maintained residual in O(n), no full
      # X %*% beta recomputation.
      r_partial <- r + xj * beta_j_old
      zj <- sum(xj * r_partial) / n
      
      idx <- omega_nbrs$idx[[j]]
      cj <- if (length(idx)) sum(omega_nbrs$val[[j]] * beta[idx]) else 0
      
      adj <- zj - lambda2 * cj
      soft <- sign(adj) * max(abs(adj) - lambda1 * weights[j], 0)
      beta_j_new <- soft / (col_ss[j] + lambda2 * omega_nbrs$diag[j])
      
      r <- r_partial - xj * beta_j_new   # incremental update, O(n)
      beta[j] <- beta_j_new
    }
    if (max(abs(beta - beta_old)) < 1e-6) break
  }
  beta
}

# Note: candidate fits during cross-validation use maxit = 500 (sufficient to
# rank candidates), while the final fit in run_replication() uses een()'s own
# default of maxit = 1000 for a fully converged estimate.
cv_een <- function(X, y, Omega, lambda1_grid = 10^seq(-3, 1, length = 6),
                   lambda2_grid = 10^seq(-3, 1, length = 6),
                   nfolds = 5, gamma = 1, lambda_ridge = LAMBDA_RIDGE_DEFAULT,
                   maxit = 500, omega_nbrs = NULL) {
  n <- nrow(X); p <- ncol(X)
  if (is.null(omega_nbrs)) omega_nbrs <- build_omega_neighbors(Omega)
  foldid <- sample(rep(seq_len(nfolds), length.out = n))
  cv_err <- matrix(NA, length(lambda1_grid), length(lambda2_grid))
  for (i in seq_along(lambda1_grid)) {
    for (j in seq_along(lambda2_grid)) {
      fold_err <- numeric(nfolds)
      for (fold in seq_len(nfolds)) {
        tr_X <- X[foldid != fold, , drop = FALSE]
        tr_y <- y[foldid != fold]
        va_X <- X[foldid == fold, , drop = FALSE]
        va_y <- y[foldid == fold]
        ridge_tr <- ridge_init(tr_X, tr_y, lambda_ridge)
        w_tr <- 1 / (abs(ridge_tr)^gamma + EPSILON_WEIGHT)
        b <- tryCatch(
          een(tr_X, tr_y, lambda1_grid[i], lambda2_grid[j], Omega,
              gamma, lambda_ridge, maxit, weights = w_tr, omega_nbrs = omega_nbrs),
          error = function(e) rep(0, ncol(X))
        )
        fold_err[fold] <- mean((va_y - va_X %*% b)^2)
      }
      cv_err[i, j] <- mean(fold_err)
    }
  }
  idx <- which(cv_err == min(cv_err, na.rm = TRUE), arr.ind = TRUE)[1, ]
  list(lambda1 = lambda1_grid[idx[1]], lambda2 = lambda2_grid[idx[2]], cv_errors = cv_err)
}

# ---- Metrics -----------------------------------------------------------------
metrics <- function(beta.hat, beta.true, tol = 1e-4) {
  detected <- abs(beta.hat) > tol
  true_pos <- abs(beta.true) > ACTIVE_TOL
  tp <- sum(detected & true_pos)
  tn <- sum(!detected & !true_pos)
  fp <- sum(detected & !true_pos)
  fn <- sum(!detected & true_pos)
  TPR <- if (sum(true_pos) > 0) tp / sum(true_pos) else NA
  FDR <- if ((tp + fp) > 0) fp / (tp + fp) else 0
  L2 <- sqrt(sum((beta.hat - beta.true)^2))
  MCC_denom <- sqrt((tp+fp) * (tp+fn) * (tn+fp) * (tn+fn))
  MCC <- if (!is.nan(MCC_denom) && MCC_denom > 0) (tp * tn - fp * fn) / MCC_denom else 0
  c(L2 = L2, TPR = TPR, FDR = FDR, MCC = MCC)
}

# ---- Single replication --------------------------------------------------------
run_replication <- function(seed, n = 150, p = 500, s = 2, sigma_bump = 1.05,
                            margin = 20L, jitter = 30L, bump_gap = 0,
                            rho = 0.7, sigma = 1) {
  data <- generate_data(n, p, s, sigma_bump, margin, jitter, bump_gap,
                        rho = rho, sigma = sigma, seed = seed)
  X <- data$X; y <- data$y; beta.true <- data$beta.true
  n_active <- length(data$active)
  Omega <- create_omega(p)
  omega_nbrs <- build_omega_neighbors(Omega)
  
  # EEN
  cv_res <- cv_een(X, y, Omega, gamma = 1, lambda_ridge = LAMBDA_RIDGE_DEFAULT,
                   omega_nbrs = omega_nbrs)
  b_een <- een(X, y, cv_res$lambda1, cv_res$lambda2, Omega,
               gamma = 1, lambda_ridge = LAMBDA_RIDGE_DEFAULT, omega_nbrs = omega_nbrs)
  m_een <- metrics(b_een, beta.true)
  
  # LASSO
  cv_lasso <- cv.glmnet(X, y, alpha = 1)
  b_lasso <- as.vector(coef(cv_lasso, s = "lambda.min"))[-1]
  m_lasso <- metrics(b_lasso, beta.true)
  
  # Elastic Net (grid search over alpha)
  alpha_seq <- seq(0.1, 0.9, by = 0.1)
  best_cvm <- Inf; best_alpha <- 0.5; best_cv_enet <- NULL
  for (a in alpha_seq) {
    cv_a <- cv.glmnet(X, y, alpha = a)
    if (min(cv_a$cvm) < best_cvm) {
      best_cvm <- min(cv_a$cvm)
      best_alpha <- a
      best_cv_enet <- cv_a
    }
  }
  b_enet <- as.vector(coef(best_cv_enet, s = "lambda.min"))[-1]
  m_enet <- metrics(b_enet, beta.true)
  
  # Adaptive LASSO (same ridge initialization as EEN)
  beta_ridge_init <- ridge_init(X, y, LAMBDA_RIDGE_DEFAULT)
  w_alasso <- 1 / (abs(beta_ridge_init)^1 + EPSILON_ALASSO)
  cv_alasso <- cv.glmnet(X, y, alpha = 1, penalty.factor = w_alasso)
  b_alasso <- as.vector(coef(cv_alasso, s = "lambda.min"))[-1]
  m_alasso <- metrics(b_alasso, beta.true)
  
  data.frame(
    Method = c("EEN", "LASSO", "Elastic Net", "Adaptive LASSO"),
    L2 = c(m_een["L2"], m_lasso["L2"], m_enet["L2"], m_alasso["L2"]),
    TPR = c(m_een["TPR"], m_lasso["TPR"], m_enet["TPR"], m_alasso["TPR"]),
    FDR = c(m_een["FDR"], m_lasso["FDR"], m_enet["FDR"], m_alasso["FDR"]),
    MCC = c(m_een["MCC"], m_lasso["MCC"], m_enet["MCC"], m_alasso["MCC"]),
    n_active = n_active,
    alpha_used = c(NA, NA, best_alpha, NA),
    Replication = seed,
    stringsAsFactors = FALSE
  )
}

# ---- Checkpointing -------------------------------------------------------------
checkpoint_dir <- "./checkpoints"
if (!dir.exists(checkpoint_dir)) dir.create(checkpoint_dir, recursive = TRUE)

get_checkpoint_path <- function(...) {
  parts <- paste0(c(...), collapse = "_")
  file.path(checkpoint_dir, paste0(parts, ".rds"))
}

load_checkpoint <- function(path) {
  if (file.exists(path)) readRDS(path) else NULL
}

save_checkpoint <- function(data, path) saveRDS(data, path)

# ---- Simulation driver -----------------------------------------------------------
run_simulation <- function(n_replications = 50, n = 150, p = 500, s = 2,
                           sigma_bump = 1.05, margin = 20L, jitter = 30L,
                           bump_gap = 0, rho = 0.7, sigma = 1,
                           use_parallel = TRUE, n_cores = NULL, verbose = TRUE,
                           checkpoint = TRUE, chunk_size = 25) {
  if (verbose) {
    cat("Simulation: n =", n, ", p =", p, ", s =", s,
        ", sigma_bump =", sigma_bump, ", bump_gap =", bump_gap, "\n")
    cat("  rho =", rho, ", sigma =", sigma, ", replications =", n_replications, "\n")
  }
  if (is.null(n_cores)) n_cores <- max(1, detectCores() - 1)
  if (verbose) cat("  cores:", n_cores, ", parallel:", use_parallel,
                   ", checkpointing:", checkpoint, "\n")
  
  cp_file <- NULL; completed_reps <- integer(0); existing_results <- NULL
  if (checkpoint) {
    cp_file <- get_checkpoint_path(
      paste0("sim_n", n_replications), paste0("s", s), paste0("rho", rho), paste0("sig", sigma)
    )
    cp_data <- load_checkpoint(cp_file)
    if (!is.null(cp_data)) {
      existing_results <- cp_data$results
      completed_reps <- cp_data$completed_reps
      if (verbose) cat("  checkpoint loaded:", length(completed_reps), "reps done\n")
    }
  }
  need_reps <- setdiff(seq_len(n_replications), completed_reps)
  if (length(need_reps) == 0) {
    if (verbose) cat("  all replications already done\n")
    return(existing_results)
  }
  if (verbose) cat("  remaining:", length(need_reps), "replications\n")
  
  if (use_parallel && length(need_reps) > 1) {
    all_new <- data.frame()
    chunks <- split(need_reps, ceiling(seq_along(need_reps) / chunk_size))
    for (ci in seq_along(chunks)) {
      chunk_reps <- chunks[[ci]]
      if (verbose) cat("  chunk", ci, "of", length(chunks), ": reps",
                       paste(range(chunk_reps), collapse = "-"), "\n")
      cl <- makeCluster(n_cores)
      registerDoParallel(cl)
      clusterExport(cl, c("generate_data", "create_omega", "build_omega_neighbors",
                          "een", "cv_een", "metrics",
                          "run_replication", "ridge_init", "LAMBDA_RIDGE_DEFAULT",
                          "EPSILON_WEIGHT", "EPSILON_ALASSO", "ACTIVE_TOL",
                          "BUMP_HALF_WIDTH", "n", "p", "s", "sigma_bump", "margin",
                          "jitter", "bump_gap", "rho", "sigma"), envir = environment())
      clusterEvalQ(cl, { library(MASS); library(glmnet); library(Matrix) })
      chunk_res <- foreach(rep = chunk_reps,
                           .packages = c("MASS", "glmnet", "Matrix"),
                           .combine = rbind) %dopar% {
                             seed <- 2026 + rep
                             rr <- run_replication(seed, n, p, s, sigma_bump, margin, jitter, bump_gap, rho, sigma)
                             rr$Replication <- rep
                             rr
                           }
      stopCluster(cl)
      all_new <- rbind(all_new, chunk_res)
      if (checkpoint) {
        combined <- rbind(existing_results, all_new)
        done <- c(completed_reps, chunk_reps)
        save_checkpoint(list(results = combined, completed_reps = done), cp_file)
        if (verbose) cat("    checkpoint saved after chunk", ci, "\n")
      }
    }
    all_results <- rbind(existing_results, all_new)
  } else {
    all_results <- existing_results
    for (rep in need_reps) {
      if (verbose && rep %% 10 == 0) cat(sprintf("  replication %d / %d\n", rep, n_replications))
      seed <- 2026 + rep
      rr <- run_replication(seed, n, p, s, sigma_bump, margin, jitter, bump_gap, rho, sigma)
      rr$Replication <- rep
      all_results <- rbind(all_results, rr)
      if (checkpoint) {
        done <- c(completed_reps, rep)
        save_checkpoint(list(results = all_results, completed_reps = done), cp_file)
      }
    }
  }
  all_results$Replication <- as.integer(all_results$Replication)
  if (checkpoint) {
    save_checkpoint(list(results = all_results, completed_reps = seq_len(n_replications)), cp_file)
    if (verbose) cat("  final checkpoint saved\n")
  }
  all_results
}

# ---- Summaries and significance ---------------------------------------------------
summarize_results <- function(results) {
  methods <- unique(results$Method)
  do.call(rbind, lapply(methods, function(m) {
    r <- results[results$Method == m, ]
    data.frame(
      Method = m, L2_mean = round(mean(r$L2, na.rm = TRUE), 3),
      L2_sd = round(sd(r$L2, na.rm = TRUE), 3),
      TPR_mean = round(mean(r$TPR, na.rm = TRUE), 3),
      TPR_sd = round(sd(r$TPR, na.rm = TRUE), 3),
      FDR_mean = round(mean(r$FDR, na.rm = TRUE), 3),
      FDR_sd = round(sd(r$FDR, na.rm = TRUE), 3),
      MCC_mean = round(mean(r$MCC, na.rm = TRUE), 3),
      MCC_sd = round(sd(r$MCC, na.rm = TRUE), 3),
      n_active_mean = round(mean(r$n_active, na.rm = TRUE), 1)
    )
  }))
}

test_significance <- function(results, baseline = "EEN") {
  comparisons <- data.frame()
  other_methods <- setdiff(unique(results$Method), baseline)
  base_data <- results[results$Method == baseline, c("Replication", "L2")]
  for (meth in other_methods) {
    meth_data <- results[results$Method == meth, c("Replication", "L2")]
    paired <- merge(base_data, meth_data, by = "Replication", suffixes = c(".base", ".meth"))
    if (nrow(paired) < 2) next
    t_res <- t.test(paired$L2.base, paired$L2.meth, paired = TRUE)
    comparisons <- rbind(comparisons, data.frame(
      Comparison = paste(baseline, "vs", meth), Mean_diff = round(t_res$estimate, 4),
      CI_lower = round(t_res$conf.int[1], 4), CI_upper = round(t_res$conf.int[2], 4),
      P_value = format.pval(t_res$p.value, digits = 3),
      Significant = ifelse(t_res$p.value < 0.05, "Yes", "No"), n_pairs = nrow(paired)
    ))
  }
  comparisons
}

# ---- Plots -----------------------------------------------------------------------
# One routine handles all four per-method metric plots (L2, FDR, TPR, MCC):
# a violin plot with mean +/- 1 SD overlaid, methods ordered by their mean
# value. Sort direction and an optional reference line are the only things
# that differ across metrics. Layer order (violin, then reference line, then
# summary markers/labels) matches the original per-metric plotting functions.
plot_metric_comparison <- function(results, metric, label,
                                   direction = c("lower", "higher"),
                                   hline = NULL) {
  direction <- match.arg(direction)
  summ <- results %>%
    group_by(Method) %>%
    summarise(mean_val = mean(.data[[metric]], na.rm = TRUE),
              sd_val = sd(.data[[metric]], na.rm = TRUE),
              .groups = "drop") %>%
    mutate(ymin = mean_val - sd_val, ymax = mean_val + sd_val)
  ord <- if (direction == "lower") order(summ$mean_val) else order(-summ$mean_val)
  summ$Method <- factor(summ$Method, levels = summ$Method[ord])
  results$Method <- factor(results$Method, levels = levels(summ$Method))
  
  p <- ggplot(results, aes(x = Method, y = .data[[metric]], fill = Method)) +
    geom_violin(alpha = 0.5, trim = TRUE, colour = NA)
  
  if (!is.null(hline)) {
    p <- p + geom_hline(yintercept = hline, linetype = "dashed", colour = "red", alpha = 0.5)
  }
  
  p +
    geom_pointrange(data = summ, aes(x = Method, y = mean_val, ymin = ymin, ymax = ymax),
                    inherit.aes = FALSE, size = 0.7, linewidth = 1, colour = "black") +
    geom_text(data = summ, aes(x = Method, y = ymax, label = round(mean_val, 3)),
              inherit.aes = FALSE, vjust = -0.8, size = 3.5) +
    theme_minimal(base_size = 11) +
    labs(x = NULL, y = paste0(label, " (mean \u00b1 1 SD)")) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1)) +
    scale_fill_manual(values = METHOD_COLORS)
}

plot_tradeoff <- function(results) {
  method_means <- aggregate(cbind(L2, FDR) ~ Method, data = results, FUN = mean)
  ggplot(method_means, aes(x = FDR, y = L2, colour = Method, label = Method)) +
    geom_point(size = 4) + geom_text(vjust = -1, hjust = 0.5, size = 3.5) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 11) +
    labs(x = "Mean FDR", y = "Mean L2 Error") +
    scale_colour_manual(values = METHOD_COLORS) +
    theme(legend.position = "none", plot.margin = margin(t = 20, r = 10, b = 10, l = 10))
}

plot_true_beta <- function(p = 500, s = 2, sigma_bump = 1.05, bump_gap = 0, seed = 2026) {
  ex <- generate_data(n = 10, p = p, s = s, sigma_bump = sigma_bump,
                      bump_gap = bump_gap, seed = seed)
  df <- data.frame(j = seq_len(p), beta = ex$beta.true, is_active = abs(ex$beta.true) > ACTIVE_TOL)
  df$sign_col <- ifelse(ex$beta.true > ACTIVE_TOL, "Positive active",
                        ifelse(ex$beta.true < -ACTIVE_TOL, "Negative active", "Inactive"))
  fine_x <- seq(1, p, length.out = 20 * p); fine_y <- rep(0, length(fine_x))
  for (k in seq_len(s)) {
    fine_y <- fine_y + ex$bump_amps[k] * exp(-0.5 * ((fine_x - ex$bump_centres[k]) / sigma_bump)^2)
  }
  fine_df <- data.frame(j = fine_x, beta = fine_y)
  ggplot(df, aes(j, beta)) +
    geom_line(data = fine_df, aes(j, beta), colour = "grey40", linewidth = 0.7) +
    geom_point(data = df[df$is_active, ], aes(j, beta, colour = sign_col), size = 2.2, alpha = 0.9) +
    scale_colour_manual(values = c("Positive active" = "steelblue", "Negative active" = "firebrick")) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    coord_cartesian(xlim = c(max(1, min(ex$bump_centres) - 25), min(p, max(ex$bump_centres) + 25))) +
    theme_bw(base_size = 11) +
    labs(x = "Predictor index j", y = expression(beta[j]^{"true"})) +
    theme(legend.position = "top")
}

plot_convergence <- function(all_summaries) {
  conv_df <- all_summaries %>%
    mutate(Replications = as.numeric(gsub(".*\\((\\d+) reps\\)", "\\1", Phase)))
  ggplot(conv_df, aes(x = Replications, y = L2_mean, colour = Method)) +
    geom_line(linewidth = 1.2) + geom_point(size = 3) +
    geom_errorbar(aes(ymin = L2_mean - L2_sd / sqrt(Replications),
                      ymax = L2_mean + L2_sd / sqrt(Replications)), width = 30, alpha = 0.5) +
    theme_minimal(base_size = 11) +
    labs(x = "Number of Replications", y = "L2 Error (mean \u00b1 SE)") +
    scale_colour_manual(values = METHOD_COLORS) +
    scale_x_continuous(breaks = c(50, 150, 300, 500))
}

save_results <- function(df, filename, dir = "simulation_results") {
  path <- file.path(dir, filename)
  write.csv(df, path, row.names = FALSE)
  cat("Saved:", path, "\n")
}

# ---- Main execution ----------------------------------------------------------------
available_cores <- detectCores()
cat("Detected", available_cores, "cores.\n")
if (!dir.exists("simulation_results")) dir.create("simulation_results")

S <- 2
SIGMA_BUMP <- 1.05
BUMP_GAP <- 0

cat(sprintf("DGP: %d adjacent bumps (bump 1 positive, bump 2 negative); expected active width ~%d predictors per bump (%d total).\n",
            S, BUMP_WIDTH_TARGET, S * BUMP_WIDTH_TARGET))

p_true_beta <- plot_true_beta(p = 500, s = S, sigma_bump = SIGMA_BUMP, bump_gap = BUMP_GAP, seed = 2026)
ggsave("simulation_results/true_beta_diagnostic.pdf", p_true_beta, width = 10, height = 4, dpi = 300)

all_summaries <- data.frame()

run_phase <- function(n_reps, phase_label, s, sigma_bump, bump_gap, available_cores) {
  cat("\n", phase_label, "\n", sep = "")
  t <- system.time({
    res <- run_simulation(n_replications = n_reps, s = s, sigma_bump = sigma_bump,
                          bump_gap = bump_gap, use_parallel = TRUE,
                          n_cores = max(1, available_cores - 1), verbose = TRUE, checkpoint = TRUE)
  })
  sm <- summarize_results(res); sm$Phase <- phase_label
  sig <- test_significance(res); sig$Phase <- phase_label
  tag <- gsub("[^0-9]", "", regmatches(phase_label, regexpr("\\d+", phase_label)))
  save_results(res, paste0("phase", tag, "_results_", n_reps, ".csv"))
  save_results(sm, paste0("phase", tag, "_summary_", n_reps, ".csv"))
  save_results(sig, paste0("phase", tag, "_significance_", n_reps, ".csv"))
  plots <- list(
    l2    = plot_metric_comparison(res, "L2",  "L2 Error", "lower"),
    fdr   = plot_metric_comparison(res, "FDR", "FDR",      "lower", hline = 0.5),
    tpr   = plot_metric_comparison(res, "TPR", "TPR",      "higher"),
    mcc   = plot_metric_comparison(res, "MCC", "MCC",      "higher"),
    trade = plot_tradeoff(res)
  )
  for (nm in names(plots))
    ggsave(paste0("simulation_results/phase", tag, "_", nm, ".pdf"), plots[[nm]], width = 8, height = 5, dpi = 300)
  cat(phase_label, "summary:\n"); print(sm)
  list(results = res, summary = sm, sig = sig, time = t, plots = plots)
}

ph1 <- run_phase(50, "Phase 1 (50 reps)", S, SIGMA_BUMP, BUMP_GAP, available_cores)
ph2 <- run_phase(150, "Phase 2 (150 reps)", S, SIGMA_BUMP, BUMP_GAP, available_cores)
ph3 <- run_phase(300, "Phase 3 (300 reps)", S, SIGMA_BUMP, BUMP_GAP, available_cores)
ph4 <- run_phase(500, "Phase 4 (500 reps)", S, SIGMA_BUMP, BUMP_GAP, available_cores)

all_summaries <- rbind(ph1$summary, ph2$summary, ph3$summary, ph4$summary)

comparison_table <- all_summaries %>%
  select(Phase, Method, L2_mean, L2_sd, TPR_mean, TPR_sd, FDR_mean, FDR_sd, MCC_mean, MCC_sd) %>%
  arrange(Phase, Method)
save_results(comparison_table, "cross_phase_comparison.csv")

p_conv <- plot_convergence(all_summaries)
ggsave("simulation_results/convergence_plot.pdf", p_conv, width = 10, height = 6, dpi = 300)

final_table <- ph4$summary %>%
  transmute(Method, L2 = paste0(L2_mean, " (", L2_sd, ")"),
            TPR = paste0(TPR_mean, " (", TPR_sd, ")"),
            FDR = paste0(FDR_mean, " (", FDR_sd, ")"),
            MCC = paste0(MCC_mean, " (", MCC_sd, ")"), n_active_mean)
save_results(final_table, "final_paper_table.csv")

cat("\nSimulation complete. All files saved to simulation_results/\n")