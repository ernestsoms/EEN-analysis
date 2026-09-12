# Empirical application: corn NIR spectra.
# Data: pcv::corn -- 80 corn samples, 700 wavelengths (1100-2498 nm, 2 nm steps).
# Response: moisture content.

library(pcv)
library(glmnet)
library(Matrix)
library(ggplot2)
library(dplyr)
library(tidyr)
library(parallel)
library(doParallel)
library(foreach)

SEED_BASE <- 2026
LAMBDA_RIDGE_DEFAULT <- 1
EPSILON_WEIGHT <- 1e-6
EPSILON_ALASSO <- 1e-6
GAMMA <- 1

METHOD_COLORS <- c(
  "EEN" = "#377eb8",
  "LASSO" = "#984ea3",
  "Elastic Net" = "#4daf4a",
  "Adaptive LASSO" = "#e41a1c"
)

# ---- Load data -------------------------------------------------------------
data(corn, package = "pcv")
X_raw <- as.matrix(corn$spectra)
y_raw <- as.numeric(corn$moisture)

wavelengths <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", colnames(X_raw))))
if (any(is.na(wavelengths)) || length(wavelengths) != ncol(X_raw)) {
  wavelengths <- seq(1100, 2498, by = 2)
  if (length(wavelengths) != ncol(X_raw)) wavelengths <- seq_len(ncol(X_raw))
}

n <- nrow(X_raw); p <- ncol(X_raw)
cat(sprintf("corn data: n = %d, p = %d wavelengths\n", n, p))

# Builds Omega = D^T D directly as a sparse banded matrix (Matrix::dgCMatrix)
# without ever forming a dense p x p matrix -- important here since p = 700.
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

# Precomputed once (Omega is fixed for the whole script): for each j, the
# handful of k != j with Omega[j, k] != 0, their values, and the diagonal
# Omega[j, j]. Bandedness means every list entry has length <= 4 regardless
# of p, which is what makes the c_j term in een()'s coordinate update O(1)
# instead of the O(p) dense sum used previously.
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
    rows <- irow[rng] + 1L
    vals <- xval[rng]
    keep <- rows != j
    nbr_idx[[j]] <- rows[keep]
    nbr_val[[j]] <- vals[keep]
  }
  list(diag = omega_diag, idx = nbr_idx, val = nbr_val)
}

Omega <- create_omega(p)
omega_nbrs <- build_omega_neighbors(Omega)

# ---- Ridge initialization --------------------------------------------------
# Used to build adaptive weights for EEN and Adaptive LASSO. Falls back to a
# marginal-correlation estimate if the ridge system is singular. Deterministic:
# no random draws, so consolidating call sites cannot change any output.
ridge_init <- function(X, y, lambda_ridge) {
  tryCatch(
    as.vector(solve(t(X) %*% X + lambda_ridge * diag(ncol(X))) %*% t(X) %*% y),
    error = function(e) as.vector(cor(X, y))
  )
}

# ---- EEN estimator ---------------------------------------------------------
# As in the simulation script: the residual r = y - X %*% beta is maintained
# incrementally (O(n) per coordinate, no full matrix-vector recompute inside
# the j-loop), and c_j is summed only over the O(1) banded neighbors of j via
# the precomputed omega_nbrs list, not over all p - 1 coefficients. This
# matters most here since p = 700; a full sweep is now O(np) rather than
# O(np^2).
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

      r_partial <- r + xj * beta_j_old
      zj <- sum(xj * r_partial) / n

      idx <- omega_nbrs$idx[[j]]
      cj <- if (length(idx)) sum(omega_nbrs$val[[j]] * beta[idx]) else 0

      adj <- zj - lambda2 * cj
      soft <- sign(adj) * max(abs(adj) - lambda1 * weights[j], 0)
      beta_j_new <- soft / (col_ss[j] + lambda2 * omega_nbrs$diag[j])

      r <- r_partial - xj * beta_j_new
      beta[j] <- beta_j_new
    }
    if (max(abs(beta - beta_old)) < 1e-6) break
  }
  beta
}

# Candidate fits during cross-validation use maxit = 500 (enough to rank
# candidates); the final fit in fit_predict_een() uses maxit = 1000 for a
# fully converged estimate.
cv_een <- function(X, y, Omega, lambda1_grid = 10^seq(-3, 1, length = 6),
                   lambda2_grid = 10^seq(-3, 1, length = 6),
                   nfolds = 5, gamma = 1, lambda_ridge = LAMBDA_RIDGE_DEFAULT,
                   maxit = 500, omega_nbrs = NULL) {
  n <- nrow(X)
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

roughness <- function(beta, Omega) as.numeric(t(beta) %*% Omega %*% beta)

standardize_fold <- function(X_train, X_test, y_train, y_test) {
  x_center <- colMeans(X_train)
  x_scale <- apply(X_train, 2, sd); x_scale[x_scale == 0] <- 1
  y_mean <- mean(y_train)
  list(X_train = scale(X_train, center = x_center, scale = x_scale),
       X_test = scale(X_test, center = x_center, scale = x_scale),
       y_train = y_train - y_mean, y_test = y_test - y_mean, y_mean = y_mean)
}

# ---- Configuration ---------------------------------------------------------
EEN_GRID_LEN <- 6; EEN_INNER_NFOLDS <- 5; EEN_MAXIT_GRID <- 500; EEN_MAXIT_FINAL <- 1000
ENET_ALPHA_SEQ <- seq(0.1, 0.9, by = 0.1)
N_REPEATS <- 20; K_FOLDS <- 10
EEN_LAMBDA1_LOG_RANGE <- c(-4, 1)
EEN_LAMBDA2_LOG_RANGE <- c(-4, 3)

# ---- Per-fold model wrappers ----------------------------------------------
fit_predict_een <- function(X_tr, y_tr, X_te, Omega, grid_len = EEN_GRID_LEN,
                            inner_folds = EEN_INNER_NFOLDS, maxit_grid = EEN_MAXIT_GRID,
                            maxit_final = EEN_MAXIT_FINAL, omega_nbrs = NULL) {
  if (is.null(omega_nbrs)) omega_nbrs <- build_omega_neighbors(Omega)
  l1_grid <- 10^seq(EEN_LAMBDA1_LOG_RANGE[1], EEN_LAMBDA1_LOG_RANGE[2], length = grid_len)
  l2_grid <- 10^seq(EEN_LAMBDA2_LOG_RANGE[1], EEN_LAMBDA2_LOG_RANGE[2], length = grid_len)
  cv <- cv_een(X_tr, y_tr, Omega, lambda1_grid = l1_grid, lambda2_grid = l2_grid,
               nfolds = inner_folds, gamma = GAMMA, lambda_ridge = LAMBDA_RIDGE_DEFAULT,
               maxit = maxit_grid, omega_nbrs = omega_nbrs)
  b <- een(X_tr, y_tr, cv$lambda1, cv$lambda2, Omega, gamma = GAMMA,
           lambda_ridge = LAMBDA_RIDGE_DEFAULT, maxit = maxit_final, omega_nbrs = omega_nbrs)
  list(pred = as.vector(X_te %*% b), beta = b,
       lambda1 = cv$lambda1, lambda2 = cv$lambda2,
       lambda1_at_boundary = cv$lambda1 %in% c(min(l1_grid), max(l1_grid)),
       lambda2_at_boundary = cv$lambda2 %in% c(min(l2_grid), max(l2_grid)))
}

fit_predict_lasso <- function(X_tr, y_tr, X_te) {
  cvfit <- cv.glmnet(X_tr, y_tr, alpha = 1, intercept = FALSE)
  b <- as.vector(coef(cvfit, s = "lambda.min"))[-1]
  list(pred = as.vector(predict(cvfit, newx = X_te, s = "lambda.min")), beta = b)
}

fit_predict_enet <- function(X_tr, y_tr, X_te, alpha_seq = ENET_ALPHA_SEQ) {
  best_cvm <- Inf; best_fit <- NULL; best_alpha <- 0.5
  for (a in alpha_seq) {
    cv_a <- cv.glmnet(X_tr, y_tr, alpha = a, intercept = FALSE)
    if (min(cv_a$cvm) < best_cvm) { best_cvm <- min(cv_a$cvm); best_fit <- cv_a; best_alpha <- a }
  }
  b <- as.vector(coef(best_fit, s = "lambda.min"))[-1]
  list(pred = as.vector(predict(best_fit, newx = X_te, s = "lambda.min")), beta = b, alpha = best_alpha)
}

fit_predict_alasso <- function(X_tr, y_tr, X_te, gamma = GAMMA) {
  ridge_init_vec <- ridge_init(X_tr, y_tr, LAMBDA_RIDGE_DEFAULT)
  w <- 1 / (abs(ridge_init_vec)^gamma + EPSILON_ALASSO)
  cvfit <- cv.glmnet(X_tr, y_tr, alpha = 1, penalty.factor = w, intercept = FALSE)
  b <- as.vector(coef(cvfit, s = "lambda.min"))[-1]
  list(pred = as.vector(predict(cvfit, newx = X_te, s = "lambda.min")), beta = b)
}

# ---- Cross-validation driver ----------------------------------------------
run_cv_comparison <- function(X, y, Omega, n_repeats = N_REPEATS, k_folds = K_FOLDS,
                              use_parallel = TRUE, n_cores = NULL) {
  n <- nrow(X)
  if (is.null(n_cores)) n_cores <- max(1, detectCores() - 1)
  if (use_parallel) {
    cl <- makeCluster(n_cores); registerDoParallel(cl)
    on.exit(stopCluster(cl), add = TRUE)
    clusterExport(cl, c("een", "cv_een", "ridge_init", "roughness", "standardize_fold",
                        "build_omega_neighbors", "omega_nbrs",
                        "fit_predict_een", "fit_predict_lasso", "fit_predict_enet",
                        "fit_predict_alasso", "LAMBDA_RIDGE_DEFAULT", "EPSILON_WEIGHT",
                        "EPSILON_ALASSO", "GAMMA", "SEED_BASE", "EEN_GRID_LEN",
                        "EEN_INNER_NFOLDS", "EEN_MAXIT_GRID", "EEN_MAXIT_FINAL",
                        "ENET_ALPHA_SEQ", "EEN_LAMBDA1_LOG_RANGE", "EEN_LAMBDA2_LOG_RANGE"),
                  envir = globalenv())
    clusterEvalQ(cl, { library(glmnet); library(Matrix) })
  }
  `%op%` <- if (use_parallel) `%dopar%` else `%do%`
  foreach(rep_id = seq_len(n_repeats), .packages = "glmnet", .combine = rbind) %op% {
    set.seed(SEED_BASE + rep_id)
    foldid <- sample(rep(seq_len(k_folds), length.out = n))
    rep_out <- data.frame()
    for (fold in seq_len(k_folds)) {
      tr <- which(foldid != fold); te <- which(foldid == fold)
      sd_fold <- standardize_fold(X[tr, , drop = FALSE], X[te, , drop = FALSE], y[tr], y[te])
      r_een <- fit_predict_een(sd_fold$X_train, sd_fold$y_train, sd_fold$X_test, Omega,
                                omega_nbrs = omega_nbrs)
      r_lasso <- fit_predict_lasso(sd_fold$X_train, sd_fold$y_train, sd_fold$X_test)
      r_enet <- fit_predict_enet(sd_fold$X_train, sd_fold$y_train, sd_fold$X_test)
      r_alasso <- fit_predict_alasso(sd_fold$X_train, sd_fold$y_train, sd_fold$X_test)
      sq_err <- function(r) mean((sd_fold$y_test - r$pred)^2)
      rep_out <- rbind(rep_out, data.frame(
        Repetition = rep_id, Fold = fold,
        Method = c("EEN", "LASSO", "Elastic Net", "Adaptive LASSO"),
        MSE = c(sq_err(r_een), sq_err(r_lasso), sq_err(r_enet), sq_err(r_alasso))
      ))
    }
    rep_out
  }
}

# ---- Main execution --------------------------------------------------------
cat("EEN empirical application: corn NIR\n")

set.seed(SEED_BASE)
t_cv <- system.time({
  cv_results <- run_cv_comparison(X_raw, y_raw, Omega, n_repeats = N_REPEATS, k_folds = K_FOLDS,
                                  use_parallel = TRUE)
})
cat(sprintf("Cross-validation completed in %.1f seconds.\n", t_cv[3]))

out_dir <- "corn_pcv_moisture_results"
if (!dir.exists(out_dir)) dir.create(out_dir)
write.csv(cv_results, file.path(out_dir, "cv_results_raw.csv"), row.names = FALSE)

per_repeat <- cv_results %>% group_by(Repetition, Method) %>%
  summarise(RMSE = sqrt(mean(MSE)), .groups = "drop")
summary_table <- per_repeat %>% group_by(Method) %>%
  summarise(RMSE_mean = round(mean(RMSE), 4), RMSE_sd = round(sd(RMSE), 4), .groups = "drop") %>%
  arrange(RMSE_mean)
cat("\n--- Cross-validated RMSE ---\n")
print(summary_table)
write.csv(summary_table, file.path(out_dir, "cv_summary.csv"), row.names = FALSE)

test_significance <- function(per_repeat_df, baseline = "EEN") {
  base_data <- per_repeat_df[per_repeat_df$Method == baseline, c("Repetition", "RMSE")]
  out <- data.frame()
  for (m in setdiff(unique(per_repeat_df$Method), baseline)) {
    m_data <- per_repeat_df[per_repeat_df$Method == m, c("Repetition", "RMSE")]
    paired <- merge(base_data, m_data, by = "Repetition", suffixes = c(".base", ".m"))
    if (nrow(paired) < 2) next
    t_res <- t.test(paired$RMSE.m, paired$RMSE.base, paired = TRUE)
    out <- rbind(out, data.frame(
      Comparison = paste(m, "vs", baseline), Mean_diff = round(t_res$estimate, 4),
      CI_lower = round(t_res$conf.int[1], 4), CI_upper = round(t_res$conf.int[2], 4),
      P_value = format.pval(t_res$p.value, digits = 3),
      Significant = ifelse(t_res$p.value < 0.05, "Yes", "No"), n_pairs = nrow(paired)
    ))
  }
  out
}
sig_table <- test_significance(per_repeat)
cat("\n--- Paired t-tests ---\n")
print(sig_table)
write.csv(sig_table, file.path(out_dir, "cv_significance.csv"), row.names = FALSE)

# Full-sample fit
full_center <- colMeans(X_raw); full_scale <- apply(X_raw, 2, sd); full_scale[full_scale == 0] <- 1
X_full_s <- scale(X_raw, center = full_center, scale = full_scale)
y_full_c <- y_raw - mean(y_raw)

set.seed(SEED_BASE)
full_een <- fit_predict_een(X_full_s, y_full_c, X_full_s, Omega, omega_nbrs = omega_nbrs)
full_lasso <- fit_predict_lasso(X_full_s, y_full_c, X_full_s)
full_enet <- fit_predict_enet(X_full_s, y_full_c, X_full_s)
full_alasso <- fit_predict_alasso(X_full_s, y_full_c, X_full_s)

cat("\n--- EEN hyperparameters ---\n")
cat(sprintf("  lambda1 = %.4g%s\n", full_een$lambda1,
            if (full_een$lambda1_at_boundary) "  [at grid boundary]" else ""))
cat(sprintf("  lambda2 = %.4g%s\n", full_een$lambda2,
            if (full_een$lambda2_at_boundary) "  [at grid boundary]" else ""))

roughness_table <- data.frame(
  Method = c("EEN", "LASSO", "Elastic Net", "Adaptive LASSO"),
  Roughness = round(c(roughness(full_een$beta, Omega), roughness(full_lasso$beta, Omega),
                      roughness(full_enet$beta, Omega), roughness(full_alasso$beta, Omega)), 4),
  n_nonzero = c(sum(abs(full_een$beta) > 1e-6), sum(abs(full_lasso$beta) > 1e-6),
                sum(abs(full_enet$beta) > 1e-6), sum(abs(full_alasso$beta) > 1e-6))
) %>% arrange(Roughness)
cat("\n--- Roughness table ---\n")
print(roughness_table)
write.csv(roughness_table, file.path(out_dir, "roughness_table.csv"), row.names = FALSE)

# Coefficient profiles
coef_df <- data.frame(
  wavelength = rep(wavelengths, 4),
  beta = c(full_een$beta, full_lasso$beta, full_enet$beta, full_alasso$beta),
  Method = rep(c("EEN", "LASSO", "Elastic Net", "Adaptive LASSO"), each = p)
)
coef_df$Method <- factor(coef_df$Method, levels = c("EEN", "LASSO", "Elastic Net", "Adaptive LASSO"))

p_coef <- ggplot(coef_df, aes(x = wavelength, y = beta, colour = Method)) +
  geom_line(linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  facet_wrap(~Method, ncol = 1, scales = "free_y") +
  theme_bw(base_size = 11) +
  labs(x = "Wavelength (nm)", y = expression(hat(beta)[j])) +
  theme(legend.position = "none") +
  scale_colour_manual(values = METHOD_COLORS)
ggsave(file.path(out_dir, "coefficient_profiles.pdf"), p_coef, width = 8, height = 9, dpi = 300)

p_rmse <- ggplot(per_repeat, aes(x = reorder(Method, RMSE, mean), y = RMSE, fill = Method)) +
  geom_boxplot(alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, colour = "black") +
  theme_minimal(base_size = 11) +
  labs(x = NULL, y = "RMSE") +
  theme(legend.position = "none") +
  scale_fill_manual(values = METHOD_COLORS)
ggsave(file.path(out_dir, "rmse_comparison.pdf"), p_rmse, width = 7, height = 5, dpi = 300)

cat("\nAll outputs saved to:", out_dir, "\n")
