## Loop controls
max_iter <- fit_controls$max_iter
threshold <- fit_controls$gate_threshold
logit <- function(p) log(p) - log(1-p)

get_Vblock_dense <- function(V_k, r, s, m) {
  rows_r <- ((r - 1) * m + 1):(r * m)
  rows_s <- ((s - 1) * m + 1):(s * m)
  V_k[rows_r, rows_s, drop = FALSE]
}

compute_ES_dense <- function(M_k, W_k, XtX_k, V_k, m, q) {
  E_mat <- matrix(0, q, q)
  S_mat <- matrix(0, q, q)
  
  for (r in seq_len(q)) {
    Mr <- M_k[, r]
    Wr <- W_k[, r]
    
    for (s in r:q) {
      Ms <- M_k[, s]
      Ws <- W_k[, s]
      
      V_rs <- get_Vblock_dense(V_k, r, s, m)
      
      Omega_rs <- outer(Wr, Ws)
      if (r == s) {
        diag(Omega_rs) <- diag(Omega_rs) + Wr * (1 - Wr)
      }
      
      MMt <- outer(Mr, Ms)
      
      val_E <- sum(XtX_k * ((MMt + V_rs) * Omega_rs))
      val_S <- sum(Mr * Ms) + sum(diag(V_rs))
      
      E_mat[r, s] <- val_E
      S_mat[r, s] <- val_S
      
      if (s != r) {
        E_mat[s, r] <- val_E
        S_mat[s, r] <- val_S
      }
    }
  }
  
  list(E_mat = E_mat, S_mat = S_mat)
}

compute_delta_one_k_dense <- function(k, Wcur, M_list, V_k,
                                      Tau_list, XtX_list, XtY_list,
                                      m, q) {
  delta_mat <- matrix(0, nrow = m, ncol = q)
  
  XtX_k <- XtX_list[[k]]
  XtY_k <- XtY_list[[k]]
  M_k   <- M_list[[k]]
  Tau_k <- Tau_list[[k]]
  W_k   <- Wcur[[k]]
  
  s_diag <- diag(XtX_k)
  
  for (r in seq_len(q)) {
    Mr <- M_k[, r]
    Wr <- W_k[, r]
    
    Vrr <- get_Vblock_dense(V_k, r, r, m)
    
    Delta1 <- -0.5 * Tau_k[r, r] * s_diag * (Mr^2 + diag(Vrr))
    Delta2 <- Mr * as.numeric(XtY_k %*% Tau_k[, r])
    
    T3 <- outer(Mr, Mr) + Vrr
    diag(T3) <- 0
    
    Delta3 <- -Tau_k[r, r] *
      rowSums(XtX_k * sweep(T3, 2, Wr, `*`))
    
    Delta4 <- numeric(m)
    
    for (s in other_r_list[[r]]) {
      Ms <- M_k[, s]
      Ws <- W_k[, s]
      
      Vrs <- get_Vblock_dense(V_k, r, s, m)
      T4 <- outer(Mr, Ms) + Vrs
      
      Delta4 <- Delta4 -
        Tau_k[r, s] * rowSums(XtX_k * sweep(T4, 2, Ws, `*`))
    }
    
    delta_mat[, r] <- Delta1 + Delta2 + Delta3 + Delta4
  }
  
  delta_mat
}

top_n_idx <- function(score, n) {
  n <- min(length(score), max(0L, as.integer(n)))
  if (n <= 0L) return(integer(0))
  order(score, decreasing = TRUE)[seq_len(n)]
}

precompute_top1_ld <- function(XtX_k, block_id_vec = NULL, ld_min = 0) {
  m <- nrow(XtX_k)
  out <- vector("list", m)
  snp_sd <- sqrt(pmax(diag(XtX_k), 0))
  
  if (is.null(block_id_vec) || length(block_id_vec) != m) {
    block_list <- list(seq_len(m))
  } else {
    block_list <- split(seq_len(m), block_id_vec)
  }
  
  for (idx_block in block_list) {
    idx_block <- as.integer(idx_block)
    if (length(idx_block) <= 1L) next
    
    denom <- outer(snp_sd[idx_block], snp_sd[idx_block], `*`)
    corr_abs <- abs(XtX_k[idx_block, idx_block, drop = FALSE] / denom)
    corr_abs[!is.finite(corr_abs)] <- 0
    diag(corr_abs) <- 0
    
    for (pos in seq_along(idx_block)) {
      jpos <- which.max(corr_abs[pos, ])
      if (length(jpos) == 1L && corr_abs[pos, jpos] > ld_min) {
        out[[idx_block[pos]]] <- idx_block[jpos]
      } else {
        out[[idx_block[pos]]] <- integer(0)
      }
    }
  }
  
  out
}

get_ld_top1_neighbors <- function(seed_idx, ld_top1_k) {
  if (length(seed_idx) == 0L) return(integer(0))
  unique(unlist(ld_top1_k[seed_idx], use.names = FALSE))
}

compute_ES_active <- function(M_C, W_C, XtX_C, V_C, Tau_k, sigmaB2, m, q) {
  c_k <- nrow(M_C)
  E_mat <- matrix(0, q, q)
  S_mat <- matrix(0, q, q)
  
  for (r in seq_len(q)) {
    rows_r <- ((r - 1) * c_k + 1):(r * c_k)
    Mr <- M_C[, r]
    Wr <- W_C[, r]
    
    for (s in r:q) {
      rows_s <- ((s - 1) * c_k + 1):(s * c_k)
      Ms <- M_C[, s]
      Ws <- W_C[, s]
      
      V_rs <- V_C[rows_r, rows_s, drop = FALSE]
      Omega_rs <- outer(Wr, Ws)
      if (r == s) {
        diag(Omega_rs) <- diag(Omega_rs) + Wr * (1 - Wr)
      }
      
      MMt <- outer(Mr, Ms)
      val_E <- sum(XtX_C * ((MMt + V_rs) * Omega_rs))
      val_S <- sum(Mr * Ms) + sum(diag(V_rs))
      
      E_mat[r, s] <- val_E
      S_mat[r, s] <- val_S
      if (s != r) {
        E_mat[s, r] <- val_E
        S_mat[s, r] <- val_S
      }
    }
  }
  
  if (c_k < m) {
    S_mat <- S_mat + (m - c_k) * sigmaB2 * solve(Tau_k)
  }
  
  list(E_mat = E_mat, S_mat = S_mat)
}

compute_delta_one_k_active <- function(k, Wcur, M_list, V_C, active_idx,
                                       Tau_list, XtX_list, XtY_list,
                                       m, q) {
  delta_mat <- matrix(-Inf, nrow = m, ncol = q)
  c_k <- length(active_idx)
  if (c_k == 0L) return(delta_mat)
  
  XtX_C <- XtX_list[[k]][active_idx, active_idx, drop = FALSE]
  XtY_C <- XtY_list[[k]][active_idx, , drop = FALSE]
  M_C <- M_list[[k]][active_idx, , drop = FALSE]
  Tau_k <- Tau_list[[k]]
  W_C <- Wcur[[k]][active_idx, , drop = FALSE]
  s_diag_C <- diag(XtX_C)
  delta_C <- matrix(0, nrow = c_k, ncol = q)
  
  for (r in seq_len(q)) {
    rows_r <- ((r - 1) * c_k + 1):(r * c_k)
    Mr <- M_C[, r]
    Wr <- W_C[, r]
    Vrr <- V_C[rows_r, rows_r, drop = FALSE]
    
    Delta1 <- -0.5 * Tau_k[r, r] * s_diag_C * (Mr^2 + diag(Vrr))
    Delta2 <- Mr * as.numeric(XtY_C %*% Tau_k[, r])
    
    T3 <- outer(Mr, Mr) + Vrr
    diag(T3) <- 0
    Delta3 <- -Tau_k[r, r] *
      rowSums(XtX_C * sweep(T3, 2, Wr, `*`))
    
    Delta4 <- numeric(c_k)
    for (s in other_r_list[[r]]) {
      rows_s <- ((s - 1) * c_k + 1):(s * c_k)
      Ms <- M_C[, s]
      Ws <- W_C[, s]
      Vrs <- V_C[rows_r, rows_s, drop = FALSE]
      T4 <- outer(Mr, Ms) + Vrs
      
      Delta4 <- Delta4 -
        Tau_k[r, s] * rowSums(XtX_C * sweep(T4, 2, Ws, `*`))
    }
    
    delta_C[, r] <- Delta1 + Delta2 + Delta3 + Delta4
  }
  
  delta_mat[active_idx, ] <- delta_C
  delta_mat
}


elbo_prev <- -Inf
tol_elbo <- 1e-7
option_prev <- NULL

## ============================================================
## Observed-data adaptive gate parameters
## ============================================================

gate_rules <- get_gate_rules(K = K, threshold = threshold)
gate_quantile <- gate_rules$gate_quantile
local_ok_rule <- gate_rules$local_ok_rule
c_W_gate <- gate_rules$c_W_gate
option2_alpha <- gate_rules$option2_alpha
option2_meta_alpha <- gate_rules$option2_meta_alpha
option2_current_z_min <- gate_rules$option2_current_z_min
option2_current_alpha <- gate_rules$option2_current_alpha
option2_oppose_alpha <- gate_rules$option2_oppose_alpha
option3_other_alpha <- gate_rules$option3_other_alpha
option3_current_alpha <- gate_rules$option3_current_alpha
min_other_signal <- gate_rules$min_other_signal
option3_margin <- gate_rules$option3_margin

## ============================================================
## Hard-switch thresholds calibrated from simulation diagnostics
## ============================================================

safe_quantile <- function(x, prob, default = NA_real_) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(default)
  as.numeric(quantile(x, probs = prob, na.rm = TRUE, names = FALSE))
}

# Initialize Tau as q × q identity matrices
Tau_list <- vector("list", K)
for (k in 1:K) {
  Tau_list[[k]] <- diag(q)   # 3x3
}
# Initialize M_list
M_list <- vector("list", K)

for (k in 1:K) {
  M_list[[k]] <- matrix(0, m, q)  # 初始化为零矩阵
}

compute_marginal_z <- function(Xk, Yk) {
  Xc <- scale(Xk, center = TRUE, scale = FALSE)
  Yc <- scale(Yk, center = TRUE, scale = FALSE)
  
  Xc <- as.matrix(Xc)
  Yc <- as.matrix(Yc)
  Xc[!is.finite(Xc)] <- 0
  Yc[!is.finite(Yc)] <- 0
  
  n_obs <- nrow(Xc)
  xtx <- colSums(Xc^2)
  xty <- crossprod(Xc, Yc)
  yty <- colSums(Yc^2)
  
  z_mat <- matrix(0, nrow = ncol(Xc), ncol = ncol(Yc))
  
  valid_snp <- is.finite(xtx) & xtx > 1e-12
  
  if (any(valid_snp)) {
    beta_hat <- matrix(0, nrow = ncol(Xc), ncol = ncol(Yc))
    fit_ss   <- matrix(0, nrow = ncol(Xc), ncol = ncol(Yc))
    
    beta_hat[valid_snp, ] <- sweep(
      xty[valid_snp, , drop = FALSE],
      1,
      xtx[valid_snp],
      "/"
    )
    
    fit_ss[valid_snp, ] <- sweep(
      xty[valid_snp, , drop = FALSE]^2,
      1,
      xtx[valid_snp],
      "/"
    )
    
    rss <- matrix(yty, nrow = ncol(Xc), ncol = ncol(Yc), byrow = TRUE) -
      fit_ss
    
    rss[!is.finite(rss)] <- 0
    rss <- pmax(rss, .Machine$double.eps)
    
    sigma2_hat <- pmax(rss / max(n_obs - 2, 1), .Machine$double.eps)
    
    se_hat <- matrix(Inf, nrow = ncol(Xc), ncol = ncol(Yc))
    se_hat[valid_snp, ] <- sqrt(
      sweep(
        sigma2_hat[valid_snp, , drop = FALSE],
        1,
        xtx[valid_snp],
        "/"
      )
    )
    
    z_mat <- beta_hat / se_hat
  }
  
  z_mat[!is.finite(z_mat)] <- 0
  z_mat
}
marginal_z_list <- vector("list", K)

for (k in seq_len(K)) {
  marginal_z_list[[k]] <- compute_marginal_z(X_list[[k]], Y_list[[k]])
}

loo_k_list <- lapply(seq_len(K), function(k) setdiff(seq_len(K), k))
other_r_list <- lapply(seq_len(q), function(r) setdiff(seq_len(q), r))
sqrt_n_k <- sqrt(n_k)

marginal_z_array <- array(NA_real_, dim = c(K, m, q))
for (k in seq_len(K)) {
  marginal_z_array[k, , ] <- marginal_z_list[[k]]
}
p_marg_fixed_array <- 2 * pnorm(-abs(marginal_z_array))
z_abs_fixed_array <- abs(marginal_z_array)

p_meta_all_fixed_mat <- matrix(NA_real_, nrow = m, ncol = q)
for (i in seq_len(m)) {
  for (r in seq_len(q)) {
    z_vec_all_ir <- marginal_z_array[, i, r]
    meta_z_all <- sum(sqrt_n_k * z_vec_all_ir) / sqrt(sum(n_k))
    p_meta_all_fixed_mat[i, r] <- 2 * pnorm(-abs(meta_z_all))
  }
}

active_controls_default <- list(
  use_active = TRUE,
  open_iter = 10L,
  drop_eps = 5e-4,
  c_min = 50L,
  c_frac_n = 1 ,
  z_frac = 0.20,
  
  ## LD neighbor augmentation
  use_ld_neighbor = TRUE,
  ld_min = 0.7,
  
  ## Diagnostics
  store_active_history = FALSE
)
active_controls <- active_controls_default
if (!is.null(fit_controls$active_controls)) {
  active_controls[names(fit_controls$active_controls)] <- fit_controls$active_controls
}

active_use <- isTRUE(active_controls$use_active)
active_open_iter <- as.integer(active_controls$open_iter)
active_drop_eps <- as.numeric(active_controls$drop_eps)
active_c_min <- as.integer(active_controls$c_min)
active_c_frac_n <- as.numeric(active_controls$c_frac_n)
active_z_frac <- as.numeric(active_controls$z_frac)
active_use_ld_neighbor <- isTRUE(active_controls$use_ld_neighbor)
active_ld_min <- as.numeric(active_controls$ld_min)

c_target_k <- pmin(
  m,
  pmax(active_c_min, floor(n_k * active_c_frac_n))
)
c_target_k <- as.integer(pmax(1L, c_target_k))
c_z_k <- as.integer(pmin(c_target_k, pmax(5L, ceiling(active_z_frac * c_target_k))))
c_w_k <- as.integer(pmax(1L, c_target_k - c_z_k))
active_min_size <- as.integer(min(m, max(1L, ceiling(0.05 * m))))

z_score_snp_list <- vector("list", K)
Z_watch_list <- vector("list", K)
C_list <- vector("list", K)
drop_count_list <- vector("list", K)
V_active_list <- vector("list", K)
active_idx_list <- vector("list", K)
active_frozen <- rep(FALSE, K)

ld_top1_list <- vector("list", K)

if (active_use_ld_neighbor) {
  for (k in seq_len(K)) {
    ld_top1_list[[k]] <- precompute_top1_ld(
      XtX_k = XtX_list[[k]],
      block_id_vec = if (exists("block_id_vec")) block_id_vec else NULL,
      ld_min = active_ld_min
    )
  }
} else {
  for (k in seq_len(K)) {
    ld_top1_list[[k]] <- vector("list", m)
  }
}

for (k in seq_len(K)) {
  z_score_snp_list[[k]] <- apply(z_abs_fixed_array[k, , , drop = FALSE], 2, max)
  Z_watch_list[[k]] <- top_n_idx(z_score_snp_list[[k]], c_z_k[k])
  
  W_score <- apply(W_list[[k]], 1, max)
  W_top <- top_n_idx(W_score, c_w_k[k])
  
  seed_idx <- sort(unique(c(W_top, Z_watch_list[[k]])))
  LD_top <- get_ld_top1_neighbors(seed_idx, ld_top1_list[[k]])
  
  C_init <- sort(unique(c(seed_idx, LD_top)))
  
  if (length(C_init) < c_target_k[k]) {
    filler_score <- W_score + 1e-6 * z_score_snp_list[[k]]
    filler <- setdiff(top_n_idx(filler_score, m), C_init)
    C_init <- sort(unique(c(C_init, head(filler, c_target_k[k] - length(C_init)))))
  }
  
  C_list[[k]] <- C_init
  active_idx_list[[k]] <- C_init
  drop_count_list[[k]] <- integer(m)
  
  inactive_idx <- setdiff(seq_len(m), C_init)
  if (length(inactive_idx) > 0L) {
    W_list[[k]][inactive_idx, ] <- 0
  }
}

## Current-iteration arrays used for option-specific W updates and summaries.
lambda_opt2_all_array <- array(NA_real_, dim = c(K, m, q))
lambda_opt3_all_array <- array(NA_real_, dim = c(K, m, q))
lambda_update_array <- array(NA_real_, dim = c(K, m, q))
W_after_update_array <- array(NA_real_, dim = c(K, m, q))

iter_summary <- data.frame(
  iter = 1:max_iter,
  ELBO = NA_real_,
  
  option1_count = NA_integer_,
  option2_count = NA_integer_,
  option3_count = NA_integer_,
  option_change = NA_integer_,
  
  tau_data = NA_real_,
  tau_z = NA_real_,
  tau_W2 = NA_real_,
  tau_W3 = NA_real_,
  tau_cross_ratio = NA_real_,
  n_local_ok = NA_integer_,
  n_option2_ok = NA_integer_,
  n_option3_ok = NA_integer_,
  
  G_active = NA_integer_,
  G_change = NA_integer_,
  mean_delta_gate_opt2 = NA_real_,
  mean_p_pool_opt2 = NA_real_,
  mean_lambda_opt2 = NA_real_,
  mean_W_after_opt2 = NA_real_,
  mean_benefit_prob_G1 = NA_real_,
  mean_benefit_prob_G0 = NA_real_,
  n_decision_improve_G1 = NA_integer_,
  n_decision_harm_G1 = NA_integer_,
  n_decision_improve_G0 = NA_integer_,
  n_decision_harm_G0 = NA_integer_,
  mean_score_gain_G1 = NA_real_,
  mean_score_gain_G0 = NA_real_,
  G_true_active_rate = NA_real_,
  stringsAsFactors = FALSE
)
## ============================================================
## Store active-set history for diagnostics
## ============================================================

active_size_history <- matrix(
  NA_integer_,
  nrow = max_iter,
  ncol = K
)
colnames(active_size_history) <- paste0("pop", seq_len(K))

store_active_history <- isTRUE(active_controls$store_active_history)

active_idx_history <- if (store_active_history) {
  vector("list", max_iter)
} else {
  NULL
}
YtY_list <- vector("list", K)
for (k in seq_len(K)) {
  YtY_list[[k]] <- crossprod(Y_list[[k]])
}

for (it in 1:max_iter) {
  logdet_V_list <- vector("list", K)
  logdet_Psi_list <- vector("list", K)
  delta_gate_array <- array(-Inf, dim = c(K, m, q))
  
  if (it > 1L && it <= active_open_iter) {
    for (k in seq_len(K)) {
      if (active_frozen[k]) next
      
      current_active <- C_list[[k]]
      W_score <- apply(W_list[[k]], 1, max)
      W_top <- top_n_idx(W_score, c_w_k[k])
      
      seed_idx <- sort(unique(c(W_top, Z_watch_list[[k]])))
      LD_top <- get_ld_top1_neighbors(seed_idx, ld_top1_list[[k]])
      
      C_new <- sort(unique(c(seed_idx, LD_top)))
      
      low_now <- rep(FALSE, m)
      if (length(current_active) > 0L) {
        low_now[current_active] <- W_score[current_active] < active_drop_eps
      }
      drop_count_list[[k]][low_now] <- drop_count_list[[k]][low_now] + 1L
      drop_count_list[[k]][!low_now] <- 0L
      
      drop_idx <- which(drop_count_list[[k]] >= 2L)
      C_new <- setdiff(C_new, drop_idx)
      
      if (length(C_new) < active_min_size) {
        C_list[[k]] <- sort(top_n_idx(z_score_snp_list[[k]], active_min_size))
        drop_count_list[[k]] <- integer(m)
        active_frozen[k] <- TRUE
        next
      }
      
      if (length(C_new) < c_target_k[k]) {
        filler_score <- W_score + 1e-6 * z_score_snp_list[[k]]
        filler <- setdiff(top_n_idx(filler_score, m), c(C_new, drop_idx))
        C_new <- sort(unique(c(C_new, head(filler, c_target_k[k] - length(C_new)))))
      }
      
      C_list[[k]] <- sort(C_new)
    }
  }
  
  ## Record current computational active sets
  for (kk in seq_len(K)) {
    active_size_history[it, kk] <- length(C_list[[kk]])
  }
  
  if (store_active_history) {
    active_idx_history[[it]] <- lapply(C_list, identity)
  }
  
  for (k in 1:K) {
    idx <- C_list[[k]]
    c_k <- length(idx)
    active_idx_list[[k]] <- idx
    
    W_full <- matrix(0, m, q)
    W_full[idx, ] <- W_list[[k]][idx, , drop = FALSE]
    W_list[[k]] <- W_full
    
    W_C <- W_full[idx, , drop = FALSE]
    wk <- as.vector(W_C)
    Tau_k <- Tau_list[[k]]
    XtX_C <- XtX_list[[k]][idx, idx, drop = FALSE]
    XtY_C <- XtY_list[[k]][idx, , drop = FALSE]
    X_k <- X_list[[k]]
    Y_k <- Y_list[[k]]
    
    K_prior_k <- kronecker(Tau_k, diag(c_k) / sigmaB2)
    K_like_k  <- kronecker(Tau_k, XtX_C)
    
    Vinv_k <- K_prior_k + sweep(sweep(K_like_k, 1, wk, `*`), 2, wk, `*`)
    
    diag_like_k <- rep(diag(Tau_k), each = c_k) * rep(diag(XtX_C), times = q)
    diag(Vinv_k) <- diag(Vinv_k) + diag_like_k * wk * (1 - wk)
    
    chol_Vinv_k <- chol(Vinv_k)
    
    logdet_active <- -2 * sum(log(diag(chol_Vinv_k)))
    logdet_tau <- 2 * sum(log(diag(chol(Tau_k))))
    logdet_inactive <- (m - c_k) * (q * log(sigmaB2) - logdet_tau)
    logdet_V_list[[k]] <- logdet_active + logdet_inactive
    
    V_C <- chol2inv(chol_Vinv_k)
    V_active_list[[k]] <- V_C
    
    RHS_mat_C <- XtY_C %*% Tau_k
    rhs_C <- as.vector(W_C * RHS_mat_C)
    
    tmp_k <- forwardsolve(t(chol_Vinv_k), rhs_C)
    m_new_C <- as.numeric(backsolve(chol_Vinv_k, tmp_k))
    
    M_C <- matrix(m_new_C, c_k, q)
    M_full <- matrix(0, m, q)
    M_full[idx, ] <- M_C
    M_list[[k]] <- M_full
    
    ES_k <- compute_ES_active(
      M_C = M_C,
      W_C = W_C,
      XtX_C = XtX_C,
      V_C = V_C,
      Tau_k = Tau_k,
      sigmaB2 = sigmaB2,
      m = m,
      q = q
    )
    
    E_mat_k <- ES_k$E_mat
    S_mat_k <- ES_k$S_mat
    
    ## ===== Update Tau_k = E[Sigma_k^{-1}] =====
    Bhat_C <- W_C * M_C
    
    XBhat_C <- X_k[, idx, drop = FALSE] %*% Bhat_C
    
    YtX_Bhat_k <- crossprod(Y_k, XBhat_C)
    
    Bt_XtY_k <- t(Bhat_C) %*% XtY_C
    
    YtY_k <- YtY_list[[k]]
    
    Psi_scale_k <- Psi0 + YtY_k -
      Bt_XtY_k -
      YtX_Bhat_k +
      E_mat_k +
      (1 / sigmaB2) * S_mat_k
    
    Psi_scale_k[!is.finite(Psi_scale_k)] <- 0
    Psi_scale_k <- Psi_scale_k + 1e-6 * diag(nrow(Psi_scale_k))
    
    chol_Psi <- chol(Psi_scale_k)
    logdet_Psi_list[[k]] <- 2 * sum(log(diag(chol_Psi)))
    
    nu_star <- nu0 + n_k[k] + m
    Tau_list[[k]] <- nu_star * chol2inv(chol_Psi)
    
    delta_gate_array[k, , ] <- compute_delta_one_k_active(
      k = k,
      Wcur = W_list,
      M_list = M_list,
      V_C = V_C,
      active_idx = idx,
      Tau_list = Tau_list,
      XtX_list = XtX_list,
      XtY_list = XtY_list,
      m = m,
      q = q
    )
  }
  
  W_temp <- vector("list", K)
  for (k in 1:K) {
    W_temp[[k]] <- matrix(0, m, q)
  }
  
  ## =====  (4) Update W (element-wise)  =====
  a <- 1
  b <- 1
  
  lambda_opt1 <- logit(pi_fix)
  
  epsS <- 1e-6
  
  ## ============================================================
  ## Build observed-data adaptive option before W update
  ## option_array[k,i,r] = 1, 2, or 3
  ## No truth is used to decide option_array.
  ## ============================================================
  
  option_array <- array(1L, dim = c(K, m, q))
  
  W_array <- array(NA_real_, dim = c(K, m, q))
  for (kk in seq_len(K)) {
    W_array[kk, , ] <- W_list[[kk]]
  }
  active_entry_mask <- array(FALSE, dim = c(K, m, q))
  for (kk in seq_len(K)) {
    active_entry_mask[kk, active_idx_list[[kk]], ] <- TRUE
  }
  
  p_data_gate_array <- plogis(delta_gate_array)
  p_pool_gate_array <- array(NA_real_, dim = c(K, m, q))
  z_abs_gate_array <- z_abs_fixed_array
  W_current_array <- W_array
  
  W_mean_loo_array <- array(NA_real_, dim = c(K, m, q))
  meta_z_loo_array <- array(NA_real_, dim = c(K, m, q))
  p_marg_array <- p_marg_fixed_array
  p_meta_loo_array <- array(NA_real_, dim = c(K, m, q))
  
  ## meta-z and meta-p using all populations, used for Option 3
  p_meta_all_mat <- p_meta_all_fixed_mat
  
  crosspheno_W_mean_array <- array(NA_real_, dim = c(K, m, q))
  
  local_ok_array <- array(FALSE, dim = c(K, m, q))
  option2_ok_array <- array(FALSE, dim = c(K, m, q))
  option3_ok_array <- array(FALSE, dim = c(K, m, q))
  
  ## ============================================================
  ## First pass: compute observed features only
  ## ============================================================
  
  for (k in 1:K) {
    for (i in active_idx_list[[k]]) {
      for (r in 1:q) {
        
        delta_tmp <- delta_gate_array[k, i, r]
        p_data <- p_data_gate_array[k, i, r]
        
        z_vec_all <- marginal_z_array[, i, r]
        W_vec_all <- W_array[, i, r]
        
        z_kir <- z_vec_all[k]
        z_abs <- z_abs_gate_array[k, i, r]
        
        p_kir <- p_marg_array[k, i, r]
        p_pool_gate_array[k, i, r] <- mean(W_vec_all)
        
        loo_idx <- loo_k_list[[k]]
        z_vec_loo <- z_vec_all[loo_idx]
        W_vec_loo <- W_vec_all[loo_idx]
        n_vec_loo <- n_k[loo_idx]
        
        W_mean_loo <- mean(W_vec_loo)
        
        W_mean_loo_array[k, i, r] <- W_mean_loo
        
        meta_z_loo <- sum(sqrt(n_vec_loo) * z_vec_loo) / sqrt(sum(n_vec_loo))
        
        p_meta_loo <- 2 * pnorm(-abs(meta_z_loo))
        p_meta_loo_array[k, i, r] <- p_meta_loo
        
        meta_z_loo_array[k, i, r] <- meta_z_loo
        
        crosspheno_values <- as.numeric(
          W_array[, i, other_r_list[[r]], drop = FALSE]
        )
        
        crosspheno_W_mean <- mean(crosspheno_values)
        
        crosspheno_W_mean_array[k, i, r] <- crosspheno_W_mean
        
      }
    }
  }
  
  ## ============================================================
  ## Adaptive thresholds from observed-data empirical distributions
  ## ============================================================
  
  ## Background-like set, using only observed quantities.
  ## This is not truth. It is used to stabilize thresholds.
  bg_mask <- (W_current_array <= 0.01) & active_entry_mask
  
  if (sum(bg_mask, na.rm = TRUE) < 100) {
    bg_mask <- active_entry_mask
  }
  
  tau_data_adapt <- max(
    safe_quantile(p_data_gate_array[bg_mask], gate_quantile, default = 0),
    safe_quantile(p_data_gate_array[active_entry_mask], 0.75, default = 0)
  )
  
  tau_z_adapt <- max(
    safe_quantile(z_abs_gate_array[bg_mask], gate_quantile, default = 0),
    safe_quantile(z_abs_gate_array[active_entry_mask], 0.75, default = 0)
  )
  
  tau_W2_adapt <- max(
    safe_quantile(W_mean_loo_array[bg_mask], gate_quantile, default = 0),
    c_W_gate
  )
  
  tau_W3_adapt <- max(
    safe_quantile(crosspheno_W_mean_array[bg_mask], gate_quantile, default = 0),
    c_W_gate
  )
  
  tau_cross_ratio <- 1 + option3_margin
  
  ## ============================================================
  ## Second pass: observed-data option decision
  ## ============================================================
  
  for (k in 1:K) {
    for (i in active_idx_list[[k]]) {
      for (r in 1:q) {
        
        p_data <- p_data_gate_array[k, i, r]
        z_abs <- z_abs_gate_array[k, i, r]
        
        if (local_ok_rule == "both") {
          local_ok <- (p_data > tau_data_adapt) && (z_abs > tau_z_adapt)
        } else {
          local_ok <- (p_data > tau_data_adapt) || (z_abs > tau_z_adapt)
        }
        
        ## ============================================================
        ## Simplified Option 2 gate:
        ## same SNP-phenotype pair replicated across populations
        ## using marginal p-values.
        ## ============================================================
        
        loo_idx <- loo_k_list[[k]]
        z_vec_loo <- marginal_z_array[loo_idx, i, r]
        p_vec_loo <- p_marg_array[loo_idx, i, r]
        
        sig_idx <- which(p_vec_loo < option2_alpha)
        N_p_loo <- length(sig_idx)
        
        ## Direction check among external populations that pass p threshold.
        ## If N_p_loo >= 2, all significant external z's should have same sign.
        if (N_p_loo >= 2) {
          sig_signs <- sign(z_vec_loo[sig_idx])
          sig_signs <- sig_signs[sig_signs != 0]
          opt2_sign_ok <- length(unique(sig_signs)) <= 1
        } else if (N_p_loo == 1) {
          ## If only one external population is significant,
          ## its sign should agree with the LOO meta-z sign.
          opt2_sign_ok <- sign(z_vec_loo[sig_idx]) == sign(meta_z_loo_array[k, i, r])
        } else {
          opt2_sign_ok <- FALSE
        }
        
        ## Strong replicated evidence:
        ## at least two external populations have p < option2_alpha.
        opt2_two_pop_ok <- (N_p_loo >= 2) && opt2_sign_ok
        
        ## Weaker replicated evidence:
        ## only one external population passes p threshold,
        ## but the LOO meta p-value is strong.
        opt2_one_pop_meta_ok <- (
          N_p_loo == 1 &&
            opt2_sign_ok &&
            p_meta_loo_array[k, i, r] < option2_meta_alpha
        )
        
        ## Optional safety guard:
        ## if current population has nominally significant opposite direction,
        ## do not borrow.
        current_z <- marginal_z_array[k, i, r]
        current_p <- p_marg_array[k, i, r]
        meta_z_loo_now <- meta_z_loo_array[k, i, r]
        
        current_support_weak <- (
          abs(current_z) >= option2_current_z_min ||
            current_p < option2_current_alpha
        )
        
        current_opposite_strong <- (
          current_p < option2_oppose_alpha &&
            sign(current_z) != 0 &&
            sign(meta_z_loo_now) != 0 &&
            sign(current_z) != sign(meta_z_loo_now)
        )
        
        samepheno_ok <- (
          (opt2_two_pop_ok || opt2_one_pop_meta_ok) &&
            current_support_weak &&
            !current_opposite_strong
        )
        
        ## ============================================================
        ## Simplified Option 3 gate:
        ## same SNP shows evidence across multiple phenotypes.
        ## ============================================================
        
        ## Evidence in other phenotypes for the same SNP
        other_pheno <- other_r_list[[r]]
        
        D_other_p <- sum(p_meta_all_mat[i, other_pheno] < option3_other_alpha)
        
        ## Current phenotype should not be completely unsupported.
        ## It can be supported either in all-population meta analysis
        ## or in the current population marginal test.
        current_pheno_has_evidence <- (
          p_meta_all_mat[i, r] < option3_current_alpha ||
            p_marg_array[k, i, r] < option3_current_alpha
        )
        
        crosspheno_ok <- (
          D_other_p >= min_other_signal &&
            current_pheno_has_evidence
        )
        
        ## ============================================================
        ## Final simplified hard switch
        ## Option 2 is prioritized because current FN mainly comes from
        ## under-using Option 2; Option 3 is secondary and more conservative.
        ## ============================================================
        
        observed_option <- 1L
        
        if (samepheno_ok) {
          observed_option <- 2L
        } else if (crosspheno_ok) {
          observed_option <- 3L
        } else {
          observed_option <- 1L
        }
        
        option_array[k, i, r] <- observed_option
        
        local_ok_array[k, i, r] <- local_ok
        option2_ok_array[k, i, r] <- samepheno_ok
        option3_ok_array[k, i, r] <- crosspheno_ok
        
        ## Option-specific lambdas used by the W update.
        delta_tmp <- delta_gate_array[k, i, r]
        
        ## Option 2 lambda: all-population same phenotype pooling
        S2_all_cf <- sum(W_array[, i, r])
        S2_all_cf <- min(max(S2_all_cf, epsS), K - epsS)
        
        lambda2_all_cf <- digamma(a + S2_all_cf) -
          digamma(b + K - S2_all_cf)
        
        ## Option 3 lambda: same SNP across populations and phenotypes
        W_snp_values <- as.numeric(W_array[, i, , drop = FALSE])
        
        S3_all_cf <- sum(W_snp_values)
        L3_all_cf <- K * q
        S3_all_cf <- min(max(S3_all_cf, epsS), L3_all_cf - epsS)
        
        lambda3_all_cf <- digamma(a + S3_all_cf) -
          digamma(b + L3_all_cf - S3_all_cf)
        
        lambda_opt2_all_array[k, i, r] <- lambda2_all_cf
        lambda_opt3_all_array[k, i, r] <- lambda3_all_cf
      }
    }
  }
  
  ## Store adaptive thresholds in iter_summary
  iter_summary$tau_data[it] <- tau_data_adapt
  iter_summary$tau_z[it] <- tau_z_adapt
  iter_summary$tau_W2[it] <- tau_W2_adapt
  iter_summary$tau_W3[it] <- tau_W3_adapt
  iter_summary$tau_cross_ratio[it] <- tau_cross_ratio
  
  iter_summary$n_local_ok[it] <- sum(local_ok_array, na.rm = TRUE)
  iter_summary$n_option2_ok[it] <- sum(option2_ok_array, na.rm = TRUE)
  iter_summary$n_option3_ok[it] <- sum(option3_ok_array, na.rm = TRUE)
  
  for (k in 1:K) {                     # population循环
    for (i in active_idx_list[[k]]) {
      for (r in 1:q) {
        
        # ==== Step 2. reuse δ_{k,ir}^{(t)} from the observed-feature pass ====
        delta_ir <- delta_gate_array[k, i, r]
        
        ## ==== Step 3. observed-data adaptive hard switch among Option 1, 2, and 3 ====
        ## option_array is determined only by observed features, not by B_true.
        
        option_kir <- option_array[k, i, r]
        
        lambda_base <- lambda_opt1
        
        if (option_kir == 1) {
          lambda_ir <- lambda_base
          
        } else if (option_kir == 2) {
          lambda_ir <- lambda_opt2_all_array[k, i, r]
          
        } else if (option_kir == 3) {
          lambda_ir <- lambda_opt3_all_array[k, i, r]
          
        } else {
          stop("Unknown option_kir")
        }
        
        eta_ir <- lambda_ir + delta_ir
        W_new_kr <- plogis(eta_ir)
        
        lambda_update_array[k, i, r] <- lambda_ir
        W_after_update_array[k, i, r] <- W_new_kr
        
        W_temp[[k]][i, r] <- W_new_kr
      }# end for r
    }# end for i
  }# end for k
  # after the loop, set W_list <- W_temp  (the new W^{(t+1)})
  for (k in 1:K) {
    W_list[[k]] <- W_temp[[k]]
  }
  
  
  
  ## ===== 计算 ELBO (完整版) =====
  
  elbo_val <- 0
  
  # 第一部分: ∑_k [0.5 log|V_k| - 0.5 nu_star log|Psi_scale_k|]
  for (k in 1:K) {
    # 计算 V_k 的对数行列式
    elbo_val <- elbo_val + 0.5 * logdet_V_list[[k]]
    
    # 计算 Psi_scale_k 的对数行列式
    nu_star_k <- nu0 + n_k[k] + m
    elbo_val <- elbo_val - 0.5 * nu_star_k * logdet_Psi_list[[k]]
  }
  
  ## (A) Bernoulli entropy for all W_{kir}
  W_all <- array(0, dim = c(m, q, K))
  for (k in 1:K) {
    W_all[, , k] <- W_list[[k]]
  }
  
  mask <- (W_all > 1e-12) & (W_all < 1 - 1e-12)
  
  elbo_val <- elbo_val -
    sum(W_all[mask] * log(W_all[mask]) +
          (1 - W_all[mask]) * log(1 - W_all[mask]))
  
  ## (B) Option 2 collapsed prior:
  ## same (i,r) across populations, only for entries assigned to Option 2
  for (i in 1:m) {
    for (r in 1:q) {
      opt2_vec <- option_array[, i, r] == 2
      L2 <- sum(opt2_vec)
      
      if (L2 > 0) {
        S2 <- sum(opt2_vec * W_all[i, r, ])
        
        elbo_val <- elbo_val +
          lgamma(a + S2) +
          lgamma(b + L2 - S2) -
          lgamma(a + b + L2) -
          lgamma(a) - lgamma(b) + lgamma(a + b)
      }
    }
  }
  
  ## (C) Option 3 collapsed prior:
  ## same SNP i across populations and phenotypes,
  ## only for entries assigned to Option 3
  for (i in 1:m) {
    S3 <- 0
    L3 <- 0
    
    for (k in 1:K) {
      for (r in 1:q) {
        if (option_array[k, i, r] == 3) {
          S3 <- S3 + W_list[[k]][i, r]
          L3 <- L3 + 1
        }
      }
    }
    
    if (L3 > 0) {
      elbo_val <- elbo_val +
        lgamma(a + S3) +
        lgamma(b + L3 - S3) -
        lgamma(a + b + L3) -
        lgamma(a) - lgamma(b) + lgamma(a + b)
    }
  }
  
  ## (D) Option 1 Bernoulli prior:
  ## fixed sparse prior for entries assigned to Option 1
  log_pi0 <- log(pi_fix)
  log_1m_pi0 <- log1p(-pi_fix)
  
  for (i in 1:m) {
    for (r in 1:q) {
      for (k in 1:K) {
        if (option_array[k, i, r] == 1) {
          wkir <- W_list[[k]][i, r]
          
          elbo_val <- elbo_val +
            wkir * log_pi0 +
            (1 - wkir) * log_1m_pi0
        }
      }
    }
  }
  
  elbo_curr <- elbo_val
  option1_count_now <- sum(option_array == 1)
  option2_count_now <- sum(option_array == 2)
  option3_count_now <- sum(option_array == 3)
  
  if (it > 1) {
    option_change_now <- sum(option_array != option_prev)
  } else {
    option_change_now <- NA_integer_
  }
  
  ## Backward-compatible names
  G_active_now <- option2_count_now + option3_count_now
  G_change_now <- option_change_now
  
  opt2_mask <- (option_array == 2)
  
  mean_delta_gate_opt2 <- if (sum(opt2_mask) > 0) {
    mean(delta_gate_array[opt2_mask])
  } else NA_real_
  
  mean_p_pool_opt2 <- if (sum(opt2_mask) > 0) {
    mean(p_pool_gate_array[opt2_mask])
  } else NA_real_
  
  mean_lambda_opt2 <- if (sum(opt2_mask) > 0) {
    mean(lambda_update_array[opt2_mask], na.rm = TRUE)
  } else NA_real_
  
  mean_W_after_opt2 <- if (sum(opt2_mask) > 0) {
    mean(W_after_update_array[opt2_mask], na.rm = TRUE)
  } else NA_real_
  
  iter_summary$ELBO[it] <- elbo_curr
  
  iter_summary$option1_count[it] <- option1_count_now
  iter_summary$option2_count[it] <- option2_count_now
  iter_summary$option3_count[it] <- option3_count_now
  iter_summary$option_change[it] <- option_change_now
  
  iter_summary$G_active[it] <- G_active_now
  iter_summary$G_change[it] <- G_change_now
  iter_summary$mean_delta_gate_opt2[it] <- mean_delta_gate_opt2
  iter_summary$mean_p_pool_opt2[it] <- mean_p_pool_opt2
  iter_summary$mean_lambda_opt2[it] <- mean_lambda_opt2
  iter_summary$mean_W_after_opt2[it] <- mean_W_after_opt2
  
  cat(sprintf(
    paste0(
      "Iter %d: ELBO=%.7f | opt1=%d | opt2=%d | opt3=%d | option_change=%s | ",
      "mean_delta_opt2=%.3f | mean_p_pool_opt2=%.3f | ",
      "mean_lambda_opt2=%.3f | mean_W_after_opt2=%.3f\n"
    ),
    it,
    elbo_curr,
    option1_count_now,
    option2_count_now,
    option3_count_now,
    ifelse(is.na(option_change_now), "NA", as.character(option_change_now)),
    mean_delta_gate_opt2,
    mean_p_pool_opt2,
    mean_lambda_opt2,
    mean_W_after_opt2
  ))
  # 收敛判断（统一处理）
  if (it > 1) {
    elbo_diff <- abs(elbo_curr - elbo_prev)
    if (elbo_diff < tol_elbo) {
      cat(sprintf("Converged at iter %d | ΔELBO = %.3e\n", it, elbo_diff))
      break
    }
  }
  elbo_prev <- elbo_curr
  option_prev <- option_array
  
}  # 结束 for (it in 1:max_iter) 循环
n_iter_done <- it


