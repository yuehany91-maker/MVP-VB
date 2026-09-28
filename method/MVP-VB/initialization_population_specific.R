precompute_XtX_XtY <- function(X_list, Y_list) {
  K <- length(X_list)
  XtX_list <- vector("list", K)
  XtY_list <- vector("list", K)
  
  for (k in seq_len(K)) {
    Xk <- X_list[[k]]
    Yk <- Y_list[[k]]
    XtX_list[[k]] <- crossprod(Xk)
    XtY_list[[k]] <- crossprod(Xk, Yk)
  }
  
  list(XtX_list = XtX_list, XtY_list = XtY_list)
}

initialize_stepwise_halo <- function(X_list, Y_list, q, K,
                                     p_enter = 1e-4,
                                     max_steps = 10,
                                     halo_weight = 0.1) {
  X_pool <- do.call(rbind, X_list)
  Y_pool <- do.call(rbind, Y_list)
  X_pool <- scale(X_pool)
  Y_pool <- scale(Y_pool)
  Cor_pool <- cor(X_pool)
  Cor_pool[!is.finite(Cor_pool)] <- 0
  diag(Cor_pool) <- 1
  
  m <- ncol(X_pool)
  W_base <- matrix(0.001, m, q)
  W_init_type <- matrix("background", nrow = m, ncol = q)
  all_selected_snps <- integer(0)
  
  for (r in seq_len(q)) {
    y_resid <- Y_pool[, r]
    current_features <- integer(0)
    
    for (step in seq_len(max_steps)) {
      scores <- abs(crossprod(X_pool, y_resid))
      if (length(current_features) > 0) scores[current_features] <- -1
      best_idx <- which.max(scores)
      
      fit <- lm(y_resid ~ X_pool[, best_idx])
      pval <- summary(fit)$coefficients[2, 4]
      
      if (is.finite(pval) && pval < p_enter) {
        current_features <- c(current_features, best_idx)
        y_resid <- residuals(fit)
      } else {
        break
      }
    }
    
    if (length(current_features) > 0) {
      W_base[current_features, r] <- 0.85
      W_init_type[current_features, r] <- "direct"
      all_selected_snps <- c(all_selected_snps, current_features)
      
      for (feat in current_features) {
        ld_neighbors <- which(abs(Cor_pool[feat, ]) > 0.6)
        ld_neighbors <- setdiff(ld_neighbors, current_features)
        mask <- W_base[ld_neighbors, r] < halo_weight
        halo_idx <- ld_neighbors[mask]
        
        if (length(halo_idx) > 0) {
          W_base[halo_idx, r] <- halo_weight
          W_init_type[halo_idx, r] <- "ld_halo"
        }
      }
    }
  }
  
  all_selected_snps <- unique(all_selected_snps)
  if (length(all_selected_snps) > 0) {
    for (idx in all_selected_snps) {
      for (r in seq_len(q)) {
        if (W_base[idx, r] < 0.8) {
          old_w <- W_base[idx, r]
          W_base[idx, r] <- max(W_base[idx, r], 0.05)
          if (old_w < 0.05 && W_init_type[idx, r] == "background") {
            W_init_type[idx, r] <- "cross_trait"
          }
        }
      }
    }
  }
  
  W_list <- vector("list", K)
  for (k in seq_len(K)) {
    W_list[[k]] <- W_base
  }
  
  list(
    W_list = W_list,
    W_base = W_base,
    W_init_type = W_init_type
  )
}

prepare_fit_inputs <- function(X_list, Y_list, q, K) {
  pre <- precompute_XtX_XtY(X_list, Y_list)
  
  if (length(X_list) != K || length(Y_list) != K) {
    stop("X_list and Y_list must both have length K.")
  }
  
  ## Population-specific initialization:
  ## build W_k^(0) only from (X_k, Y_k), without pooling populations.
  W_list <- vector("list", K)
  W_init_type_list <- vector("list", K)
  
  for (k in seq_len(K)) {
    init_k <- initialize_stepwise_halo(
      X_list = list(X_list[[k]]),
      Y_list = list(Y_list[[k]]),
      q = q,
      K = 1L
    )
    
    W_list[[k]] <- init_k$W_base
    W_init_type_list[[k]] <- init_k$W_init_type
  }
  
  init <- list(
    W_list = W_list,
    W_base_list = W_list,
    W_init_type_list = W_init_type_list
  )
  
  c(pre, init)
}
