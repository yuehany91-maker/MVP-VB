## ============================================================
## Pair-level elementwise multivariate Lasso with EBIC tuning
##
## Model for each population k:
##   min_B  0.5 * ||Y - X B||_F^2 + lambda * sum_{j,r} |B[j,r]|
##
## Group penalty is completely disabled:
##   Pen.G = 0 and lam.G = 0.
##
## EBIC uses the pair-level model space:
##   N = n * q  scalar responses
##   P = m * q  candidate SNP-phenotype coefficients
##   d = number of nonzero estimated coefficients
##
##   EBIC_gamma(lambda)
##     = N * log(RSS_lambda / N)
##       + d * log(N)
##       + 2 * gamma * log choose(P, d).
##
## B_true is never used in fitting, lambda-grid construction, or EBIC.
## It is accessed only after fitting to compute simulation metrics.
##
## Required server entry point:
##   fit_MVLasso_from_common_data(common_data, scenario_name, paths,
##                                config, n_cores = 1)
##
## Returned metrics are population-specific only. No AUC and no
## cross-population average/overall row are computed.
## ============================================================

suppressPackageStartupMessages({
  library(MSGLasso)
})

## ------------------------------------------------------------
## Small helpers
## ------------------------------------------------------------

env_flag <- function(name, default = FALSE) {
  x <- tolower(trimws(Sys.getenv(
    name,
    unset = if (isTRUE(default)) "true" else "false"
  )))
  x %in% c("1", "true", "t", "yes", "y")
}

validate_scalar_integer <- function(x, name, lower = 1L, upper = Inf) {
  if (length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != floor(x) || x < lower || x > upper) {
    stop(name, " must be one integer in [", lower, ", ", upper, "].")
  }
  as.integer(x)
}

## ------------------------------------------------------------
## Corrected MSGLasso 2.1 low-level wrapper
##
## This keeps the package algorithm and compiled C routine, but fixes
## the R-wrapper naming bug involving Pen.L/Pen.G versus Pen_L/Pen_G.
## ------------------------------------------------------------

MSGLasso_fixed <- function(X.m,
                           Y.m,
                           grp.WTs,
                           Pen.L,
                           Pen.G,
                           PQ.grps,
                           GR.grps,
                           grp_Norm0,
                           lam1,
                           lam.G,
                           Beta0 = NULL) {
  X.m <- as.matrix(X.m)
  Y.m <- as.matrix(Y.m)
  grp.WTs <- as.matrix(grp.WTs)
  Pen.L <- as.matrix(Pen.L)
  Pen.G <- as.matrix(Pen.G)
  PQ.grps <- as.matrix(PQ.grps)
  GR.grps <- as.matrix(GR.grps)
  grp_Norm0 <- as.matrix(grp_Norm0)
  lam.G <- as.matrix(lam.G)
  
  N <- nrow(X.m)
  P <- ncol(X.m)
  Q <- ncol(Y.m)
  G <- nrow(grp.WTs)
  R <- ncol(grp.WTs)
  
  X.v <- as.vector(t(X.m))
  Y.v <- as.vector(t(Y.m))
  grpWTs <- as.vector(t(grp.WTs))
  Pen_L <- as.vector(t(Pen.L))
  Pen_G <- as.vector(t(Pen.G))
  
  gmax <- ncol(PQ.grps)
  if (is.null(gmax) || gmax <= 1L) {
    stop("PQ.grps must have at least two columns.")
  }
  PQgrps <- as.vector(t(PQ.grps))
  
  cmax <- ncol(GR.grps)
  if (is.null(cmax) || cmax <= 1L) {
    stop("GR.grps must have at least two columns.")
  }
  GRgrps <- as.vector(t(GR.grps))
  
  lambda1 <- as.numeric(lam1)
  lambdaG <- as.vector(t(lam.G))
  
  Beta.m <- matrix(0, nrow = P, ncol = Q)
  Beta.v <- as.vector(Beta.m)
  grp_Norm.v <- as.vector(t(grp_Norm0))
  
  iter.count <- 0L
  RSS <- 0
  Edebug <- rep(0, N * Q)
  
  if (is.null(Beta0)) {
    junk <- .C(
      "MSGLasso",
      as.integer(N),
      as.integer(P),
      as.integer(Q),
      as.integer(G),
      as.integer(R),
      as.double(X.v),
      as.double(Y.v),
      as.double(grpWTs),
      as.integer(Pen_L),
      as.integer(Pen_G),
      as.integer(gmax),
      as.integer(PQgrps),
      as.integer(cmax),
      as.integer(GRgrps),
      as.double(lambda1),
      as.double(lambdaG),
      grp_Norm = as.double(grp_Norm.v),
      Beta.out = as.double(Beta.v),
      n.iter = as.integer(iter.count),
      RSS = as.double(RSS),
      Edebug = as.double(Edebug),
      PACKAGE = "MSGLasso"
    )
  } else {
    Beta0 <- as.matrix(Beta0)
    if (!identical(dim(Beta0), c(P, Q))) {
      stop("Beta0 has the wrong dimension.")
    }
    
    Beta.ini.v <- as.vector(t(Beta0))
    
    junk <- .C(
      "MSGLasso_Ini",
      as.integer(N),
      as.integer(P),
      as.integer(Q),
      as.integer(G),
      as.integer(R),
      as.double(X.v),
      as.double(Y.v),
      as.double(grpWTs),
      as.integer(Pen_L),
      as.integer(Pen_G),
      as.integer(gmax),
      as.integer(PQgrps),
      as.integer(cmax),
      as.integer(GRgrps),
      as.double(lambda1),
      as.double(lambdaG),
      Beta.ini = as.double(Beta.ini.v),
      grp_Norm = as.double(grp_Norm.v),
      Beta.out = as.double(Beta.v),
      n.iter = as.integer(iter.count),
      RSS = as.double(RSS),
      Edebug = as.double(Edebug),
      PACKAGE = "MSGLasso"
    )
  }
  
  Beta.result <- matrix(junk$Beta.out, nrow = P, byrow = TRUE)
  grp_Norm.result <- matrix(junk$grp_Norm, nrow = G, byrow = TRUE)
  E <- matrix(junk$Edebug, nrow = N, ncol = Q, byrow = TRUE)
  
  rss.v <- colSums(E^2)
  rss <- sum(rss.v)
  
  list(
    Beta = Beta.result,
    grpNorm = grp_Norm.result,
    E = E,
    rss.v = rss.v,
    rss = rss,
    iter = as.integer(junk$n.iter)
  )
}

## ------------------------------------------------------------
## Dummy singleton group metadata required by the package C API.
##
## The metadata does not create an active group penalty because both
## Pen.G and lam.G are set to zero in every fit.
## ------------------------------------------------------------

build_MVLasso_dummy_groups <- function(m, q) {
  G <- m
  R <- q
  gmax <- 1L
  cmax <- max(m, q)
  
  G.Starts <- 0:(m - 1L)
  G.Ends <- 0:(m - 1L)
  R.Starts <- 0:(q - 1L)
  R.Ends <- 0:(q - 1L)
  
  PQgrps <- MSGLasso::FindingPQGrps(
    m,
    q,
    G,
    R,
    gmax,
    G.Starts,
    G.Ends,
    R.Starts,
    R.Ends
  )$PQgrps
  
  GRgrps <- MSGLasso::FindingGRGrps(
    m,
    q,
    G,
    R,
    cmax,
    G.Starts,
    G.Ends,
    R.Starts,
    R.Ends
  )$GRgrps
  
  grpWTs <- MSGLasso::Cal_grpWTs(
    m,
    q,
    G,
    R,
    gmax,
    PQgrps
  )$grpWTs
  
  list(
    G = G,
    R = R,
    gmax = gmax,
    cmax = cmax,
    PQgrps = PQgrps,
    GRgrps = GRgrps,
    grpWTs = grpWTs
  )
}

## ------------------------------------------------------------
## Standardization
##
## X columns: centered and scaled to sample SD 1.
## Y columns: centered only.
## No true-effect information is used.
## ------------------------------------------------------------

standardize_for_MVLasso <- function(X, Y) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  if (nrow(X) != nrow(Y)) {
    stop("X and Y must have the same number of rows.")
  }
  
  if (any(!is.finite(X))) {
    stop("X contains a non-finite value.")
  }
  
  if (any(!is.finite(Y))) {
    stop("Y contains a non-finite value.")
  }
  
  x_center <- colMeans(X)
  Xs <- sweep(X, 2, x_center, "-")
  
  x_scale <- apply(Xs, 2, sd)
  bad_x <- !is.finite(x_scale) | x_scale <= 0
  
  if (any(bad_x)) {
    warning(
      sum(bad_x),
      " constant X columns were set to zero after centering."
    )
    x_scale[bad_x] <- 1
  }
  
  Xs <- sweep(Xs, 2, x_scale, "/")
  Xs[, bad_x] <- 0
  
  y_center <- colMeans(Y)
  Ys <- sweep(Y, 2, y_center, "-")
  
  list(
    X = Xs,
    Y = Ys,
    x_center = x_center,
    x_scale = x_scale,
    y_center = y_center,
    constant_x = which(bad_x)
  )
}

## ------------------------------------------------------------
## Pair-level EBIC
## ------------------------------------------------------------

compute_pair_EBIC <- function(rss, n, m, q, d, gamma = 1) {
  N_eff <- as.double(n) * as.double(q)
  P_eff <- as.double(m) * as.double(q)
  
  if (d < 0 || d > P_eff) {
    stop("d must be between 0 and m*q.")
  }
  
  rss_safe <- max(as.numeric(rss), .Machine$double.xmin)
  
  ebic <- N_eff * log(rss_safe / N_eff) +
    as.double(d) * log(N_eff) +
    2 * as.numeric(gamma) * lchoose(P_eff, as.double(d))
  
  list(
    EBIC = as.numeric(ebic),
    N_EBIC = N_eff,
    P_EBIC = P_eff
  )
}

## ------------------------------------------------------------
## Native pair-level selection metrics
## ------------------------------------------------------------

make_MVLasso_metric_row <- function(scenario_name,
                                    population,
                                    selected,
                                    truth,
                                    lambda,
                                    ebic,
                                    gamma,
                                    N_EBIC,
                                    P_EBIC,
                                    d_lambda) {
  pred <- as.logical(as.vector(selected))
  true <- as.logical(as.vector(truth))
  
  if (length(pred) != length(true)) {
    stop("selected and truth have different lengths.")
  }
  
  if (anyNA(pred) || anyNA(true)) {
    stop("selected and truth must not contain NA values.")
  }
  
  TP <- sum(pred & true)
  FP <- sum(pred & !true)
  TN <- sum(!pred & !true)
  FN <- sum(!pred & true)
  
  Precision <- if (TP + FP > 0L) TP / (TP + FP) else 0
  Recall <- if (TP + FN > 0L) TP / (TP + FN) else NA_real_
  Specificity <- if (TN + FP > 0L) TN / (TN + FP) else NA_real_
  
  F1 <- if (is.finite(Recall) && Precision + Recall > 0) {
    2 * Precision * Recall / (Precision + Recall)
  } else {
    0
  }
  
  data.frame(
    scenario_name = scenario_name,
    population = population,
    method = "MVLasso_EBIC",
    tune_method = "pair_level_EBIC",
    lambda = lambda,
    EBIC = ebic,
    EBIC_gamma = gamma,
    N_EBIC = N_EBIC,
    P_EBIC = P_EBIC,
    d_lambda = d_lambda,
    TP = TP,
    FP = FP,
    TN = TN,
    FN = FN,
    Precision = Precision,
    Recall = Recall,
    Specificity = Specificity,
    F1 = F1,
    n_selected = TP + FP,
    n_true = TP + FN,
    stringsAsFactors = FALSE
  )
}

## ------------------------------------------------------------
## Fit one population
## ------------------------------------------------------------

fit_MVLasso_one_population <- function(X,
                                       Y,
                                       gamma = 1,
                                       grid_len = 50L,
                                       lambda_min_ratio = 1e-3,
                                       lambda_max_multiplier = 1.05,
                                       max_lambda_expansions = 30L,
                                       beta_tol = 1e-12,
                                       verbose = TRUE) {
  std <- standardize_for_MVLasso(X, Y)
  Xs <- std$X
  Ys <- std$Y
  
  n <- nrow(Xs)
  m <- ncol(Xs)
  q <- ncol(Ys)
  
  grid_len <- validate_scalar_integer(
    grid_len,
    "grid_len",
    lower = 2L,
    upper = 10000L
  )
  
  max_lambda_expansions <- validate_scalar_integer(
    max_lambda_expansions,
    "max_lambda_expansions",
    lower = 0L,
    upper = 100L
  )
  
  if (!is.finite(gamma) || gamma < 0 || gamma > 1) {
    stop("gamma must be in [0, 1].")
  }
  
  if (!is.finite(lambda_min_ratio) ||
      lambda_min_ratio <= 0 || lambda_min_ratio >= 1) {
    stop("lambda_min_ratio must be in (0, 1).")
  }
  
  if (!is.finite(lambda_max_multiplier) || lambda_max_multiplier <= 0) {
    stop("lambda_max_multiplier must be positive.")
  }
  
  if (!is.finite(beta_tol) || beta_tol < 0) {
    stop("beta_tol must be nonnegative.")
  }
  
  groups <- build_MVLasso_dummy_groups(m = m, q = q)
  
  grp.WTs <- groups$grpWTs
  PQ.grps <- groups$PQgrps
  GR.grps <- groups$GRgrps
  
  ## Elementwise L1 penalty is active for every pair.
  Pen.L <- matrix(1L, nrow = m, ncol = q)
  
  ## Group penalty is fully disabled.
  Pen.G <- matrix(0L, nrow = groups$G, ncol = groups$R)
  lam.G <- matrix(0, nrow = groups$G, ncol = groups$R)
  grp_Norm0 <- matrix(0, nrow = groups$G, ncol = groups$R)
  
  fit_at_lambda <- function(lambda,
                            Beta0 = NULL,
                            grp_norm_start = grp_Norm0) {
    MSGLasso_fixed(
      X.m = Xs,
      Y.m = Ys,
      grp.WTs = grp.WTs,
      Pen.L = Pen.L,
      Pen.G = Pen.G,
      PQ.grps = PQ.grps,
      GR.grps = GR.grps,
      grp_Norm0 = grp_norm_start,
      lam1 = lambda,
      lam.G = lam.G,
      Beta0 = Beta0
    )
  }
  
  ## Initial scale uses only observed X and Y.
  lambda_max <- max(abs(crossprod(Xs, Ys) / n), na.rm = TRUE)
  
  if (!is.finite(lambda_max) || lambda_max <= 0) {
    lambda_max <- 1
  }
  
  lambda_max <- lambda_max * lambda_max_multiplier
  
  ## Expand upward until the path includes the all-zero model.
  zero_fit <- NULL
  expansion_count <- 0L
  
  repeat {
    zero_fit <- fit_at_lambda(lambda_max)
    d_zero <- sum(abs(zero_fit$Beta) > beta_tol)
    
    if (d_zero == 0L) {
      break
    }
    
    if (expansion_count >= max_lambda_expansions) {
      stop(
        "Could not find an all-zero model after ",
        max_lambda_expansions,
        " lambda expansions."
      )
    }
    
    lambda_max <- lambda_max * 2
    expansion_count <- expansion_count + 1L
  }
  
  lambda_min <- lambda_max * lambda_min_ratio
  lambda_grid <- exp(seq(
    log(lambda_max),
    log(lambda_min),
    length.out = grid_len
  ))
  
  if (verbose) {
    cat("\nMVLasso pair-level EBIC tuning\n")
    cat("n =", n, "| m =", m, "| q =", q, "\n")
    cat("N_EBIC = n*q =", n * q, "\n")
    cat("P_EBIC = m*q =", m * q, "\n")
    cat("EBIC gamma =", gamma, "\n")
    cat("Group penalty disabled: Pen.G = 0 and lam.G = 0\n")
    cat("beta_tol =", beta_tol, "\n")
    cat("lambda_max expansions =", expansion_count, "\n")
    cat(
      "lambda range =",
      signif(max(lambda_grid), 6),
      "to",
      signif(min(lambda_grid), 6),
      "| grid length =",
      length(lambda_grid),
      "\n"
    )
  }
  
  ebic_rows <- vector("list", length(lambda_grid))
  
  Beta_prev <- NULL
  grp_prev <- grp_Norm0
  
  best_fit <- NULL
  best_idx <- NA_integer_
  best_ebic <- Inf
  
  for (ii in seq_along(lambda_grid)) {
    lambda_i <- lambda_grid[ii]
    
    fit_i <- fit_at_lambda(
      lambda = lambda_i,
      Beta0 = Beta_prev,
      grp_norm_start = grp_prev
    )
    
    if (any(!is.finite(fit_i$Beta))) {
      stop("Non-finite coefficient at lambda index ", ii, ".")
    }
    
    selected_i <- abs(fit_i$Beta) > beta_tol
    d_i <- sum(selected_i)
    rss_i <- as.numeric(fit_i$rss)
    
    ebic_i <- compute_pair_EBIC(
      rss = rss_i,
      n = n,
      m = m,
      q = q,
      d = d_i,
      gamma = gamma
    )
    
    ebic_rows[[ii]] <- data.frame(
      lambda_index = ii,
      lambda = lambda_i,
      RSS = rss_i,
      RSS_per_scalar_response = rss_i / ebic_i$N_EBIC,
      d_lambda = d_i,
      n = n,
      m = m,
      q = q,
      N_EBIC = ebic_i$N_EBIC,
      P_EBIC = ebic_i$P_EBIC,
      EBIC_gamma = gamma,
      EBIC = ebic_i$EBIC,
      iterations = fit_i$iter,
      stringsAsFactors = FALSE
    )
    
    ## Lambda grid is decreasing. Strict improvement keeps the larger,
    ## sparser lambda when EBIC values are exactly tied.
    if (ebic_i$EBIC < best_ebic) {
      best_ebic <- ebic_i$EBIC
      best_idx <- ii
      best_fit <- fit_i
    }
    
    Beta_prev <- fit_i$Beta
    grp_prev <- fit_i$grpNorm
  }
  
  ebic_table <- do.call(rbind, ebic_rows)
  ebic_table$selected <- seq_len(nrow(ebic_table)) == best_idx
  
  selected <- abs(best_fit$Beta) > beta_tol
  
  if (verbose) {
    cat("\nSelected pair-level EBIC model\n")
    cat("lambda =", signif(lambda_grid[best_idx], 8), "\n")
    cat("EBIC =", signif(best_ebic, 10), "\n")
    cat("nonzero pairs =", sum(selected), "\n")
    cat("RSS =", signif(best_fit$rss, 10), "\n")
  }
  
  list(
    fit = best_fit,
    Beta = as.matrix(best_fit$Beta),
    selected = matrix(selected, nrow = m, ncol = q),
    lambda = lambda_grid[best_idx],
    EBIC = best_ebic,
    RSS = best_fit$rss,
    d_lambda = sum(selected),
    gamma = gamma,
    N_EBIC = n * q,
    P_EBIC = m * q,
    beta_tol = beta_tol,
    lambda_grid = lambda_grid,
    ebic_table = ebic_table,
    lambda_max_expansions = expansion_count,
    standardization = std
  )
}

## ------------------------------------------------------------
## Validate X/Y before fitting. B_true is intentionally omitted here.
## ------------------------------------------------------------

validate_MVLasso_XY <- function(common_data) {
  X_list <- common_data$X_list
  Y_list <- common_data$Y_list
  
  if (is.null(X_list) || is.null(Y_list)) {
    stop("common_data must contain X_list and Y_list.")
  }
  
  if (length(X_list) < 1L || length(X_list) != length(Y_list)) {
    stop("X_list and Y_list must have the same positive length.")
  }
  
  m <- ncol(as.matrix(X_list[[1]]))
  q <- ncol(as.matrix(Y_list[[1]]))
  
  for (k in seq_along(X_list)) {
    Xk <- as.matrix(X_list[[k]])
    Yk <- as.matrix(Y_list[[k]])
    
    if (nrow(Xk) != nrow(Yk)) {
      stop("X/Y row mismatch in population ", k, ".")
    }
    
    if (ncol(Xk) != m || ncol(Yk) != q) {
      stop("All populations must share the same m and q.")
    }
  }
  
  invisible(list(K = length(X_list), m = m, q = q))
}

## ------------------------------------------------------------
## Server entry point
## ------------------------------------------------------------

fit_MVLasso_from_common_data <- function(common_data,
                                         scenario_name,
                                         paths,
                                         config,
                                         n_cores = 1) {
  dims <- validate_MVLasso_XY(common_data)
  
  K <- dims$K
  m <- dims$m
  q <- dims$q
  
  gamma <- as.numeric(Sys.getenv("MVLASSO_EBIC_GAMMA", unset = "1"))
  grid_len_raw <- suppressWarnings(as.numeric(Sys.getenv(
    "MVLASSO_GRID_LEN",
    unset = "50"
  )))
  grid_len <- validate_scalar_integer(
    grid_len_raw,
    "MVLASSO_GRID_LEN",
    lower = 2L,
    upper = 10000L
  )
  
  lambda_min_ratio <- as.numeric(Sys.getenv(
    "MVLASSO_LAMBDA_MIN_RATIO",
    unset = "1e-3"
  ))
  lambda_max_multiplier <- as.numeric(Sys.getenv(
    "MVLASSO_LAMBDA_MAX_MULTIPLIER",
    unset = "1.05"
  ))
  max_exp_raw <- suppressWarnings(as.numeric(Sys.getenv(
    "MVLASSO_MAX_LAMBDA_EXPANSIONS",
    unset = "30"
  )))
  max_lambda_expansions <- validate_scalar_integer(
    max_exp_raw,
    "MVLASSO_MAX_LAMBDA_EXPANSIONS",
    lower = 0L,
    upper = 100L
  )
  
  beta_tol <- as.numeric(Sys.getenv(
    "MVLASSO_BETA_TOL",
    unset = "1e-12"
  ))
  verbose <- env_flag("MVLASSO_VERBOSE", default = TRUE)
  
  if (!is.finite(gamma) || gamma < 0 || gamma > 1) {
    stop("MVLASSO_EBIC_GAMMA must be in [0, 1].")
  }
  
  if (!is.finite(beta_tol) || beta_tol < 0) {
    stop("MVLASSO_BETA_TOL must be nonnegative.")
  }
  
  cat("\n============================================================\n")
  cat("Running pair-level elementwise MVLasso with EBIC\n")
  cat("Scenario:", scenario_name, "\n")
  cat("K =", K, "| m =", m, "| q =", q, "\n")
  cat("EBIC: N = n*q, P = m*q, gamma =", gamma, "\n")
  cat("Group penalty: disabled\n")
  cat("AUC: not computed\n")
  cat("Population averaging: not computed\n")
  cat("============================================================\n")
  
  fit_list <- vector("list", K)
  B_hat_list <- vector("list", K)
  selected_list <- vector("list", K)
  tuning_rows <- vector("list", K)
  ebic_path_list <- vector("list", K)
  
  start_time <- Sys.time()
  
  ## Fitting stage: only X and Y are used in this entire loop.
  for (k in seq_len(K)) {
    cat("\n------------------------------------------------------------\n")
    cat("MVLasso population", k, "of", K, "\n")
    cat("------------------------------------------------------------\n")
    
    fit_k <- fit_MVLasso_one_population(
      X = common_data$X_list[[k]],
      Y = common_data$Y_list[[k]],
      gamma = gamma,
      grid_len = grid_len,
      lambda_min_ratio = lambda_min_ratio,
      lambda_max_multiplier = lambda_max_multiplier,
      max_lambda_expansions = max_lambda_expansions,
      beta_tol = beta_tol,
      verbose = verbose
    )
    
    if (!identical(dim(fit_k$Beta), c(m, q))) {
      stop("MVLasso coefficient dimension mismatch in population ", k, ".")
    }
    
    fit_list[[k]] <- fit_k
    B_hat_list[[k]] <- fit_k$Beta
    selected_list[[k]] <- 1L * fit_k$selected
    
    tuning_rows[[k]] <- data.frame(
      scenario_name = scenario_name,
      population = k,
      method = "MVLasso_EBIC",
      tune_method = "pair_level_EBIC",
      lambda = fit_k$lambda,
      EBIC = fit_k$EBIC,
      EBIC_gamma = fit_k$gamma,
      RSS = fit_k$RSS,
      d_lambda = fit_k$d_lambda,
      N_EBIC = fit_k$N_EBIC,
      P_EBIC = fit_k$P_EBIC,
      beta_tol = fit_k$beta_tol,
      lambda_max_expansions = fit_k$lambda_max_expansions,
      stringsAsFactors = FALSE
    )
    
    ebic_path_list[[k]] <- cbind(
      scenario_name = scenario_name,
      population = k,
      fit_k$ebic_table
    )
    
    invisible(gc())
  }
  
  end_time <- Sys.time()
  
  ## Evaluation stage: B_true is first accessed here, after every lambda
  ## and selected support has already been fixed.
  B_true_list <- common_data$B_true_list
  
  if (is.null(B_true_list) || length(B_true_list) != K) {
    stop("B_true_list is required for simulation metrics and must have length K.")
  }
  
  metric_rows <- vector("list", K)
  
  for (k in seq_len(K)) {
    B_true_k <- as.matrix(B_true_list[[k]])
    
    if (!identical(dim(B_true_k), c(m, q))) {
      stop("B_true dimension mismatch in population ", k, ".")
    }
    
    truth_k <- B_true_k != 0
    tune_k <- tuning_rows[[k]]
    
    metric_rows[[k]] <- make_MVLasso_metric_row(
      scenario_name = scenario_name,
      population = k,
      selected = selected_list[[k]],
      truth = truth_k,
      lambda = tune_k$lambda,
      ebic = tune_k$EBIC,
      gamma = tune_k$EBIC_gamma,
      N_EBIC = tune_k$N_EBIC,
      P_EBIC = tune_k$P_EBIC,
      d_lambda = tune_k$d_lambda
    )
  }
  
  metrics_by_population <- do.call(rbind, metric_rows)
  rownames(metrics_by_population) <- NULL
  
  tuning_table <- do.call(rbind, tuning_rows)
  rownames(tuning_table) <- NULL
  
  names(ebic_path_list) <- paste0("population_", seq_len(K))
  
  cat("\n--- MVLasso EBIC tuning by population ---\n")
  print(tuning_table)
  
  cat("\n--- MVLasso pair-level metrics by population ---\n")
  print(metrics_by_population)
  
  list(
    method = "MVLasso_EBIC",
    scenario_name = scenario_name,
    config = common_data$config,
    
    ## Compatibility field for an existing server pipeline.
    ## It is the native 0/1 selected support, not a continuous AUC score.
    score_list = selected_list,
    selected_list = selected_list,
    B_hat_list = B_hat_list,
    fit_list = fit_list,
    
    tuning_table = tuning_table,
    ebic_path_list = ebic_path_list,
    metrics_by_population = metrics_by_population,
    
    supports_auc = FALSE,
    has_overall_metrics = FALSE,
    tune_method = "pair_level_EBIC",
    group_penalty = FALSE,
    EBIC_definition = list(
      N = "n*q",
      P = "m*q",
      gamma = gamma,
      formula = paste0(
        "N*log(RSS/N) + d*log(N) + ",
        "2*gamma*lchoose(P,d)"
      )
    ),
    beta_tol = beta_tol,
    runtime_sec = as.numeric(difftime(
      end_time,
      start_time,
      units = "secs"
    ))
  )
}

