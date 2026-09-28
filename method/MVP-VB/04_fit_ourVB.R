fit_ourVB <- function(
    common_data,
    scenario_name,
    max_iter = 200,
    gate_threshold = 0.5,
    active_controls = NULL,
    core_script = Sys.getenv(
      "CORE_SCRIPT",
      unset = file.path(
        getwd(),
        "04_fit_ourVB_core_dense_exact_from_current.R"
      )
    )
) {
  
  core_script <- normalizePath(
    core_script,
    mustWork = FALSE
  )
  
  if (!file.exists(core_script)) {
    stop(
      "Cannot find static VB core script: ",
      core_script
    )
  }
  
  cat("ourVB core script:", core_script, "\n")
  
  fit_controls <- list(
    max_iter = max_iter,
    gate_threshold = gate_threshold,
    active_controls = active_controls
  )
  
  env <- new.env(parent = globalenv())
  
  env$fit_controls <- fit_controls
  env$X_list <- common_data$X_list
  env$Y_list <- common_data$Y_list
  env$Sigma_true <- common_data$Sigma_true
  env$n_k <- common_data$n_k
  env$K <- common_data$K
  env$m <- common_data$m
  env$q <- common_data$q
  env$snp_key_vec <- common_data$snp_key_vec
  env$block_id_vec <- common_data$block_id_vec
  
  prepared <- prepare_fit_inputs(
    X_list = common_data$X_list,
    Y_list = common_data$Y_list,
    q = common_data$q,
    K = common_data$K
  )
  
  list2env(prepared, envir = env)
  
  env$sigmaB2 <- 5.0
  env$pi_fix <- 0.005
  env$nu0 <- common_data$q + 2
  env$Psi0 <- diag(common_data$q)
  
  runtime_sec <- system.time({
    source(
      core_script,
      local = env,
      encoding = "UTF-8"
    )
  })[["elapsed"]]
  
  n_iter_done <- env$n_iter_done
  iter_summary <- env$iter_summary
  
  if (
    !is.null(iter_summary) &&
    nrow(iter_summary) >= n_iter_done
  ) {
    iter_summary <- iter_summary[
      seq_len(n_iter_done),
      ,
      drop = FALSE
    ]
  }
  
  converged <- isTRUE(n_iter_done < max_iter)
  
  list(
    method = "ourVB",
    scenario_name = scenario_name,
    config = common_data$config,
    W_list = env$W_list,
    M_list = env$M_list,
    Tau_list = env$Tau_list,
    option_array = env$option_array,
    iter_summary = iter_summary,
    n_iter_done = n_iter_done,
    converged = converged,
    runtime_sec = runtime_sec,
    gate_threshold = gate_threshold,
    
    active_idx_list =
      if (exists(
        "active_idx_list",
        envir = env,
        inherits = FALSE
      )) {
        env$active_idx_list
      } else {
        NULL
      },
    
    C_list =
      if (exists(
        "C_list",
        envir = env,
        inherits = FALSE
      )) {
        env$C_list
      } else {
        NULL
      },
    
    Z_watch_list =
      if (exists(
        "Z_watch_list",
        envir = env,
        inherits = FALSE
      )) {
        env$Z_watch_list
      } else {
        NULL
      },
    
    active_controls =
      if (exists(
        "active_controls",
        envir = env,
        inherits = FALSE
      )) {
        env$active_controls
      } else {
        active_controls
      },
    
    c_target_k =
      if (exists(
        "c_target_k",
        envir = env,
        inherits = FALSE
      )) {
        env$c_target_k
      } else {
        NULL
      },
    
    c_z_k =
      if (exists(
        "c_z_k",
        envir = env,
        inherits = FALSE
      )) {
        env$c_z_k
      } else {
        NULL
      },
    
    c_w_k =
      if (exists(
        "c_w_k",
        envir = env,
        inherits = FALSE
      )) {
        env$c_w_k
      } else {
        NULL
      },
    
    c_ld_k =
      if (exists(
        "c_ld_k",
        envir = env,
        inherits = FALSE
      )) {
        env$c_ld_k
      } else {
        NULL
      },
    
    active_size_history =
      if (exists(
        "active_size_history",
        envir = env,
        inherits = FALSE
      )) {
        env$active_size_history
      } else {
        NULL
      },
    
    active_idx_history =
      if (exists(
        "active_idx_history",
        envir = env,
        inherits = FALSE
      )) {
        env$active_idx_history
      } else {
        NULL
      }
  )
}