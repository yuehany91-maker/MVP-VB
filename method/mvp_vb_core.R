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
      
