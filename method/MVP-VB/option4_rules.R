get_gate_rules <- function(K, threshold) {
  option2_current_z_min <- 0.75
  
  list(
    gate_quantile = 0.90,
    local_ok_rule = "either",
    c_W_gate = threshold,
    
    ## Option 2: same SNP-same phenotype replication across populations.
    option2_alpha = 0.10,
    
    ## When K=2, LOO meta evidence is essentially one external population.
    ## Therefore do not make the one-external-population rule stricter than option2_alpha.
    option2_meta_alpha = if (K <= 2) 0.10 else 0.05,
    
    ## Require at least very weak current-population support before borrowing.
    ## |z| > 0.75 corresponds to two-sided p approximately 0.45.
    option2_current_z_min = option2_current_z_min,
    option2_current_alpha = 2 * pnorm(-option2_current_z_min),
    
    ## Opposite-direction guard.
    ## Use 0.05 rather than 0.10 to avoid blocking borrowing due to weak noisy opposite signs.
    option2_oppose_alpha = 0.05,
    
    ## Option 3: same SNP cross-phenotype pleiotropy.
    option3_other_alpha = 1e-4,
    option3_current_alpha = 0.10,
    
    ## One other phenotype with strong evidence is enough to define pleiotropy.
    min_other_signal = 1L,
    
    option3_margin = 0.05
  )
}
