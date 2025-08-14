#Permute the values and calculate a null distribution.  
eJTK <- function(timepoints, values, phase_vec, phase_shift_vec, asymmetry_vec, shape, rate){
  #alpha_hat <- mean(null$pval)^2 / var(null$pval)
  #beta_hat <- var(null$pval) / mean(null$pval)
  
  
  
  
  # for (i in 1:5){print(i)
  #   temp <- jtk_n(timepoints, 
  #                 sample(values, size = length(values), replace = F), 
  #                 phase_vec,
  #                 phase_shift_vec,
  #                 asymmetry_vec)
  #   null <- bind_rows(null, temp)
  # }
  
  observed <- jtk_n(timepoints, values, phase_vec, phase_shift_vec, asymmetry_vec)
  if (observed$tau <= 0){p_val <- 1} else {
  p_val <- pgamma(observed$tau, shape = shape,rate = rate, lower.tail = FALSE)}
  #p_val <- sum(null$tau >= observed$tau) / nsim
  return(data.frame(tau = observed$tau, phase = observed$phase, shift = observed$shift, asymmetry = observed$asymmetry, pval = p_val))
}
