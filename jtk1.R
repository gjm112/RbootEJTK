jtk1 <- function(timepoints, values, phase = 12, phase_shift = 0, asymmetry = 0){
  
  #Reference function
  valuesref <-  (-1/asymmetry * (timepoints-phase_shift) %% phase + 1)*((timepoints-phase_shift) %% phase < asymmetry) + (1/(phase-asymmetry) * (timepoints-phase_shift) %% phase - (asymmetry/(phase-asymmetry)))*((timepoints-phase_shift) %% phase > asymmetry)
  #return(valuesref))
  
  #Compute Kendall's tau
  tau <- cor(values,valuesref, method = "kendall")
  pval <- cor.test(values,valuesref, method = "kendall", alternative = "greater")$p.value
  return(c(tau, pval))
}
