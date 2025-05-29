##########################################
# JTK_Cycle 
##########################################
# # install 'devtools' in R(>3.0.2)
# install.packages("devtools")
# # install MetaCycle
# devtools::install_github('gangwug/MetaCycle')
library(tidyverse)
x <- seq(0,46,1)
y <- cos(2*pi*x/10) + rnorm(length(x),0,1)
plot(x, y, type = "l")

dat <- data.frame(time = x, value = y)
timepoints <- dat$time
values <- dat$value
#Function takes in a times series 
jtk1 <- function(timepoints, values, phase = 12, phase_shift = 0){
    #Reference function  
    valuesref <- cos(2*pi*(timepoints - phase_shift)/phase)
    
    #Compute Kendall's tau
    tau <- cor(values,valuesref, method = "kendall")
    pval <- cor.test(values,valuesref, method = "kendall", alternative = "greater")$p.value
    return(c(tau, pval))
}

#1. Turn this into a function
phase_vec <- seq(4,20,1)
phase_shift_vec <- seq(0,8,4)
results <- data.frame()
for (p in phase_vec){
  for (s in phase_shift_vec){
  temp <- jtk1(timepoints, values, p, s)
  results <- bind_rows(results, data.frame(phase = p, shift = s , tau = temp[1], pval = temp[2]))
  }
}

results %>% arrange(pval)

#Permute the values and calculate a null distribution.  
eJTK <- function(timepoints, values, phase_vec, phase_shift_vec, nsim = 10){
    out <- replicate(nsim,jtk1(timepoints, 
       values, 
       sample(phase_vec,1), 
       sample(phase_shift_vec,1)))
  out <- data.frame(t(out))
  names(out) <- c("tau","pval")
    return(out)
}
  
eJTK <- function(timepoints, values, phase_vec, phase_shift_vec, nsim = 10){
  out <- replicate(nsim,jtk1(timepoints, 
       values, 
       sample(phase_vec,1), 
       sample(phase_shift_vec,1)))
  out <- data.frame(t(out))
  names(out) <- c("tau","pval")
    return(out)
}

#eJTK
ejtk_null <- eJTK(timepoints, 
                  values, 
                  phase_vec = seq(4,20,4), 
                  phase_shift_vec = seq(0,8,4), 
                  nsim = 1000000)

hist(ejtk_null$pval)
results %>% arrange(pval)

for (i in 1:nrow(results)){
  results$epval[i] <- mean(results$pval[i] >= ejtk_null$pval)
}

results %>% arrange(epval)





#vashr
#install.packages("devtools")
library(devtools)
install_github("mengyin/vashr",build_vignettes=TRUE)

library(vashr)

browseVignettes("vashr")
#https://pmc.ncbi.nlm.nih.gov/articles/PMC5181563/pdf/btw483.pdf

#Metacycle: https://github.com/gangwug/MetaCycle


