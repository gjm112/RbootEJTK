
##########################################
# JTK_Cycle 
##########################################
# # install 'devtools' in R(>3.0.2)
# install.packages("devtools")
# # install MetaCycle
# devtools::install_github('gangwug/MetaCycle')
library(tidyverse)
x <- seq(0,48,1)
#y <- cos(2*pi*x/12) + rnorm(length(x),0,10)
asymmetry <- 6 
phase_shift <- 2
phase <- 8
x <- (x + phase_shift)%%phase
y <- cos(2*pi*(x - phase_shift)/(2*asymmetry)) * ((x - phase_shift) <= asymmetry) + 
  cos(2*pi*(x - phase_shift - phase)/(2*(phase - asymmetry))) * ((x - phase_shift) > asymmetry) + 
  rnorm(length(x), 0, 0)
x <- seq(0,48,1) + phase_shift

#y <- rcauchy(length(x))
plot(x, y, type = "l", xlim = c(0, 10))
abline(v = 2)
plot(x, y, type = "l", xlim = c(12, 24))
plot(x, y, type = "l", xlim = c(24, 36))

dat <- data.frame(time = x, value = y)
timepoints <- dat$time
values <- dat$value
#Function takes in a times series 
jtk1 <- function(timepoints, values, phase = 12, phase_shift = 0, asymmetry = 0){
  t <- timepoints%% phase
    #Reference function
    valuesref <- cos(2*pi*(t - phase_shift)/(2*asymmetry)) * (t <= asymmetry) + 
      cos(2*pi*(t - phase_shift - phase)/(2*(phase - asymmetry))) * (t > asymmetry)
    #return(valuesref))
    
    #Compute Kendall's tau
    tau <- cor(values,valuesref, method = "kendall")
    pval <- cor.test(values,valuesref, method = "kendall", alternative = "greater")$p.value
    return(c(tau, pval))
}

jtk1(x, y, 12, 0, 4)

plot(x, jtk1(x, y, 12, 0, 4), type = "l")
points(x,y, type = "l")

#1. Turn this into a function
phase_vec <- seq(4,20,4)
phase_shift_vec <- seq(0,12,4)
asymmetry_vec <- seq(2,12,2)

jtk_n <- function(timepoints, values, phase_vec, phase_shift_vec, asymmetry_vec){
  results <- data.frame()
  for (p in phase_vec){
    for (s in phase_shift_vec){
      for (a in asymmetry_vec){
        if (a < p){
      temp <- jtk1(timepoints, values, p, s, a)
      results <- bind_rows(results, data.frame(phase = p, shift = s , asymmetry = a, tau = temp[1], pval = temp[2]))
        }
      }
    }
  }
  results <- results |> arrange(pval)
  return(head(results,1))
}

jtk_n(timepoints, values, phase_vec, phase_shift_vec, asymmetry_vec)

#Permute the values and calculate a null distribution.  
eJTK <- function(timepoints, values, phase_vec, phase_shift_vec, asymmetry_vec, nsim = 1000){
  null <- data.frame()
  for (i in 1:nsim){
    temp <- jtk_n(timepoints, 
          sample(values, size = length(values), replace = F), 
          phase_vec,
          phase_shift_vec,
          asymmetry_vec)
    null <- bind_rows(null, temp)
  }
  
  actual <- jtk_n(timepoints, values, phase_vec, phase_shift_vec, asymmetry_vec)
  p_val <- sum(null$tau >= actual$tau) / nsim
  return(data.frame(phase = actual$phase, shift = actual$shift, asymmetry = actual$asymmetry, pval = p_val))
}
  
eJTK(timepoints, values, phase_vec, phase_shift_vec, asymmetry_vec)

#bootJTK
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
# The following initializes usage of Bioc devel
#BiocManager::install(version='devel')
#BiocManager::install("limma")
library(limma)

#Step 1 of BooteJTK
reps <- 4
x <- seq(0,24,1)
y <- matrix(NA, nrow = 4, ncol = length(x))
for (q in 1:reps){
x <- seq(0,24,1)
#y <- cos(2*pi*x/12) + rnorm(length(x),0,10)
asymmetry <- 4
phase_shift <- 4
phase <- 8
x <- (x - phase_shift) %% phase 
y[q,] <- cos(2*pi*(x - phase_shift)/(2*asymmetry)) * (x <= asymmetry) + 
  cos(2*pi*(x - phase_shift - phase)/(2*(phase - asymmetry))) * (x > asymmetry) + 
  rnorm(length(x), 0, 0)
x <- seq(0,24,1)
}

plot(x,y[1,], type = "l", xlim = c(0,8))


#y <- matrix(rpois(reps*length(x),10), ncol = reps)
xbar <- apply(y, 1, mean)
greg <- vooma(y, plot = TRUE)
#variances
1/greg$weights[,1] 
#standard errors
sehat <- sqrt(1/greg$weights[,1])/sqrt(reps)
test <- vash(sehat,reps)
sd.post <- test$sd.post

#Now bootstrap
boot <- rnorm(length(x),xbar, sd.post)
boot

fuck <- eJTK(x,boot,phase_vec, phase_shift_vec, asymmetry_vec)

plot(test$sd.post,sehat, asp = 1)
abline(a = 0, b = 1)


#boostrap






#100 time points, 4 replicates
reps <- 1000

y <- matrix(rpois(reps*100,100), ncol = reps)
greg <- vooma(y, plot = TRUE)
plot(cbind(apply(y,1,function(x){var(x)}),(1/greg$weights[,1])))
apply(y,1,function(x){var(x)})
(1/greg$weights[,1])


cbind(apply(y,1,function(x){var(x)}),(1/greg$weights[,1]))
plot(cbind(apply(y,1,function(x){var(x)}),(1/greg$weights[,1])))


(cbind(1/greg$weights[,1],apply((y),1,var)))
abline(a = 0, b = 1)
abline(v = 100)







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


