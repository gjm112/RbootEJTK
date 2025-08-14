##########################################
####Data Prep
##########################################
library(tidyverse)
library(DiscoRhythm)
test <- as.data.frame(t(read.csv("/Users/gregorymatthews/Dropbox/fruit-flies-circadian-rhythms//MX7590~2(data)_clean_greg.csv", header = FALSE)))
names(test) <- test[1,]
test <- test[-1,]

#table(test$treatment)

#Converting to numeric
test[,9:ncol(test)] <- apply(test[,9:ncol(test)],2,as.numeric)

test <- test %>% filter(treatment != "pool - pool")
test <- 
  test %>% mutate(
    replicate = substr(treatment, nchar(treatment), nchar(treatment)),
    trt = substr(treatment, 1, nchar(treatment) - 4),
    index = 1:192,
    time = rep(seq(0, 46, 2), 8),
    c_time = rep(seq(0, 22, 2), 16),
    day = rep(rep(c(1, 2), each = 12), 8)
  )
test$trt[test$trt == "clock856"] <- "clk856"
names(test) <- gsub(" ", "", names(test))


#Pull out ONE 
sub <- test %>% select(treatment, zymostenol, replicate, trt, index, time, c_time, day) %>% filter(trt == "iso") %>% group_by(time) %>% summarize(values = mean(zymostenol), c_time = head(c_time,1)) 
sub 



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
phase <- p <- 10
phase_shift <- s <- 4
asymmetry <- a <- 6
y <- (-1/a * (x-s) %% p + 1)*((x-s) %% p < a) + (1/(p-a) * (x-s) %% p - (a/(p-a)))*((x-s) %% p > a) + rnorm(length(x), 0, 0)
plot(x, y, type = "l")


#These get passed to the functions
dat <- data.frame(time = x, value = y)
timepoints <- dat$time
values <- dat$value

timepoints <- sub$time
values <- sub$values




source("./jtk1.R")
source("./jtk_n.R")
source("./eJTK.R")
#1. Turn this into a function
#As I understand it, booteJTK requires three parameters: period length, asymmetry and phase. We may need to play around with these but a good stating place would be
#Period: 24 hrs
#Phase: every 2 hours from 0-22 (12 total)
#Asymmetry: every 2 hours from 2-22 (11 total)
phase_vec <- 24
phase_shift_vec <- seq(0,22,2)
asymmetry_vec <- seq(2,22,2)
#jtk1(x, y, 16, 0, 6)
#jtk_n(timepoints, values, phase_vec, phase_shift_vec, asymmetry_vec)
eJTK(timepoints, values, phase_vec, phase_shift_vec, asymmetry_vec)

#bootJTK
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
# The following initializes usage of Bioc devel
#BiocManager::install(version='devel')
#BiocManager::install("limma")
library(limma)

#Step 1 of BooteJTK
#Take the real data and average it across reps
reps <- 4
x <- seq(0,24,1)
y <- matrix(NA, nrow = reps, ncol = length(x))
for (q in 1:reps){
x <- seq(0,24,1)
#y <- cos(2*pi*x/12) + rnorm(length(x),0,10)
a <- asymmetry <- 4
s <- phase_shift <- 4
p <- phase <- 8
y[q,] <- (-1/a * (x-s) %% p + 1)*((x-s) %% p < a) + (1/(p-a) * (x-s) %% p - (a/(p-a)))*((x-s) %% p > a) + rnorm(length(x), 0, 10)
}

plot(x,y[1,], type = "l", xlim = c(0,8))

#y <- matrix(rpois(reps*length(x),10), ncol = reps)
xbar <- apply(y, 2, mean)
voom <- vooma(t(y), plot = TRUE)
#variances
1/voom$weights[,1] 
#standard errors
sehat <- sqrt(1/voom$weights[,1])/sqrt(reps)
test <- vash(sehat,reps)
sd.post <- test$sd.post
#cbind(sd.post, sehat)
#Now bootstrap
boot_results <- list()
for (j in 1:10){print(j)
boot <- rnorm(length(x),xbar, sd.post)
boot_results[[j]] <- eJTK(x,boot,phase_vec, phase_shift_vec, asymmetry_vec)
}

boot_results_df <- do.call(rbind,boot_results)
boot_results_df
apply(boot_results_df[,1:4],2,mean)

# plot(test$sd.post,sehat, asp = 1)
# abline(a = 0, b = 1)


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


