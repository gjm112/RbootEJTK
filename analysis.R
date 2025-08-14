source("./jtk1.R")
source("./jtk_n.R")
source("./eJTK.R")
library(tidyverse)
library(limma)
library(DiscoRhythm)
library(vashr)

all <- data.frame()
##########################################
####Data Prep
##########################################
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

phase_vec <- 24
phase_shift_vec <- seq(0,22,12)
asymmetry_vec <- seq(2,22,12)

#All the markers to test
markers <- test %>% select(zymostenol:`134`) %>% names()
test_long <- test |> pivot_longer(cols = zymostenol:`134`, names_to = "markers", values_to = "values")


#Pull out ONE 
for (i in 1:1){
  for (t in c("clk856","iso","to","trf")){
    sub <- test_long  %>% select(treatment,
                                 markers,
                                 values,
                                 replicate,
                                 trt,
                                 index,
                                 time,
                                 c_time,
                                 day) %>% filter(trt == t &  markers == markers[i])


##########################################
# JTK_Cycle 
##########################################
    
timepoints <- sub %>% group_by(time) %>% summarize(values = mean(values)) %>% pull(time)
values <- sub %>% group_by(time) %>% summarize(values = mean(values)) %>% pull(values)

#Step 1 of BooteJTK
#Take the real data and average it across reps
reps <- 2
xbar <- sub %>% group_by(time) %>% summarize(mean = mean(values)) %>% pull(mean)
voom <- vooma(matrix(sub$values, ncol = 2))
#variances
#1/voom$weights[,1] 
#standard errors
sehat <- sqrt(1/voom$weights[,1])/sqrt(reps)
vashout <- vash(sehat,reps)
sd.post <- vashout$sd.post
#cbind(sd.post, sehat)
#Now bootstrap
boot_results <- function(xbar, sd.post, timepoints){
  boot <- rnorm(length(xbar),xbar, sd.post)
  out <- eJTK(timepoints,boot,phase_vec, phase_shift_vec, asymmetry_vec)
  return(out)
}

c(1:10) |> map(\(x) boot_results(xbar, sd.post, timepoints))

# boot_results <- list()
# for (j in 1:100){print(j)
#   boot <- rnorm(length(xbar),xbar, sd.post)
#   boot_results[[j]] <- eJTK(timepoints,boot,phase_vec, phase_shift_vec, asymmetry_vec)
# }



boot_results_df <- do.call(rbind,boot_results)
boot_results_df <- boot_results_df %>% summarize(across(everything(),mean))
boot_results_df$trt <- t
boot_results_df$marker <- markers[i]

all <- bind_rows(boot_results_df, all)
print(all)
  }
}


#apply(boot_results_df[,1:4],2,mean)
