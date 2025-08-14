source("./jtk1.R")
source("./jtk_n.R")
source("./eJTK.R")
library(tidyverse)
library(limma)
library(DiscoRhythm)
library(vashr)
library(furrr)

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
phase_shift_vec <- seq(0,22,2)
asymmetry_vec <- seq(2,22,2)

timepoints <- seq(0,46,2)

nsim <- 1000
start <- Sys.time()
#Create the null distribution.  
#Null is a function of phase_vec, phase_shift_vec, asymmetry_vec, and timepoints
null <- map_dfr(1:nsim, ~ {
  jtk_n(
    timepoints,
    sample(1:length(timepoints), size = length(timepoints), replace = FALSE),
    phase_vec,
    phase_shift_vec,
    asymmetry_vec
  )
})
end <- Sys.time()
end - start

library(MASS)
#Estimate parameters
null$tau[null$tau <=0] <- 0.00001
fit <- fitdistr(null$tau, "gamma")

#All the markers to test
#markers <- test %>% select(zymostenol:`134`) %>% names()
#test_long <- test |> pivot_longer(cols = zymostenol:`134`, names_to = "markers", values_to = "values")
#Works!   
#test_long %>% split(list(test_long$trt,test_long$markers)) %>% map(\(x) lm(values ~ time, data = x))
##########################################
# Do eJTK on bootstrap data sets
##########################################
wrapperfun <- function(sub){
  print(sub$trt[1])
  print(sub$markers[1])
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
  #sehat <- sqrt(1/voom$weights[,1])/sqrt(reps) #We had this
  sehat <- sqrt(1/voom$weights[,1]) #They don't divide by reps here?
  vashout <- vash(sehat,reps)
  sd.post <- vashout$sd.post
  #cbind(sd.post, sehat)
  #Now bootstrap
  boot_results <- function(xbar, sd.post, timepoints){
    boot <- rnorm(length(xbar),xbar, sd.post)
    out <- eJTK(timepoints,boot,phase_vec, phase_shift_vec, asymmetry_vec, shape = fit$estimate["shape"], rate = fit$estimate["rate"])
    return(out)
  }
  
  nboot <- 25
  results <- c(1:nboot) |> map(\(x) boot_results(xbar, sd.post, timepoints), .progress = TRUE) |> bind_rows() |> mutate(trt = sub$trt[1], marker = sub$markers[1], id = 1:nboot)
  
  return(results)
}


all <- test_long %>% group_split(trt,markers) %>% map(\(x) wrapperfun(x))
save(all, file = "/Users/gregorymatthews/Dropbox/RbootEJTK_git/all.RData")

all_stack <- bind_rows(all)
summ <- all_stack %>% mutate(arctanhtau = atanh(ifelse(tau > 0.99, .99, tau))) %>% group_by(trt, marker) %>%  summarize(
  mn_phase = mean(phase),
  mn_shift = mean(shift),
  mn_asymmetry = mean(asymmetry),
  mn_arctanhtau = mean(arctanhtau)
) %>% mutate(pvalue = pgamma(tanh(mn_arctanhtau), shape = fit$estimate["shape"],rate = fit$estimate["rate"], lower.tail = FALSE)) %>% arrange(pvalue) %>% view()



sub <-test_long %>% filter(trt == t, markers == markers[1]) 



    # boot_results <- list()
    # for (j in 1:100){print(j)
    #   boot <- rnorm(length(xbar),xbar, sd.post)
    #   boot_results[[j]] <- eJTK(timepoints,boot,phase_vec, phase_shift_vec, asymmetry_vec)
    # }
    
    
    
    # boot_results_df <- do.call(rbind,boot_results)
    # boot_results_df <- boot_results_df %>% summarize(across(everything(),mean))
    # boot_results_df$trt <- t
    # boot_results_df$marker <- markers[i]
    
    #all <- bind_rows(results, all)
    #print(all)


#apply(boot_results_df[,1:4],2,mean)
