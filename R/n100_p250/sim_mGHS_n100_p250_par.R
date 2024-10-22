require(mvtnorm)
require(pROC)
require(mGHS)
library(StatPerMeCo)
library(foreach)
library(doMC)
registerDoMC(25)

source("utils_simulations.R")
# source("C:/Users/HP/Dropbox/Thesis/mGHS_Project/simulations/utils_simulations.R")

get_hyp_ta = function(mu, var) {
  a = mu^2 * (1 - mu) / var - mu
  b = (a - mu * a) / mu
  
  return(c(a, b))
}

set.seed(10)
K = 4
p = 250
n = rep(100, K)

N_iter = 50

B = 10000
bn = 5000

## (2.17, 1.78)
hyp_ta = get_hyp_ta(0.6, 0.05)

## ::::::::::::::::::::::::::
## coupled

load("c_p250_Coupled.RData")
tm = true_model(list(c1, c2, c3, c4))

RES = matrix(0, N_iter, 7*8 + 1) 

cpar = foreach (is = 1:N_iter, .combine = rbind, .packages = c("Rcpp", "mvtnorm", "pROC", "mGHS", "StatPerMeCo")) %dopar% {
  
  print(":::::::::: iteration: ")
  print(is)
  Y1 = rmvnorm(n[1], rep(0, p), solve(c1))
  Y2 = rmvnorm(n[2], rep(0, p), solve(c2))
  Y3 = rmvnorm(n[3], rep(0, p), solve(c3))
  Y4 = rmvnorm(n[4], rep(0, p), solve(c4))
  
  S = list()
  S[[1]] = t(Y1) %*% Y1
  S[[2]] = t(Y2) %*% Y2
  S[[3]] = t(Y3) %*% Y3
  S[[4]] = t(Y4) %*% Y4
  
  print("\n")
  print(":::: mGHS: ")
  stuck = 1
  ta_reasonable = 1 
  
  while ((stuck == 1) | (ta_reasonable == 1)) {
    res_mGHS = mGHS(B, bn, n, S, p, hyp_ta, chain = 0, ping = 1000)
    
    stuck = res_mGHS$stuck
    
    tmp_mGHS = tmp_mGHS2 = matrix(0, 4, 7)
    if (res_mGHS$stuck == 0) {
      tmp_mGHS = summary_mGHS(res_mGHS, tm, list(c1, c2, c3, c4))
      
      tmp = summary_mGHS3(res_mGHS, tm, list(c1, c2, c3, c4), res_mGHS$ta)
      tmp_mGHS2 = tmp$res
      ta_reasonable = tmp$ta_reasonable
    }
  }
  
  print("\n")
  print("stuck iter: ") 
  print(is)
  print(": ")
  print(res_mGHS$stuck)
  
  RES[is, ] = c(tmp_mGHS[1, ], tmp_mGHS[2, ], tmp_mGHS[3, ], tmp_mGHS[4, ],
                tmp_mGHS2[1, ], tmp_mGHS2[2, ], tmp_mGHS2[3, ], tmp_mGHS2[4, ], res_mGHS$stuck)
}

res_mean_final = apply(cpar, 2, mean)
res_sd_final = apply(cpar, 2, sd)

res_g1_mean = matrix(c(res_mean_final[c(1:7, 29:35)]), 2, 7, byrow = T)
colnames(res_g1_mean) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g1_mean) = c("mGHS", "mGHS2")

res_g2_mean = matrix(c(res_mean_final[c(8:14, 36:42)]), 2, 7, byrow = T)
colnames(res_g2_mean) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g2_mean) = c("mGHS", "mGHS2")

res_g3_mean = matrix(c(res_mean_final[c(15:21, 43:49)]), 2, 7, byrow = T)
colnames(res_g3_mean) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g3_mean) = c("mGHS", "mGHS2")

res_g4_mean = matrix(c(res_mean_final[c(22:28, 50:56)]), 2, 7, byrow = T)
colnames(res_g4_mean) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g4_mean) = c("mGHS", "mGHS2")

res_g1_sd = matrix(c(res_sd_final[c(1:7, 29:35)]), 2, 7, byrow = T)
colnames(res_g1_sd) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g1_sd) = c("mGHS", "mGHS2")

res_g2_sd = matrix(c(res_sd_final[c(8:14, 36:42)]), 2, 7, byrow = T)
colnames(res_g2_sd) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g2_sd) = c("mGHS", "mGHS2")

res_g3_sd = matrix(c(res_sd_final[c(15:21, 43:49)]), 2, 7, byrow = T)
colnames(res_g3_sd) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g3_sd) = c("mGHS", "mGHS2")

res_g4_sd = matrix(c(res_sd_final[c(22:28, 50:56)]), 2, 7, byrow = T)
colnames(res_g4_sd) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g4_sd) = c("mGHS", "mGHS2")

res_Coupled_g1_mean = res_g1_mean
res_Coupled_g2_mean = res_g2_mean
res_Coupled_g3_mean = res_g3_mean
res_Coupled_g4_mean = res_g4_mean

res_Coupled_g1_sd = res_g1_sd
res_Coupled_g2_sd = res_g2_sd
res_Coupled_g3_sd = res_g3_sd
res_Coupled_g4_sd = res_g4_sd

res_Coupled_mean = res_Coupled_g1_mean + res_Coupled_g2_mean + res_Coupled_g3_mean + res_Coupled_g4_mean
res_Coupled_sd = res_Coupled_g1_sd + res_Coupled_g2_sd + res_Coupled_g3_sd + res_Coupled_g4_sd

stuck_index = sum(cpar[, 57])

save(cpar, stuck_index, res_Coupled_mean, res_Coupled_g1_mean, res_Coupled_g2_mean, res_Coupled_g3_mean, res_Coupled_g4_mean, file = "sim_mGHS_n100_p250_mean_Coupled_last.RData")
save(res_Coupled_sd, res_Coupled_g1_sd, res_Coupled_g2_sd, res_Coupled_g3_sd, res_Coupled_g4_sd, file = "sim_mGHS_n100_p250_sd_Coupled_last.RData")


## ::::::::::::::::::::::::::
## Independent

load("c_p250_Independent.RData")
tm = true_model(list(c1, c2, c3, c4))

RES = matrix(0, N_iter, 7*8 + 1) 

cpar = foreach (is = 1:N_iter, .combine = rbind, .packages = c("Rcpp", "mvtnorm", "pROC", "mGHS", "StatPerMeCo")) %dopar% {
  
  print(":::::::::: iteration: ")
  print(is)
  Y1 = rmvnorm(n[1], rep(0, p), solve(c1))
  Y2 = rmvnorm(n[2], rep(0, p), solve(c2))
  Y3 = rmvnorm(n[3], rep(0, p), solve(c3))
  Y4 = rmvnorm(n[4], rep(0, p), solve(c4))
  
  S = list()
  S[[1]] = t(Y1) %*% Y1
  S[[2]] = t(Y2) %*% Y2
  S[[3]] = t(Y3) %*% Y3
  S[[4]] = t(Y4) %*% Y4
  
  print("\n")
  print(":::: mGHS: ")
  stuck = 1
  ta_reasonable = 1 
  
  while ((stuck == 1) | (ta_reasonable == 1)) {
    res_mGHS = mGHS(B, bn, n, S, p, hyp_ta, chain = 0, ping = 1000)
    
    stuck = res_mGHS$stuck
    
    tmp_mGHS = tmp_mGHS2 = matrix(0, 4, 7)
    if (res_mGHS$stuck == 0) {
      tmp_mGHS = summary_mGHS(res_mGHS, tm, list(c1, c2, c3, c4))
      
      tmp = summary_mGHS3(res_mGHS, tm, list(c1, c2, c3, c4), res_mGHS$ta)
      tmp_mGHS2 = tmp$res
      ta_reasonable = tmp$ta_reasonable
    }
  }
  
  print("\n")
  print("stuck iter: ") 
  print(is)
  print(": ")
  print(res_mGHS$stuck)
  
  RES[is, ] = c(tmp_mGHS[1, ], tmp_mGHS[2, ], tmp_mGHS[3, ], tmp_mGHS[4, ],
                tmp_mGHS2[1, ], tmp_mGHS2[2, ], tmp_mGHS2[3, ], tmp_mGHS2[4, ], res_mGHS$stuck)
}

res_mean_final = apply(cpar, 2, mean)
res_sd_final = apply(cpar, 2, sd)

res_g1_mean = matrix(c(res_mean_final[c(1:7, 29:35)]), 2, 7, byrow = T)
colnames(res_g1_mean) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g1_mean) = c("mGHS", "mGHS2")

res_g2_mean = matrix(c(res_mean_final[c(8:14, 36:42)]), 2, 7, byrow = T)
colnames(res_g2_mean) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g2_mean) = c("mGHS", "mGHS2")

res_g3_mean = matrix(c(res_mean_final[c(15:21, 43:49)]), 2, 7, byrow = T)
colnames(res_g3_mean) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g3_mean) = c("mGHS", "mGHS2")

res_g4_mean = matrix(c(res_mean_final[c(22:28, 50:56)]), 2, 7, byrow = T)
colnames(res_g4_mean) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g4_mean) = c("mGHS", "mGHS2")

res_g1_sd = matrix(c(res_sd_final[c(1:7, 29:35)]), 2, 7, byrow = T)
colnames(res_g1_sd) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g1_sd) = c("mGHS", "mGHS2")

res_g2_sd = matrix(c(res_sd_final[c(8:14, 36:42)]), 2, 7, byrow = T)
colnames(res_g2_sd) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g2_sd) = c("mGHS", "mGHS2")

res_g3_sd = matrix(c(res_sd_final[c(15:21, 43:49)]), 2, 7, byrow = T)
colnames(res_g3_sd) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g3_sd) = c("mGHS", "mGHS2")

res_g4_sd = matrix(c(res_sd_final[c(22:28, 50:56)]), 2, 7, byrow = T)
colnames(res_g4_sd) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g4_sd) = c("mGHS", "mGHS2")

res_Independent_g1_mean = res_g1_mean
res_Independent_g2_mean = res_g2_mean
res_Independent_g3_mean = res_g3_mean
res_Independent_g4_mean = res_g4_mean

res_Independent_g1_sd = res_g1_sd
res_Independent_g2_sd = res_g2_sd
res_Independent_g3_sd = res_g3_sd
res_Independent_g4_sd = res_g4_sd

res_Independent_mean = res_Independent_g1_mean + res_Independent_g2_mean + res_Independent_g3_mean + res_Independent_g4_mean
res_Independent_sd = res_Independent_g1_sd + res_Independent_g2_sd + res_Independent_g3_sd + res_Independent_g4_sd

stuck_index = sum(cpar[, 57])

save(cpar, stuck_index, res_Independent_mean, res_Independent_g1_mean, res_Independent_g2_mean, res_Independent_g3_mean, res_Independent_g4_mean, file = "sim_mGHS_n100_p250_mean_Independent_last.RData")
save(res_Independent_sd, res_Independent_g1_sd, res_Independent_g2_sd, res_Independent_g3_sd, res_Independent_g4_sd, file = "sim_mGHS_n100_p250_sd_Independent_last.RData")


## ::::::::::::::::::::::::::
## Full

load("c_p250_Full.RData")
tm = true_model(list(c1, c2, c3, c4))

RES = matrix(0, N_iter, 7*8 + 1) 

cpar = foreach (is = 1:N_iter, .combine = rbind, .packages = c("Rcpp", "mvtnorm", "pROC", "mGHS", "StatPerMeCo")) %dopar% {
  
  print(":::::::::: iteration: ")
  print(is)
  Y1 = rmvnorm(n[1], rep(0, p), solve(c1))
  Y2 = rmvnorm(n[2], rep(0, p), solve(c2))
  Y3 = rmvnorm(n[3], rep(0, p), solve(c3))
  Y4 = rmvnorm(n[4], rep(0, p), solve(c4))
  
  S = list()
  S[[1]] = t(Y1) %*% Y1
  S[[2]] = t(Y2) %*% Y2
  S[[3]] = t(Y3) %*% Y3
  S[[4]] = t(Y4) %*% Y4
  
  print("\n")
  print(":::: mGHS: ")
  stuck = 1
  ta_reasonable = 1 
  
  while ((stuck == 1) | (ta_reasonable == 1)) {
    res_mGHS = mGHS(B, bn, n, S, p, hyp_ta, chain = 0, ping = 1000)
    
    stuck = res_mGHS$stuck
    
    tmp_mGHS = tmp_mGHS2 = matrix(0, 4, 7)
    if (res_mGHS$stuck == 0) {
      tmp_mGHS = summary_mGHS(res_mGHS, tm, list(c1, c2, c3, c4))
      
      tmp = summary_mGHS3(res_mGHS, tm, list(c1, c2, c3, c4), res_mGHS$ta)
      tmp_mGHS2 = tmp$res
      ta_reasonable = tmp$ta_reasonable
    }
  }
  
  print("\n")
  print("stuck iter: ") 
  print(is)
  print(": ")
  print(res_mGHS$stuck)
  
  RES[is, ] = c(tmp_mGHS[1, ], tmp_mGHS[2, ], tmp_mGHS[3, ], tmp_mGHS[4, ],
                tmp_mGHS2[1, ], tmp_mGHS2[2, ], tmp_mGHS2[3, ], tmp_mGHS2[4, ], res_mGHS$stuck)
}

res_mean_final = apply(cpar, 2, mean)
res_sd_final = apply(cpar, 2, sd)

res_g1_mean = matrix(c(res_mean_final[c(1:7, 29:35)]), 2, 7, byrow = T)
colnames(res_g1_mean) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g1_mean) = c("mGHS", "mGHS2")

res_g2_mean = matrix(c(res_mean_final[c(8:14, 36:42)]), 2, 7, byrow = T)
colnames(res_g2_mean) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g2_mean) = c("mGHS", "mGHS2")

res_g3_mean = matrix(c(res_mean_final[c(15:21, 43:49)]), 2, 7, byrow = T)
colnames(res_g3_mean) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g3_mean) = c("mGHS", "mGHS2")

res_g4_mean = matrix(c(res_mean_final[c(22:28, 50:56)]), 2, 7, byrow = T)
colnames(res_g4_mean) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g4_mean) = c("mGHS", "mGHS2")

res_g1_sd = matrix(c(res_sd_final[c(1:7, 29:35)]), 2, 7, byrow = T)
colnames(res_g1_sd) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g1_sd) = c("mGHS", "mGHS2")

res_g2_sd = matrix(c(res_sd_final[c(8:14, 36:42)]), 2, 7, byrow = T)
colnames(res_g2_sd) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g2_sd) = c("mGHS", "mGHS2")

res_g3_sd = matrix(c(res_sd_final[c(15:21, 43:49)]), 2, 7, byrow = T)
colnames(res_g3_sd) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g3_sd) = c("mGHS", "mGHS2")

res_g4_sd = matrix(c(res_sd_final[c(22:28, 50:56)]), 2, 7, byrow = T)
colnames(res_g4_sd) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g4_sd) = c("mGHS", "mGHS2")

res_Full_g1_mean = res_g1_mean
res_Full_g2_mean = res_g2_mean
res_Full_g3_mean = res_g3_mean
res_Full_g4_mean = res_g4_mean

res_Full_g1_sd = res_g1_sd
res_Full_g2_sd = res_g2_sd
res_Full_g3_sd = res_g3_sd
res_Full_g4_sd = res_g4_sd

res_Full_mean = res_Full_g1_mean + res_Full_g2_mean + res_Full_g3_mean + res_Full_g4_mean
res_Full_sd = res_Full_g1_sd + res_Full_g2_sd + res_Full_g3_sd + res_Full_g4_sd

stuck_index = sum(cpar[, 57])

save(cpar, stuck_index, res_Full_mean, res_Full_g1_mean, res_Full_g2_mean, res_Full_g3_mean, res_Full_g4_mean, file = "sim_mGHS_n100_p250_mean_Full_last.RData")
save(res_Full_sd, res_Full_g1_sd, res_Full_g2_sd, res_Full_g3_sd, res_Full_g4_sd, file = "sim_mGHS_n100_p250_sd_Full_last.RData")


## ::::::::::::::::::::::::::
## P2020

load("c_p250_P2020.RData")
tm = true_model(list(c1, c2, c3, c4))

RES = matrix(0, N_iter, 7*8 + 1) 

cpar = foreach (is = 1:N_iter, .combine = rbind, .packages = c("Rcpp", "mvtnorm", "pROC", "mGHS", "StatPerMeCo")) %dopar% {
  
  print(":::::::::: iteration: ")
  print(is)
  Y1 = rmvnorm(n[1], rep(0, p), solve(c1))
  Y2 = rmvnorm(n[2], rep(0, p), solve(c2))
  Y3 = rmvnorm(n[3], rep(0, p), solve(c3))
  Y4 = rmvnorm(n[4], rep(0, p), solve(c4))
  
  S = list()
  S[[1]] = t(Y1) %*% Y1
  S[[2]] = t(Y2) %*% Y2
  S[[3]] = t(Y3) %*% Y3
  S[[4]] = t(Y4) %*% Y4
  
  print("\n")
  print(":::: mGHS: ")
  stuck = 1
  ta_reasonable = 1 
  
  while ((stuck == 1) | (ta_reasonable == 1)) {
    res_mGHS = mGHS(B, bn, n, S, p, hyp_ta, chain = 0, ping = 1000)
    
    stuck = res_mGHS$stuck
    
    tmp_mGHS = tmp_mGHS2 = matrix(0, 4, 7)
    if (res_mGHS$stuck == 0) {
      tmp_mGHS = summary_mGHS(res_mGHS, tm, list(c1, c2, c3, c4))
      
      tmp = summary_mGHS3(res_mGHS, tm, list(c1, c2, c3, c4), res_mGHS$ta)
      tmp_mGHS2 = tmp$res
      ta_reasonable = tmp$ta_reasonable
    }
  }  
  
  print("\n")
  print("stuck iter: ") 
  print(is)
  print(": ")
  print(res_mGHS$stuck)
  
  RES[is, ] = c(tmp_mGHS[1, ], tmp_mGHS[2, ], tmp_mGHS[3, ], tmp_mGHS[4, ],
                tmp_mGHS2[1, ], tmp_mGHS2[2, ], tmp_mGHS2[3, ], tmp_mGHS2[4, ], res_mGHS$stuck)
}

res_mean_final = apply(cpar, 2, mean)
res_sd_final = apply(cpar, 2, sd)

res_g1_mean = matrix(c(res_mean_final[c(1:7, 29:35)]), 2, 7, byrow = T)
colnames(res_g1_mean) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g1_mean) = c("mGHS", "mGHS2")

res_g2_mean = matrix(c(res_mean_final[c(8:14, 36:42)]), 2, 7, byrow = T)
colnames(res_g2_mean) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g2_mean) = c("mGHS", "mGHS2")

res_g3_mean = matrix(c(res_mean_final[c(15:21, 43:49)]), 2, 7, byrow = T)
colnames(res_g3_mean) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g3_mean) = c("mGHS", "mGHS2")

res_g4_mean = matrix(c(res_mean_final[c(22:28, 50:56)]), 2, 7, byrow = T)
colnames(res_g4_mean) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g4_mean) = c("mGHS", "mGHS2")

res_g1_sd = matrix(c(res_sd_final[c(1:7, 29:35)]), 2, 7, byrow = T)
colnames(res_g1_sd) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g1_sd) = c("mGHS", "mGHS2")

res_g2_sd = matrix(c(res_sd_final[c(8:14, 36:42)]), 2, 7, byrow = T)
colnames(res_g2_sd) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g2_sd) = c("mGHS", "mGHS2")

res_g3_sd = matrix(c(res_sd_final[c(15:21, 43:49)]), 2, 7, byrow = T)
colnames(res_g3_sd) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g3_sd) = c("mGHS", "mGHS2")

res_g4_sd = matrix(c(res_sd_final[c(22:28, 50:56)]), 2, 7, byrow = T)
colnames(res_g4_sd) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g4_sd) = c("mGHS", "mGHS2")

res_P2020_g1_mean = res_g1_mean
res_P2020_g2_mean = res_g2_mean
res_P2020_g3_mean = res_g3_mean
res_P2020_g4_mean = res_g4_mean

res_P2020_g1_sd = res_g1_sd
res_P2020_g2_sd = res_g2_sd
res_P2020_g3_sd = res_g3_sd
res_P2020_g4_sd = res_g4_sd

res_P2020_mean = res_P2020_g1_mean + res_P2020_g2_mean + res_P2020_g3_mean + res_P2020_g4_mean
res_P2020_sd = res_P2020_g1_sd + res_P2020_g2_sd + res_P2020_g3_sd + res_P2020_g4_sd

stuck_index = sum(cpar[, 57])

save(cpar, stuck_index, res_P2020_mean, res_P2020_g1_mean, res_P2020_g2_mean, res_P2020_g3_mean, res_P2020_g4_mean, file = "sim_mGHS_n100_p250_mean_P2020_last.RData")
save(res_P2020_sd, res_P2020_g1_sd, res_P2020_g2_sd, res_P2020_g3_sd, res_P2020_g4_sd, file = "sim_mGHS_n100_p250_sd_P2020_last.RData")


