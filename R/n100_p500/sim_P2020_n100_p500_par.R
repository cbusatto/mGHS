require(mvtnorm)
require(pROC)
require(mGHS)
library(StatPerMeCo)
library(foreach)
library(doMC)
registerDoMC(10)

source("utils_simulations.R")
source('BIC_GemBag.R')
Rcpp::sourceCpp('GemBag-algo.cpp')

load("c_p500_P2020.RData")

get_hyp_ta = function(mu, var) {
  a = mu^2 * (1 - mu) / var - mu
  b = (a - mu * a) / mu
  
  return(c(a, b))
}

set.seed(10)
K = 4
p = 500

N_iter = 10

B = 10000
bn = 5000

tm = true_model(list(c1, c2, c3, c4))

## (2.17, 1.78)
hyp_ta = get_hyp_ta(0.6, 0.05)

## ::::::::::::::::::::::::::
## n = 50

n = rep(100, K)

RES = matrix(0, N_iter, 7*12 + 1)

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
  
  S_GBAG = list()
  S_GBAG[[1]] = cov(Y1)
  S_GBAG[[2]] = cov(Y2)
  S_GBAG[[3]] = cov(Y3)
  S_GBAG[[4]] = cov(Y4)
  
  ## ::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ## multiple Graphical Horseshoe
  print("\n")
  print(":::: mGHS: ")
  
  res_mGHS = mGHS(B, bn, n, S, p, hyp_ta, 1500)
  
  print("\n")
  print("stuck iter: ") 
  print(is)
  print(": ")
  print(res_mGHS$stuck)
  
  # ------------------------ Joint Estimation of Multiple Graphical Models ------------------------ 
  
  print("\n")
  print(":::: GBAG: ")
  
  # Set of hyperparameter
  v0_l = c(0.25, 0.5, 0.75, 1) * sqrt(1 / n / log(p))
  v1_l = c(2.5, 5, 7.5, 10) * sqrt(1 / n / log(p))
  p1 = 0.4
  
  # Tuning by BIC
  hyper = Tune_GemBag(v0_l, v1_l, S_GBAG, n, maxiter = 20, p1, p_2 = 0.8)
  v0 = hyper$v0
  v1 = hyper$v1
  
  # Output:
  #   Theta: estimated precision matrices
  #   P: estimated posterior inclusion probabilities
  #   W: estimated covariance matrices
  res_GBAG = GemBag(S_l = S_GBAG, n = n, v_0 = v0, v_1 = v1, tau = v0, p_1 = p1, p_2 = 0.8, maxiter = 20)
  names(res_GBAG) = c('Theta', 'P', 'W')
  
  ## ::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ## Performance indexes
  
  ## mGHS
  
  tmp_mGHS = tmp_mGHS2 = matrix(0, 4, 7)
  if (res_mGHS$stuck == 0) {
    tmp_mGHS = summary_mGHS(res_mGHS, tm, list(c1, c2, c3, c4))
    tmp_mGHS2 = summary_mGHS2(res_mGHS, tm, list(c1, c2, c3, c4), res_mGHS$ta)
  }
  
  ## GBAG
  tmp_GBAG = summary_GBAG(res_GBAG, tm, list(c1, c2, c3, c4))
  
  RES[is, ] = c(tmp_mGHS[1, ], tmp_mGHS[2, ], tmp_mGHS[3, ], tmp_mGHS[4, ],
                tmp_mGHS2[1, ], tmp_mGHS2[2, ], tmp_mGHS2[3, ], tmp_mGHS2[4, ],
                tmp_GBAG[1, ], tmp_GBAG[2, ], tmp_GBAG[3, ], tmp_GBAG[4, ], res_mGHS$stuck)
}

res_mean_final = apply(cpar, 2, mean)
res_sd_final = apply(cpar, 2, sd)

res_g1_mean = matrix(c(res_mean_final[c(1:7, 29:35, 57:63)]), 3, 7, byrow = T)
colnames(res_g1_mean) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g1_mean) = c("mGHS", "mGHS2", "GBAG")

res_g2_mean = matrix(c(res_mean_final[c(8:14, 36:42, 64:70)]), 3, 7, byrow = T)
colnames(res_g2_mean) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g2_mean) = c("mGHS", "mGHS2", "GBAG")

res_g3_mean = matrix(c(res_mean_final[c(15:21, 43:49, 71:77)]), 3, 7, byrow = T)
colnames(res_g3_mean) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g3_mean) = c("mGHS", "mGHS2", "GBAG")

res_g4_mean = matrix(c(res_mean_final[c(22:28, 50:56, 78:84)]), 3, 7, byrow = T)
colnames(res_g4_mean) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g4_mean) = c("mGHS", "mGHS2", "GBAG")

res_g1_sd = matrix(c(res_sd_final[c(1:7, 29:35, 57:63)]), 3, 7, byrow = T)
colnames(res_g1_sd) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g1_sd) = c("mGHS", "mGHS2", "GBAG")

res_g2_sd = matrix(c(res_sd_final[c(8:14, 36:42, 64:70)]), 3, 7, byrow = T)
colnames(res_g2_sd) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g2_sd) = c("mGHS", "mGHS2", "GBAG")

res_g3_sd = matrix(c(res_sd_final[c(15:21, 43:49, 71:77)]), 3, 7, byrow = T)
colnames(res_g3_sd) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g3_sd) = c("mGHS", "mGHS2", "GBAG")

res_g4_sd = matrix(c(res_sd_final[c(22:28, 50:56, 78:84)]), 3, 7, byrow = T)
colnames(res_g4_sd) = c("Acc", "MCC", "TPR", "FPR", "AUC", "Fr Loss", "Fr2 Loss")
rownames(res_g4_sd) = c("mGHS", "mGHS2", "GBAG")

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

stuck_index = sum(cpar[, 85])

save(stuck_index, res_P2020_mean, res_P2020_g1_mean, res_P2020_g2_mean, res_P2020_g3_mean, res_P2020_g4_mean, file = "sim_n100_p500_mean_P2020.RData")
save(res_P2020_sd, res_P2020_g1_sd, res_P2020_g2_sd, res_P2020_g3_sd, res_P2020_g4_sd, file = "sim_n100_p500_sd_P2020.RData")

