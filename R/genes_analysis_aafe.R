load("expression_data_monocytes_LYZ_region_FDR_5.RData")
load("res_genes_pred.RData")

# p = 381
# K = 4

i1 = 1:191
i2 = 192:p

X_test = list(Y1_test, Y2_test, Y3_test, Y4_test)

bn = 5000
n_test = n - n_train

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## mGHS

Omega_mGHS = array(0, c(p, p, K))
Sigma_mGHS = array(0, c(p, p, K))

res_mGHS_aafe = rep(0, K)
for (b in 1:bn) {
  
  for (k in 1:K) {
    diag(Omega_mGHS[, , k]) = res_mGHS$Omega_jj[, (b-1)* K + k]
    
    cs = 1
    for (j in 1:(p-1)) {
      for (js in (j+1):p) {
        Omega_mGHS[js, j, k] = Omega_mGHS[j, js, k] = res_mGHS$Omega_ij[cs, (b-1)* K + k]
        cs = cs + 1
      }
    }
    
    Sigma_mGHS[, , k] = solve(Omega_mGHS[, , k])
  }
  
  AAFE_mGHS = rep(0, K)
  for (k in 1:K) {
    O11 = Sigma_mGHS[i1, i1, k]
    O12 = Sigma_mGHS[i2, i1, k]
    
    for (i in 1:n_test[k]) {
      y_pred = O12 %*% (solve(O11) %*% X_test[[k]][i, i1])
      AAFE_mGHS[k] = AAFE_mGHS[k] + sum(abs(X_test[[k]][i, i2] - y_pred))
      # AAFE_mGHS[k] = AAFE_mGHS[k] + sum((X_test[[k]][i, i2] - y_pred)^2)
    }
    
    AAFE_mGHS[k] = AAFE_mGHS[k] / (190 * n_test[k])
  }
  
  print(b)
  
  res_mGHS_aafe = res_mGHS_aafe + AAFE_mGHS
}

res_mGHS_aafe = res_mGHS_aafe / bn

mean(res_mGHS_aafe)
save(res_mGHS_aafe, file = "res_bikeshare_aafe.RData")

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## GHS

Omega_GHS = array(0, c(p, p, K))
Sigma_GHS = array(0, c(p, p, K))

res_GHS_aafe = rep(0, K)
for (b in 1:bn) {
  
  Omega_GHS[, , 1] = res_GHS1[[1]][, , b]
  Omega_GHS[, , 2] = res_GHS2[[1]][, , b]
  Omega_GHS[, , 3] = res_GHS3[[1]][, , b]
  Omega_GHS[, , 4] = res_GHS4[[1]][, , b]
  
  Sigma_GHS[, , 1] = solve(Omega_GHS[, , 1])
  Sigma_GHS[, , 2] = solve(Omega_GHS[, , 2])
  Sigma_GHS[, , 3] = solve(Omega_GHS[, , 3])
  Sigma_GHS[, , 4] = solve(Omega_GHS[, , 4])
  
  AAFE_GHS = rep(0, K)
  O111 = Sigma_GHS[i1, i1, 1]
  O121 = Sigma_GHS[i2, i1, 1]
  
  O112 = Sigma_GHS[i1, i1, 2]
  O122 = Sigma_GHS[i2, i1, 2]
  
  O113 = Sigma_GHS[i1, i1, 3]
  O123 = Sigma_GHS[i2, i1, 3]
  
  O114 = Sigma_GHS[i1, i1, 4]
  O124 = Sigma_GHS[i2, i1, 4]
  
  for (i in 1:n_test[1]) {
    y_pred = O121 %*% (solve(O111) %*% X_test[[1]][i, i1])
    AAFE_GHS[1] = AAFE_GHS[1] + sum(abs(X_test[[1]][i, i2] - y_pred))
  }
  
  for (i in 1:n_test[2]) {
    y_pred = O122 %*% (solve(O112) %*% X_test[[2]][i, i1])
    AAFE_GHS[2] = AAFE_GHS[2] + sum(abs(X_test[[2]][i, i2] - y_pred))
  }
  
  for (i in 1:n_test[3]) {
    y_pred = O123 %*% (solve(O113) %*% X_test[[3]][i, i1])
    AAFE_GHS[3] = AAFE_GHS[3] + sum(abs(X_test[[3]][i, i2] - y_pred))
  }
  
  for (i in 1:n_test[4]) {
    y_pred = O124 %*% (solve(O114) %*% X_test[[4]][i, i1])
    AAFE_GHS[4] = AAFE_GHS[4] + sum(abs(X_test[[4]][i, i2] - y_pred))
  }
  
  for (k in 1:K) AAFE_GHS[k] = AAFE_GHS[k] / (190 * n_test[k])
  
  print(b)
  res_GHS_aafe = res_GHS_aafe + AAFE_GHS
}
res_GHS_aafe = res_GHS_aafe / bn

mean(res_GHS_aafe)
save(res_mGHS_aafe, res_GHS_aafe, file = "res_bikeshare_aafe.RData")

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## GBAG

Sigma_GBAG = array(0, c(p, p, K))
for (k in 1:K) Sigma_GBAG[, , k] = solve(res_GBAG[[1]][[k]])

AAFE_GBAG = rep(0, K)
for (k in 1:K) {
  O11 = Sigma_GBAG[i1, i1, k]
  O12 = Sigma_GBAG[i2, i1, k]
  
  for (i in 1:n_test[k]) {
    y_pred = O12 %*% (solve(O11) %*% X_test[[k]][i, i1])
    AAFE_GBAG[k] = AAFE_GBAG[k] + sum(abs(X_test[[k]][i, i2] - y_pred))
    # AAFE_GBAG[k] = AAFE_GBAG[k] + sum((X_test[[k]][i, i2] - y_pred)^2)
  }
  
  AAFE_GBAG[k] = AAFE_GBAG[k] / (190 * n_test[k])
}
mean(AAFE_GBAG)
res_GBAG_aafe = AAFE_GBAG

save(res_mGHS_aafe, res_GHS_aafe, res_GBAG_aafe, file = "res_bikeshare_aafe.RData")

## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## jointGHS

# Sigma_jointGHS = array(0, c(p, p, K))
# for (k in 1:K) Sigma_jointGHS[, , k] = solve(res_jointGHS[[2]][[k]])
# # for (k in 1:K) Sigma_jointGHS[, , k] = solve(cov2cor(res_jointGHS[[2]][[k]]))
# 
# AAFE_jointGHS = rep(0, K)
# for (k in 1:K) {
#   O11 = Sigma_jointGHS[i1, i1, k]
#   O12 = Sigma_jointGHS[i2, i1, k]
#   
#   for (i in 1:n_test[k]) {
#     y_pred = O12 %*% (solve(O11) %*% X_test[[k]][i, i1])
#     AAFE_jointGHS[k] = AAFE_jointGHS[k] + sum(abs(X_test[[k]][i, i2] - y_pred))
#     # AAFE_jointGHS[k] = AAFE_jointGHS[k] + sum((X_test[[k]][i, i2] - y_pred)^2)
#   }
#   
#   AAFE_jointGHS[k] = AAFE_jointGHS[k] / (190 * n_test[k])
# }
# mean(AAFE_jointGHS)
# res_jointGHS_aafe = AAFE_jointGHS
# 
# save(res_mGHS_aafe, res_GHS_aafe, res_GBAG_aafe, res_jointGHS_aafe, file = "res_bikeshare_aafe.RData")

