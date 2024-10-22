load("expression_data_monocytes_LYZ_region_FDR_5.RData")

require(mGHS)
require(GHS)
require(jointGHS)
source('BIC_GemBag.R')
Rcpp::sourceCpp('GemBag-algo.cpp')

get_hyp_ta = function(mu, var) {
  a = mu^2 * (1 - mu) / var - mu
  b = (a - mu * a) / mu
  
  return(c(a, b))
}

set.seed(10)

n = c(366, 260, 321, 413)
n_train = round(0.8 * n)

ind_train1 = sample(1:n[1], n_train[1])
ind_train2 = sample(1:n[2], n_train[2])
ind_train3 = sample(1:n[3], n_train[3])
ind_train4 = sample(1:n[4], n_train[4])

ind_test1 = setdiff(1:n[1], ind_train1)
ind_test2 = setdiff(1:n[2], ind_train2)
ind_test3 = setdiff(1:n[3], ind_train3)
ind_test4 = setdiff(1:n[4], ind_train4)

p = 381
K = 4

B = 10000
bn = 5000

hyp_ta = get_hyp_ta(0.6, 0.05)

Y1 = scale(expr_all$IFNG)
Y2 = scale(expr_all$LPS2)
Y3 = scale(expr_all$LPS24)
Y4 = scale(expr_all$UNSTIM)

Y1_train = Y1[ind_train1, ]
Y2_train = Y2[ind_train2, ]
Y3_train = Y3[ind_train3, ]
Y4_train = Y4[ind_train4, ]

Y1_test = Y1[ind_test1, ]
Y2_test = Y2[ind_test2, ]
Y3_test = Y3[ind_test3, ]
Y4_test = Y4[ind_test4, ]

S_train = list()
S_train[[1]] = t(Y1_train) %*% Y1_train
S_train[[2]] = t(Y2_train) %*% Y2_train
S_train[[3]] = t(Y3_train) %*% Y3_train
S_train[[4]] = t(Y4_train) %*% Y4_train

S_GBAG_train = list()
S_GBAG_train[[1]] = cov(Y1_train)
S_GBAG_train[[2]] = cov(Y2_train)
S_GBAG_train[[3]] = cov(Y3_train)
S_GBAG_train[[4]] = cov(Y4_train)

S_jointGHS_train = list()
S_jointGHS_train[[1]] = Y1_train
S_jointGHS_train[[2]] = Y2_train
S_jointGHS_train[[3]] = Y3_train
S_jointGHS_train[[4]] = Y4_train

print("\n")
print(":::: jointGHS: ")

res_jointGHS = jointGHS::jointGHS(S_jointGHS_train, scale = T, AIC_selection = T, AIC_eps = 5,
                                  eps = 1e-3, stop_overflow = TRUE) 

save(n, n_train, p, K, res_jointGHS,
     Y1_train, Y2_train, Y3_train, Y4_train, 
     Y1_test, Y2_test, Y3_test, Y4_test, file = "res_genes_pred.RData")

print("\n")
print(":::: mGHS: ")

res_mGHS = mGHS(B, bn, n_train, S_train, p, hyp_ta, chain = 1, ping = 1000)

save(n, n_train, p, K, res_mGHS, res_jointGHS,
     Y1_train, Y2_train, Y3_train, Y4_train, 
     Y1_test, Y2_test, Y3_test, Y4_test, file = "res_genes_pred.RData")

print("\n")
print(":::: GHS: ")

res_GHS1 = GHS_est(S_train[[1]], n_train[1], bn, B-bn)
res_GHS2 = GHS_est(S_train[[2]], n_train[2], bn, B-bn)
res_GHS3 = GHS_est(S_train[[3]], n_train[3], bn, B-bn)
res_GHS4 = GHS_est(S_train[[4]], n_train[4], bn, B-bn)

save(n, n_train, p, K, res_mGHS, res_jointGHS, res_GHS1, res_GHS2, res_GHS3, res_GHS4,
     Y1_train, Y2_train, Y3_train, Y4_train, 
     Y1_test, Y2_test, Y3_test, Y4_test, file = "res_genes_pred.RData")

print("\n")
print(":::: GBAG: ")

v0_l = c(0.25, 0.5, 0.75, 1) * sqrt(1 / n_train / log(p))
v1_l = c(2.5, 5, 7.5, 10) * sqrt(1 / n_train / log(p))
p1 = 0.4

hyper = Tune_GemBag(v0_l, v1_l, S_GBAG_train, n_train, maxiter = 20, p1, p_2 = 0.8)
v0 = hyper$v0
v1 = hyper$v1

res_GBAG = GemBag(S_l = S_GBAG_train, n = n_train, v_0 = v0, v_1 = v1, tau = v0, p_1 = p1, p_2 = 0.8, maxiter = 20)
names(res_GBAG) = c('Theta', 'P', 'W')

save(n, n_train, p, K, res_mGHS, res_jointGHS, res_GBAG, res_GHS1, res_GHS2, res_GHS3, res_GHS4,
     Y1_train, Y2_train, Y3_train, Y4_train, 
     Y1_test, Y2_test, Y3_test, Y4_test, file = "res_genes_pred.RData")

