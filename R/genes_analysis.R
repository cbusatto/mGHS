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
p = 381
K = 4

B = 10000
bn = 5000

hyp_ta = get_hyp_ta(0.6, 0.05)

Y1 = scale(expr_all$IFNG)
Y2 = scale(expr_all$LPS2)
Y3 = scale(expr_all$LPS24)
Y4 = scale(expr_all$UNSTIM)

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

S_jointGHS = list()
S_jointGHS[[1]] = Y1
S_jointGHS[[2]] = Y2
S_jointGHS[[3]] = Y3
S_jointGHS[[4]] = Y4


print("\n")
print(":::: jointGHS: ")

res_jointGHS = jointGHS::jointGHS(S_jointGHS, scale = T, AIC_selection = T, AIC_eps = 5,
                                  eps = 1e-3, stop_overflow = TRUE) 

save(n, p, K, res_jointGHS, file = "res_genes.RData")


print("\n")
print(":::: mGHS: ")

set.seed(22)
res_mGHS = mGHS(B, bn, n, S, p, hyp_ta, chain = 1, ping = 1000)

save(n, p, K, res_jointGHS, res_mGHS, file = "res_genes.RData")

print("\n")
print(":::: GHS: ")

res_GHS1 = GHS_est(S[[1]], n[1], bn, B-bn)
res_GHS2 = GHS_est(S[[2]], n[2], bn, B-bn)
res_GHS3 = GHS_est(S[[3]], n[3], bn, B-bn)
res_GHS4 = GHS_est(S[[4]], n[4], bn, B-bn)

save(n, p, K, res_jointGHS, res_mGHS, res_GHS1, res_GHS2, res_GHS3, res_GHS4, file = "res_genes.RData")

print("\n")
print(":::: GBAG: ")

v0_l = c(0.25, 0.5, 0.75, 1) * sqrt(1 / n / log(p))
v1_l = c(2.5, 5, 7.5, 10) * sqrt(1 / n / log(p))
p1 = 0.4

hyper = Tune_GemBag(v0_l, v1_l, S_GBAG, n, maxiter = 20, p1, p_2 = 0.8)
v0 = hyper$v0
v1 = hyper$v1

res_GBAG = GemBag(S_l = S_GBAG, n = n, v_0 = v0, v_1 = v1, tau = v0, p_1 = p1, p_2 = 0.8, maxiter = 20)
names(res_GBAG) = c('Theta', 'P', 'W')

save(n, p, K, res_jointGHS, res_mGHS, res_GBAG, res_GHS1, res_GHS2, res_GHS3, res_GHS4, file = "res_genes.RData")
