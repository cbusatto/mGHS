load("bikeshare_data.RData")

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
K = 6
p = 239

B = 10000
bn = 5000

## set hyperaparameters for cut-model
hyp_ta = get_hyp_ta(0.6, 0.05)


S = list()
S[[1]] = t(X1c) %*% X1c
S[[2]] = t(X2c) %*% X2c
S[[3]] = t(X3c) %*% X3c
S[[4]] = t(X1m) %*% X1m
S[[5]] = t(X2m) %*% X2m
S[[6]] = t(X3m) %*% X3m

S_GBAG = list()
S_GBAG[[1]] = cov(X1c)
S_GBAG[[2]] = cov(X2c)
S_GBAG[[3]] = cov(X3c)
S_GBAG[[4]] = cov(X1m)
S_GBAG[[5]] = cov(X2m)
S_GBAG[[6]] = cov(X3m)

S_jointGHS = list()
S_jointGHS[[1]] = X1c
S_jointGHS[[2]] = X2c
S_jointGHS[[3]] = X3c
S_jointGHS[[4]] = X1m
S_jointGHS[[5]] = X2m
S_jointGHS[[6]] = X3m


print("\n")
print(":::: mGHS: ")

res_mGHS = mGHS_chain(B, bn, n, S, p, hyp_ta, chain = 1, 1500)

save(n, p, K, res_mGHS, file = "res_bikeshare.RData")

print("\n")
print(":::: GHS: ")

res_GHS1 = GHS_est(S[[1]], n[1], bn, B-bn)
res_GHS2 = GHS_est(S[[2]], n[2], bn, B-bn)
res_GHS3 = GHS_est(S[[3]], n[3], bn, B-bn)
res_GHS4 = GHS_est(S[[4]], n[4], bn, B-bn)
res_GHS5 = GHS_est(S[[5]], n[5], bn, B-bn)
res_GHS6 = GHS_est(S[[6]], n[6], bn, B-bn)

save(n, p, K, res_mGHS, res_GHS1, res_GHS2, res_GHS3, res_GHS4, res_GHS5, res_GHS6, file = "res_bikeshare.RData")

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

save(n, p, K, res_mGHS, res_GBAG, res_GHS1, res_GHS2, res_GHS3, res_GHS4, file = "res_bikeshare.RData")

print("\n")
print(":::: jointGHS: ")

res_jointGHS = jointGHS(S_jointGHS, epsilon = 1e-5, AIC_selection = T, AIC_eps = 0.1)

save(n, p, K, res_mGHS, res_GBAG, res_jointGHS, res_GHS1, res_GHS2, res_GHS3, res_GHS4, file = "res_bikeshare.RData")

