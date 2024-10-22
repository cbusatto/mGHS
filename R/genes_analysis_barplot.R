###########################################################################################################
## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## mGHS
## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
###########################################################################################################

load("res_genes.RData")


mean_lambda = matrix(0, p*(p-1)/2, K)
for (b in 1:5000) mean_lambda = mean_lambda + 1 / (1+res_mGHS$lambda[, ((b-1)*K+1):(b*K)])
mean_lambda = mean_lambda / 5000

t = 0.8
sum(1 - mean_lambda[, 1] > t)
sum(1 - mean_lambda[, 2] > t)
sum(1 - mean_lambda[, 3] > t)
sum(1 - mean_lambda[, 4] > t)

ta_est = apply(res_mGHS$ta[, c(1:5000)], 1, mean)
sum(1 - mean_lambda[, 1] > ta_est[1])
sum(1 - mean_lambda[, 2] > ta_est[2])
sum(1 - mean_lambda[, 3] > ta_est[3])
sum(1 - mean_lambda[, 4] > ta_est[4])

adj_mGHS = array(0, c(p, p, K))
cs = 1
for (k in 1:K) {
  cs = 1
  for (j in 1:(p-1)) {
    for (js in (j+1):p) {
      if (1-mean_lambda[cs, k] > t) {
        # if (1-mean_lambda[cs, k] > ta_est[k]) {
        adj_mGHS[js, j, k] = adj_mGHS[j, js, k] = 1
      }
      
      cs = cs + 1
    }
  }
}

mat_mGHS = matrix(0, 4, 381)

### IFN
tmp = apply(adj_mGHS[, , 1], 1, sum)
mat_mGHS[1, ] = tmp[order(tmp)]
id_IFN_mGHS = genes_id[order(tmp)]

### LPS2
tmp = apply(adj_mGHS[, , 2], 1, sum)
mat_mGHS[2, ] = tmp[order(tmp)]
id_LPS2_mGHS = genes_id[order(tmp)]

### LPS24
tmp = apply(adj_mGHS[, , 3], 1, sum)
mat_mGHS[3, ] = tmp[order(tmp)]
id_LPS24_mGHS = genes_id[order(tmp)]

### UNSTIM
tmp = apply(adj_mGHS[, , 4], 1, sum)
mat_mGHS[4, ] = tmp[order(tmp)]
id_UNSTIM_mGHS = genes_id[order(tmp)]


## select top-30 hub-genes
id_res_mGHS = matrix(c(id_IFN_mGHS, id_LPS2_mGHS, id_LPS24_mGHS, id_UNSTIM_mGHS), nrow = 4, byrow = T)
id_res_mGHS[, 352:381]
mat_mGHS[, 352:381]

ID = "LGALS3"
ID = "COX6A1"
ID = "CYP27A1"
ID = "LYZ"
ID = "YEATS4"
ID = "AFMID"
ID = "AIRE"

par(mfrow = c(2, 2))
plot(1:381, mat_mGHS[1, ], type = "h", col = "black", xaxt='n', xlab = "", ylab = "number of connections", main = expression(paste("IFN-", gamma, sep = "")))

segments(which(id_IFN_mGHS == "LYZ"), y0 = 0, y1 = 31, lty = 2)
axis(1, at=which(id_IFN_mGHS == "LYZ"), lwd.ticks = -1, labels="LYZ", las = 2)

segments(which(id_IFN_mGHS == "LGALS3"), y0 = 0, y1 = 31, lty = 2)
axis(1, at=which(id_IFN_mGHS == "LGALS3"), lwd.ticks = -1, labels="LGALS3", las = 2)

segments(which(id_IFN_mGHS == "CYP27A1"), y0 = 0, y1 = 31, lty = 2)
axis(1, at=which(id_IFN_mGHS == "CYP27A1"), lwd.ticks = -1, labels="CYP27A1", las = 2)

segments(which(id_IFN_mGHS == "RELB"), y0 = 0, y1 = 31, lty = 2)
axis(1, at=which(id_IFN_mGHS == "RELB"), lwd.ticks = -1, labels="RELB", las = 2)

segments(which(id_IFN_mGHS == "SLC3A2"), y0 = 0, y1 = 31, lty = 2)
axis(1, at=which(id_IFN_mGHS == "SLC3A2"), lwd.ticks = -1, labels="SLC3A2", las = 2)

segments(which(id_IFN_mGHS == "TNK2"), y0 = 0, y1 = 31, lty = 2)
axis(1, at=which(id_IFN_mGHS == "TNK2"), lwd.ticks = -1, labels="SLC3A2", las = 2)



plot(1:381, mat_mGHS[2, ], type = "h", col = "black", xaxt='n', xlab = "", ylab = "number of connections", main = "LPS2h")

segments(which(id_LPS2_mGHS == "LYZ"), y0 = 0, y1 = 31, lty = 2)
axis(1, at=which(id_LPS2_mGHS == "LYZ"), lwd.ticks = -1, labels="LYZ", las = 2)

segments(which(id_LPS2_mGHS == "LGALS3"), lwd.ticks = -1, y0 = 0, y1 = 31, lty = 2)
axis(1, at=which(id_LPS2_mGHS == "LGALS3"), lwd.ticks = -1, labels="LGALS3", las = 2)

segments(which(id_LPS2_mGHS == "CYP27A1"), y0 = 0, y1 = 31, lty = 2)
axis(1, at=which(id_LPS2_mGHS == "CYP27A1"), lwd.ticks = -1, labels="CYP27A1", las = 2)

segments(which(id_LPS2_mGHS == "RELB"), y0 = 0, y1 = 31, lty = 2)
axis(1, at=which(id_LPS2_mGHS == "RELB"), lwd.ticks = -1, labels="RELB", las = 2)

segments(which(id_LPS2_mGHS == "SLC3A2"), y0 = 0, y1 = 31, lty = 2)
axis(1, at=which(id_LPS2_mGHS == "SLC3A2"), lwd.ticks = -1, labels="SLC3A2", las = 2)

segments(which(id_LPS2_mGHS == "TNK2"), y0 = 0, y1 = 31, lty = 2)
axis(1, at=which(id_LPS2_mGHS == "TNK2"), lwd.ticks = -1, labels="SLC3A2", las = 2)



plot(1:381, mat_mGHS[3, ], type = "h", col = "black", xaxt='n', xlab = "", ylab = "number of connections", main = "LPS24h")

segments(which(id_LPS24_mGHS == "LYZ"), y0 = 0, y1 = 31, lty = 2)
axis(1, at=which(id_LPS24_mGHS == "LYZ"), lwd.ticks = -1, labels="LYZ", las = 2)

segments(which(id_LPS24_mGHS == "LGALS3"), y0 = 0, y1 = 31, lty = 2)
axis(1, at=which(id_LPS24_mGHS == "LGALS3"), lwd.ticks = -1, labels="LGALS3", las = 2)

segments(which(id_LPS24_mGHS == "CYP27A1"), y0 = 0, y1 = 31, lty = 2)
axis(1, at=which(id_LPS24_mGHS == "CYP27A1"), lwd.ticks = -1, labels="CYP27A1", las = 2)

segments(which(id_LPS24_mGHS == "RELB"), y0 = 0, y1 = 31, lty = 2)
axis(1, at=which(id_LPS24_mGHS == "RELB"), lwd.ticks = -1, labels="RELB", las = 2)

segments(which(id_LPS24_mGHS == "SLC3A2"), y0 = 0, y1 = 31, lty = 2)
axis(1, at=which(id_LPS24_mGHS == "SLC3A2"), lwd.ticks = -1, labels="SLC3A2", las = 2)

segments(which(id_LPS24_mGHS == "TNK2"), y0 = 0, y1 = 31, lty = 2)
axis(1, at=which(id_LPS24_mGHS == "TNK2"), lwd.ticks = -1, labels="SLC3A2", las = 2)



plot(1:381, mat_mGHS[4, ], type = "h", col = "black", xaxt='n', xlab = "", ylab = "number of connections", main = "UNSTIM")

segments(which(id_UNSTIM_mGHS == "LYZ"), y0 = 0, y1 = 31, lty = 2)
axis(1, at=which(id_UNSTIM_mGHS == "LYZ"), lwd.ticks = -1, labels="LYZ", las = 2)

segments(which(id_UNSTIM_mGHS == "LGALS3"), y0 = 0, y1 = 31, lty = 2)
axis(1, at=which(id_UNSTIM_mGHS == "LGALS3"), lwd.ticks = -1, labels="LGALS3", las = 2)

segments(which(id_UNSTIM_mGHS == "CYP27A1"), y0 = 0, y1 = 31, lty = 2)
axis(1, at=which(id_UNSTIM_mGHS == "CYP27A1"), lwd.ticks = -1, labels="CYP27A1", las = 2)

segments(which(id_UNSTIM_mGHS == "RELB"), y0 = 0, y1 = 31, lty = 2)
axis(1, at=which(id_UNSTIM_mGHS == "RELB"), lwd.ticks = -1, labels="RELB", las = 2)

segments(which(id_UNSTIM_mGHS == "SLC3A2"), y0 = 0, y1 = 31, lty = 2)
axis(1, at=which(id_UNSTIM_mGHS == "SLC3A2"), lwd.ticks = -1, labels="SLC3A2", las = 2)

segments(which(id_UNSTIM_mGHS == "TNK2"), y0 = 0, y1 = 31, lty = 2)
axis(1, at=which(id_UNSTIM_mGHS == "TNK2"), lwd.ticks = -1, labels="SLC3A2", las = 2)




