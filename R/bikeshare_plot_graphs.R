load("~/res_bikeshare.RData")

p = 239
K = 6

mean_R = matrix(0, K, K)
for (b in 1:5000) mean_R = mean_R + res_mGHS$R[, ((b-1)*K+1):(b*K)]
mean_R = mean_R / 5000
mean_R

library(igraph)


###########################################################################################################
## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## PLOT GRAPHS WITH mGHS
## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
###########################################################################################################



adj_mGHS = array(0, c(p, p, K))
Prob = array(0, c(p, p, K))

mean_lambda = matrix(0, p*(p-1)/2, K)
for (b in 1:5000) mean_lambda = mean_lambda + 1 / (1+res_mGHS$lambda[, ((b-1)*K+1):(b*K)])
mean_lambda = mean_lambda / 5000

t = 0.8
cs = 1
for (k in 1:K) {
  cs = 1
  for (j in 1:(p-1)) {
    for (js in (j+1):p) {
      if (1-mean_lambda[cs, k] > t) {
        adj_mGHS[js, j, k] = adj_mGHS[j, js, k] = 1
      }
      
      cs = cs + 1
    }
  }
}

nadj_mGHS_c = adj_mGHS[, , 1] + adj_mGHS[, , 2] + adj_mGHS[, , 3]
nadj_mGHS_m = adj_mGHS[, , 4] + adj_mGHS[, , 5] + adj_mGHS[, , 6]

w_mGHS_c = nadj_mGHS_c[lower.tri(nadj_mGHS_c)]
w_mGHS_c[(w_mGHS_c == 2)] = 1

w_mGHS_m = nadj_mGHS_m[lower.tri(nadj_mGHS_m)]
w_mGHS_m[(w_mGHS_m == 2)] = 1


c1 = graph_from_adjacency_matrix(adj_mGHS[, , 1], mode = "undirected", diag = F)
V(c1)$name <- paste0("Node", 1:p)
g1 = adj_mGHS[, , 1]
g1 = g1[lower.tri(g1)]
E(c1)$weight = w_mGHS_c[g1 != 0]
E(c1)$color[E(c1)$weight == 3] <- 'black'
E(c1)$color[E(c1)$weight == 1] <- 'darkgrey'
# E(c1)$lty = w_mGHS_c[g1 != 0]

c2 = graph_from_adjacency_matrix(adj_mGHS[, , 2], mode = "undirected", diag = F)
V(c2)$name <- paste0("Node", 1:p)
g2 = adj_mGHS[, , 2]
g2 = g2[lower.tri(g2)]
E(c2)$weight = w_mGHS_c[g2 != 0]
E(c2)$color[E(c2)$weight == 3] <- 'black'
E(c2)$color[E(c2)$weight == 1] <- 'darkgrey'

c3 = graph_from_adjacency_matrix(adj_mGHS[, , 3], mode = "undirected", diag = F)
V(c3)$name <- paste0("Node", 1:p)
g3 = adj_mGHS[, , 3]
g3 = g3[lower.tri(g3)]
E(c3)$weight = w_mGHS_c[g3 != 0]
E(c3)$color[E(c3)$weight == 3] <- 'black'
E(c3)$color[E(c3)$weight == 1] <- 'darkgrey'


m1 = graph_from_adjacency_matrix(adj_mGHS[, , 4], mode = "undirected", diag = F)
V(m1)$name <- paste0("Node", 1:p)
g4 = adj_mGHS[, , 4]
g4 = g4[lower.tri(g4)]
E(m1)$weight = w_mGHS_m[g4 != 0]
E(m1)$color[E(m1)$weight == 3] <- 'black'
E(m1)$color[E(m1)$weight == 1] <- 'darkgrey'

m2 = graph_from_adjacency_matrix(adj_mGHS[, , 5], mode = "undirected", diag = F)
V(m2)$name <- paste0("Node", 1:p)
g5 = adj_mGHS[, , 5]
g5 = g5[lower.tri(g5)]
E(m2)$weight = w_mGHS_m[g5 != 0]
E(m2)$color[E(m2)$weight == 3] <- 'black'
E(m2)$color[E(m2)$weight == 1] <- 'darkgrey'

m3 = graph_from_adjacency_matrix(adj_mGHS[, , 6], mode = "undirected", diag = F)
V(m3)$name <- paste0("Node", 1:p)
g6 = adj_mGHS[, , 6]
g6 = g6[lower.tri(g6)]
E(m3)$weight = w_mGHS_m[g6 != 0]
E(m3)$color[E(m3)$weight == 3] <- 'black'
E(m3)$color[E(m3)$weight == 1] <- 'darkgrey'


v_size_c1 = adj_mGHS[, , 1]
v_size_c2 = adj_mGHS[, , 2]
v_size_c3 = adj_mGHS[, , 3]

v_size_m1 = adj_mGHS[, , 4]
v_size_m2 = adj_mGHS[, , 5]
v_size_m3 = adj_mGHS[, , 6]

v_size_c1 = apply(v_size_c1, 1, function(x) sum(x == 1))
v_size_c2 = apply(v_size_c2, 1, function(x) sum(x == 1))
v_size_c3 = apply(v_size_c3, 1, function(x) sum(x == 1))

v_size_m1 = apply(v_size_m1, 1, function(x) sum(x == 1))
v_size_m2 = apply(v_size_m2, 1, function(x) sum(x == 1))
v_size_m3 = apply(v_size_m3, 1, function(x) sum(x == 1))

set.seed(10)
layg1 <- layout_on_sphere(m1)

xlim <- range(c(layg1[,1], layg1[,1]))
ylim <- range(c(layg1[,2], layg1[,2]))

par(mfrow=c(2, 3))
plot(c1, vertex.size = v_size_c1, layout=layg1, xlim=xlim, ylim=ylim, rescale=FALSE, vertex.label = NA, vertex.color = "black", edge.color = E(c1)$color, main = "casual 2016")
plot(c2, vertex.size = v_size_c2, layout=layg1, xlim=xlim, ylim=ylim, rescale=FALSE, vertex.label = NA, vertex.color = "black", edge.color = E(c2)$color, main = "casual 2017")
plot(c3, vertex.size = v_size_c3, layout=layg1, xlim=xlim, ylim=ylim, rescale=FALSE, vertex.label = NA, vertex.color = "black", edge.color = E(c3)$color, main = "casual 2018")
plot(m1, vertex.size = v_size_m1, layout=layg1, xlim=xlim, ylim=ylim, rescale=FALSE, vertex.label = NA, vertex.color = "black", edge.color = E(m1)$color, main = "member 2016")
plot(m2, vertex.size = v_size_m2, layout=layg1, xlim=xlim, ylim=ylim, rescale=FALSE, vertex.label = NA, vertex.color = "black", edge.color = E(m2)$color, main = "member 2017")
plot(m3, vertex.size = v_size_m3, layout=layg1, xlim=xlim, ylim=ylim, rescale=FALSE, vertex.label = NA, vertex.color = "black", edge.color = E(m3)$color, main = "member 2018")



###########################################################################################################
## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## PLOT GRAPHS WITH GHS
## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
###########################################################################################################


mean_lambda_GHS = matrix(0, p*(p-1)/2, K)
for (b in 1:5000) {
  mean_lambda_GHS[, 1] = mean_lambda_GHS[, 1] + 1 / (1+res_GHS1[[2]][, b])
  mean_lambda_GHS[, 2] = mean_lambda_GHS[, 2] + 1 / (1+res_GHS2[[2]][, b])
  mean_lambda_GHS[, 3] = mean_lambda_GHS[, 3] + 1 / (1+res_GHS3[[2]][, b])
  mean_lambda_GHS[, 4] = mean_lambda_GHS[, 4] + 1 / (1+res_GHS4[[2]][, b])
  mean_lambda_GHS[, 5] = mean_lambda_GHS[, 5] + 1 / (1+res_GHS5[[2]][, b])
  mean_lambda_GHS[, 6] = mean_lambda_GHS[, 6] + 1 / (1+res_GHS6[[2]][, b])
  
  if (b %% 500 == 0) print(b)
}

mean_lambda_GHS = mean_lambda_GHS / 5000

adj_GHS = array(0, c(p, p, K))
cs = 1
for (k in 1:K) {
  cs = 1
  for (j in 1:(p-1)) {
    for (js in (j+1):p) {
      if (1-mean_lambda_GHS[cs, k] > t) {
        adj_GHS[js, j, k] = adj_GHS[j, js, k] = 1
      }
      
      cs = cs + 1
    }
  }
}

nadj_GHS_c = adj_GHS[, , 1] + adj_GHS[, , 2] + adj_GHS[, , 3]
nadj_GHS_m = adj_GHS[, , 4] + adj_GHS[, , 5] + adj_GHS[, , 6]

w_GHS_c = nadj_GHS_c[lower.tri(nadj_GHS_c)]
w_GHS_c[(w_GHS_c == 2)] = 1

w_GHS_m = nadj_GHS_m[lower.tri(nadj_GHS_m)]
w_GHS_m[(w_GHS_m == 2)] = 1


c1 = graph_from_adjacency_matrix(adj_GHS[, , 1], mode = "undirected", diag = F)
V(c1)$name <- paste0("Node", 1:p)
g1 = adj_GHS[, , 1]
g1 = g1[lower.tri(g1)]
E(c1)$weight = w_GHS_c[g1 != 0]
E(c1)$color[E(c1)$weight == 3] <- 'black'
E(c1)$color[E(c1)$weight == 1] <- 'darkgrey'
# E(c1)$lty = w_GHS_c[g1 != 0]

c2 = graph_from_adjacency_matrix(adj_GHS[, , 2], mode = "undirected", diag = F)
V(c2)$name <- paste0("Node", 1:p)
g2 = adj_GHS[, , 2]
g2 = g2[lower.tri(g2)]
E(c2)$weight = w_GHS_c[g2 != 0]
E(c2)$color[E(c2)$weight == 3] <- 'black'
E(c2)$color[E(c2)$weight == 1] <- 'darkgrey'

c3 = graph_from_adjacency_matrix(adj_GHS[, , 3], mode = "undirected", diag = F)
V(c3)$name <- paste0("Node", 1:p)
g3 = adj_GHS[, , 3]
g3 = g3[lower.tri(g3)]
E(c3)$weight = w_GHS_c[g3 != 0]
E(c3)$color[E(c3)$weight == 3] <- 'black'
E(c3)$color[E(c3)$weight == 1] <- 'darkgrey'


m1 = graph_from_adjacency_matrix(adj_GHS[, , 4], mode = "undirected", diag = F)
V(m1)$name <- paste0("Node", 1:p)
g4 = adj_GHS[, , 4]
g4 = g4[lower.tri(g4)]
E(m1)$weight = w_GHS_m[g4 != 0]
E(m1)$color[E(m1)$weight == 3] <- 'black'
E(m1)$color[E(m1)$weight == 1] <- 'darkgrey'

m2 = graph_from_adjacency_matrix(adj_GHS[, , 5], mode = "undirected", diag = F)
V(m2)$name <- paste0("Node", 1:p)
g5 = adj_GHS[, , 5]
g5 = g5[lower.tri(g5)]
E(m2)$weight = w_GHS_m[g5 != 0]
E(m2)$color[E(m2)$weight == 3] <- 'black'
E(m2)$color[E(m2)$weight == 1] <- 'darkgrey'

m3 = graph_from_adjacency_matrix(adj_GHS[, , 6], mode = "undirected", diag = F)
V(m3)$name <- paste0("Node", 1:p)
g6 = adj_GHS[, , 6]
g6 = g6[lower.tri(g6)]
E(m3)$weight = w_GHS_m[g6 != 0]
E(m3)$color[E(m3)$weight == 3] <- 'black'
E(m3)$color[E(m3)$weight == 1] <- 'darkgrey'


v_size_c1 = adj_GHS[, , 1]
v_size_c2 = adj_GHS[, , 2]
v_size_c3 = adj_GHS[, , 3]

v_size_m1 = adj_GHS[, , 4]
v_size_m2 = adj_GHS[, , 5]
v_size_m3 = adj_GHS[, , 6]

v_size_c1 = apply(v_size_c1, 1, function(x) sum(x == 1))
v_size_c2 = apply(v_size_c2, 1, function(x) sum(x == 1))
v_size_c3 = apply(v_size_c3, 1, function(x) sum(x == 1))

v_size_m1 = apply(v_size_m1, 1, function(x) sum(x == 1))
v_size_m2 = apply(v_size_m2, 1, function(x) sum(x == 1))
v_size_m3 = apply(v_size_m3, 1, function(x) sum(x == 1))

# set.seed(10)
# layg1 <- layout_on_sphere(m1)
# 
# xlim <- range(c(layg1[,1], layg1[,1]))
# ylim <- range(c(layg1[,2], layg1[,2]))

par(mfrow=c(2, 3))
plot(c1, vertex.size = v_size_c1, layout=layg1, xlim=xlim, ylim=ylim, rescale=FALSE, vertex.label = NA, vertex.color = "black", edge.color = E(c1)$color, main = "casual 2016")
plot(c2, vertex.size = v_size_c2, layout=layg1, xlim=xlim, ylim=ylim, rescale=FALSE, vertex.label = NA, vertex.color = "black", edge.color = E(c2)$color, main = "casual 2017")
plot(c3, vertex.size = v_size_c3, layout=layg1, xlim=xlim, ylim=ylim, rescale=FALSE, vertex.label = NA, vertex.color = "black", edge.color = E(c3)$color, main = "casual 2018")
plot(m1, vertex.size = v_size_m1, layout=layg1, xlim=xlim, ylim=ylim, rescale=FALSE, vertex.label = NA, vertex.color = "black", edge.color = E(m1)$color, main = "member 2016")
plot(m2, vertex.size = v_size_m2, layout=layg1, xlim=xlim, ylim=ylim, rescale=FALSE, vertex.label = NA, vertex.color = "black", edge.color = E(m2)$color, main = "member 2017")
plot(m3, vertex.size = v_size_m3, layout=layg1, xlim=xlim, ylim=ylim, rescale=FALSE, vertex.label = NA, vertex.color = "black", edge.color = E(m3)$color, main = "member 2018")



###########################################################################################################
## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## PLOT INTERSECTIONS
## ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
###########################################################################################################


Omega_mGHS = array(0, c(p, p, K))

cs = 1
for (k in 1:K) {
  diag(Omega_mGHS[, , k]) = res_mGHS$Omega_jj[, k]
  cs = 1
  for (j in 1:(p-1)) {
    for (js in (j+1):p) {
      Omega_mGHS[js, j, k] = Omega_mGHS[j, js, k] = res_mGHS$Omega_ij[cs, k]
      cs = cs + 1
    }
  }
}


Omega_GHS = array(0, c(p, p, K))
Omega_GHS[, , 1] = apply(res_GHS1[[1]], c(1, 2), mean)
Omega_GHS[, , 2] = apply(res_GHS2[[1]], c(1, 2), mean)
Omega_GHS[, , 3] = apply(res_GHS3[[1]], c(1, 2), mean)
Omega_GHS[, , 4] = apply(res_GHS4[[1]], c(1, 2), mean)
Omega_GHS[, , 5] = apply(res_GHS5[[1]], c(1, 2), mean)
Omega_GHS[, , 6] = apply(res_GHS6[[1]], c(1, 2), mean)


Omega_mGHS_c = (Omega_mGHS[, , 1] + Omega_mGHS[, , 2] + Omega_mGHS[, , 3]) / 3
Omega_mGHS_m = (Omega_mGHS[, , 4] + Omega_mGHS[, , 5] + Omega_mGHS[, , 6]) / 3

Omega_GHS_c = (Omega_GHS[, , 1] + Omega_GHS[, , 2] + Omega_GHS[, , 3]) / 3
Omega_GHS_m = (Omega_GHS[, , 4] + Omega_GHS[, , 5] + Omega_GHS[, , 6]) / 3


adj_mGHS = array(0, c(p, p, K))
w_mGHS = matrix(0, 0.5 * p * (p-1), K)

cs = 1
for (k in 1:K) {
  cs = 1
  for (j in 1:(p-1)) {
    for (js in (j+1):p) {
      if (1-mean_lambda[cs, k] > t) {
        
        adj_mGHS[js, j, k] = adj_mGHS[j, js, k] = 1
        w_mGHS[cs, k] = res_mGHS$Omega_ij[cs, k]
      }
      
      cs = cs + 1
    }
  }
}

adj_GHS = array(0, c(p, p, K))
w_GHS = matrix(0, 0.5 * p * (p-1), K)

for (k in 1:K) {
  cs = 1
  for (j in 1:(p-1)) {
    for (js in (j+1):p) {
      if (1-mean_lambda_GHS[cs, k] > t) {
        
        adj_GHS[js, j, k] = adj_GHS[j, js, k] = 1
        w_GHS[cs, k] = Omega_GHS[js, j, k]
      }
      
      cs = cs + 1
    }
  }
}



adj_mGHS_c = matrix(0, p, p)
adj_mGHS_c[((adj_mGHS[, , 1] == 1) & (adj_mGHS[, , 2] == 1) & (adj_mGHS[, , 3] == 1))] = 1

adj_mGHS_m = matrix(0, p, p)
adj_mGHS_m[((adj_mGHS[, , 4] == 1) & (adj_mGHS[, , 5] == 1) & (adj_mGHS[, , 6] == 1))] = 1

Omega_mGHS_c2 = matrix(0, p, p)
Omega_mGHS_m2 = matrix(0, p, p)
Omega_mGHS_c2[((adj_mGHS[, , 1] == 1) & (adj_mGHS[, , 2] == 1) & (adj_mGHS[, , 3] == 1))] = Omega_mGHS_c[((adj_mGHS[, , 1] == 1) & (adj_mGHS[, , 2] == 1) & (adj_mGHS[, , 3] == 1))]
Omega_mGHS_m2[((adj_mGHS[, , 4] == 1) & (adj_mGHS[, , 5] == 1) & (adj_mGHS[, , 6] == 1))] = Omega_mGHS_m[((adj_mGHS[, , 4] == 1) & (adj_mGHS[, , 5] == 1) & (adj_mGHS[, , 6] == 1))]

sum(adj_mGHS_c == 1)
sum(adj_mGHS_m == 1)


adj_GHS_c = matrix(0, p, p)
adj_GHS_c[((adj_GHS[, , 1] == 1) & (adj_GHS[, , 2] == 1) & (adj_GHS[, , 3] == 1))] = 1

adj_GHS_m = matrix(0, p, p)
adj_GHS_m[((adj_GHS[, , 4] == 1) & (adj_GHS[, , 5] == 1) & (adj_GHS[, , 6] == 1))] = 1

Omega_GHS_c2 = matrix(0, p, p)
Omega_GHS_m2 = matrix(0, p, p)
Omega_GHS_c2[((adj_GHS[, , 1] == 1) & (adj_GHS[, , 2] == 1) & (adj_GHS[, , 3] == 1))] = Omega_GHS_c[((adj_GHS[, , 1] == 1) & (adj_GHS[, , 2] == 1) & (adj_GHS[, , 3] == 1))]
Omega_GHS_m2[((adj_GHS[, , 4] == 1) & (adj_GHS[, , 5] == 1) & (adj_GHS[, , 6] == 1))] = Omega_GHS_m[((adj_GHS[, , 4] == 1) & (adj_GHS[, , 5] == 1) & (adj_GHS[, , 6] == 1))]

sum(adj_GHS_c == 1)
sum(adj_GHS_m == 1)


c_mGHS = graph_from_adjacency_matrix(adj_mGHS_c, mode = "undirected", diag = F)
w_mGHS_c = Omega_mGHS_c2[lower.tri(Omega_mGHS_c2)]
w_mGHS_c = w_mGHS_c[w_mGHS_c != 0.0]
E(c_mGHS)$weight = w_mGHS_c
V(c_mGHS)$name <- paste0("Node", 1:p)
E(c_mGHS)$color<-as.character(cut(w_mGHS_c, 
                                  breaks=c(-11, -0.1, 0.1, 5), 
                                  # labels=c("darkgrey", "#F0F0F0", "red")))
                                  labels=c("blue", "lightgrey", "red")))

m_mGHS = graph_from_adjacency_matrix(adj_mGHS_m, mode = "undirected", diag = F)
w_mGHS_m = Omega_mGHS_m2[lower.tri(Omega_mGHS_m2)]
w_mGHS_m = w_mGHS_m[w_mGHS_m != 0.0]
E(m_mGHS)$weight = w_mGHS_m
V(m_mGHS)$name <- paste0("Node", 1:p)
E(m_mGHS)$color<-as.character(cut(w_mGHS_m, 
                                  breaks=c(-11, -0.1, 0.1, 5), 
                                  # labels=c("darkgrey", "#F0F0F0", "red")))
                                  labels=c("blue", "lightgrey", "red")))

c_GHS = graph_from_adjacency_matrix(adj_GHS_c, mode = "undirected", diag = F)
w_GHS_c = Omega_GHS_c2[lower.tri(Omega_GHS_c2)]
w_GHS_c = w_GHS_c[w_GHS_c != 0.0]
E(c_GHS)$weight = w_GHS_c
V(c_GHS)$name <- paste0("Node", 1:p)
E(c_GHS)$color<-as.character(cut(w_GHS_c, 
                                 breaks=c(-11, -0.1, 0.1, 5), 
                                 # labels=c("darkgrey", "#F0F0F0", "red")))
                                 labels=c("blue", "lightgrey", "red")))

m_GHS = graph_from_adjacency_matrix(adj_GHS_m, mode = "undirected", diag = F)
w_GHS_m = Omega_GHS_m2[lower.tri(Omega_GHS_m2)]
w_GHS_m = w_GHS_m[w_GHS_m != 0.0]
E(m_GHS)$weight = w_GHS_m
V(m_GHS)$name <- paste0("Node", 1:p)
E(m_GHS)$color<-as.character(cut(w_GHS_m, 
                                 breaks=c(-11, -0.1, 0.1, 5), 
                                 # labels=c("darkgrey", "#F0F0F0", "red")))
                                 labels=c("blue", "lightgrey", "red")))

v_size_mGHS_m = apply(adj_mGHS_m, 1, function(x) sum(x == 1))
v_size_GHS_m = apply(adj_GHS_m, 1, function(x) sum(x == 1))
v_size_mGHS_c = apply(adj_mGHS_c, 1, function(x) sum(x == 1))
v_size_GHS_c = apply(adj_GHS_c, 1, function(x) sum(x == 1))

sum(adj_mGHS_m != 0)
sum(adj_GHS_m != 0)

# set.seed(10)
# layg1 <- layout_on_sphere(c_mGHS)
# layg1 = layout_on_sphere(m_mGHS)
# layg1 = layout_with_fr(c_mGHS, weight = NA)
# layg1 = layout_with_fr(m_mGHS, weight = NA)


# xlim <- range(c(layg1[,1], layg1[,1]))
# ylim <- range(c(layg1[,2], layg1[,2]))

par(mfrow=c(1, 1))
plot(m_mGHS, vertex.size = v_size_mGHS_m * 3, layout=layg1, xlim=xlim, ylim=ylim, rescale=FALSE, vertex.label = NA, edge.color = "darkgrey", vertex.color = "black", main = "")
plot(c_mGHS, vertex.size = v_size_mGHS_c * 3, layout=layg1, xlim=xlim, ylim=ylim, rescale=FALSE, vertex.label = NA, edge.color = "darkgrey", vertex.color = "black", main = "")
plot(m_GHS, vertex.size = v_size_GHS_m * 3, layout=layg1, xlim=xlim, ylim=ylim, rescale=FALSE, vertex.label = NA, edge.color = "darkgrey", vertex.color = "black", main = "")
plot(c_GHS, vertex.size = v_size_GHS_c * 3, layout=layg1, xlim=xlim, ylim=ylim, rescale=FALSE, vertex.label = NA, edge.color = "darkgrey", vertex.color = "black", main = "")


par(mfrow=c(2, 2), mar=c(0,2,2,0))
plot(m_mGHS, vertex.size = v_size_mGHS_m * 3, layout=layg1, xlim=xlim, ylim=ylim, rescale=FALSE, vertex.label = NA, edge.color = "darkgrey", vertex.color = "black", main = "member mGHS")
plot(c_mGHS, vertex.size = v_size_mGHS_c * 3, layout=layg1, xlim=xlim, ylim=ylim, rescale=FALSE, vertex.label = NA, edge.color = "darkgrey", vertex.color = "black", main = "casual mGHS")
plot(m_GHS, vertex.size = v_size_GHS_m * 3, layout=layg1, xlim=xlim, ylim=ylim, rescale=FALSE, vertex.label = NA, edge.color = "darkgrey", vertex.color = "black", main = "member GHS")
plot(c_GHS, vertex.size = v_size_GHS_c * 3, layout=layg1, xlim=xlim, ylim=ylim, rescale=FALSE, vertex.label = NA, edge.color = "darkgrey", vertex.color = "black", main = "casual GHS")

par(mfrow=c(2, 2), mar=c(0,2,2,0))
plot(m_mGHS, vertex.size = v_size_mGHS_m * 4, layout=layg1, xlim=xlim, ylim=ylim, rescale=FALSE, vertex.label = NA, edge.color = E(m_mGHS)$color, vertex.color = "black", main = "member mGHS")
plot(c_mGHS, vertex.size = v_size_mGHS_c * 4, layout=layg1, xlim=xlim, ylim=ylim, rescale=FALSE, vertex.label = NA, edge.color = E(c_mGHS)$color, vertex.color = "black", main = "casual mGHS")
plot(m_GHS, vertex.size = v_size_GHS_m * 4, layout=layg1, xlim=xlim, ylim=ylim, rescale=FALSE, vertex.label = NA, edge.color = E(m_GHS)$color, vertex.color = "black", main = "member GHS")
plot(c_GHS, vertex.size = v_size_GHS_c * 4, layout=layg1, xlim=xlim, ylim=ylim, rescale=FALSE, vertex.label = NA, edge.color = E(c_GHS)$color, vertex.color = "black", main = "casual GHS")
