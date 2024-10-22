BIC_GemBag <- function(Theta_l, S_l, P_l, n_l) {
  ############################
  # Calculate BIC for GemBag #
  ############################
  
  k <- length(n_l)
  BIC <- 0
  for (i in 1:k){
    Theta <- Theta_l[[i]]
    S <- S_l[[i]]
    P <- P_l[[i]]
    n <- n_l[i]
    
    p_n <- nrow(Theta)
    P1 <- diag(1,p_n)
    P1[P>0.5] <- 1
    if(min(eigen(Theta*P1)$values) > 0){
      BIC <- BIC + n*(sum(diag(S%*%(Theta*P1)))-log(det(Theta*P1)))+(log(n))*(sum(P>0.5)/2)
    }else{
      return(Inf)
    }
  }
  return(BIC)
}

Tune_GemBag<-function(v0_l, v1_l, S_l, n_l, maxiter, p_1, p_2) {
  ####################################################
  # Select tuning parameters for GemBag based on BIC #
  ####################################################
  
  BIC <- matrix(NA, nrow=length(v0_l), ncol=length(v1_l))
  i = 1
  for (v0 in v0_l) {
    j = 1
    for (v1 in v1_l) {
      tmp <- GemBag(S_l=S_l, n=n_l, v_0=v0, v_1=v1, maxiter=maxiter, p_1=p_1, p_2=p_2, tau=v0)
      names(tmp) <- c('Theta', 'P', 'W')
      BIC[i, j] <- BIC_GemBag(tmp$Theta, S_l, tmp$P, n_l)
      j <- j + 1
    }
    i <- i + 1
  }
  idx <- which(BIC == min(BIC, na.rm = T), arr.ind = TRUE)
  return(list(v0 = v0_l[idx[1]], v1 = v1_l[idx[2]], val = BIC[idx[1], idx[2]], vals = BIC))
}
