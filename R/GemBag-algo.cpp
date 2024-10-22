#include <RcppArmadillo.h>
#include <math.h>
#include <cstring>
#include <R.h>
#include <Rinternals.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace std;
using namespace arma;
using namespace Rcpp;

// Soft thresholding function
int sgn(double v) {
  return (v < 0) ? -1 : 1;
}

double SoftThreshold(double a,double b,double lambda, double c = 0){
  return -c+sgn(c-b/a)*fmax(abs(c-b/a)-lambda/a,0.0);
}

Mat<double> Cor_desc(
  Mat<double> Theta_11_inv,
  double W_22,
  Mat<double> Theta_12,
  Mat<double> S_12,
  double v_0,
  double v_1,
  Mat<double> P_12,
  int n 
){
  int p = Theta_11_inv.n_rows;
  Mat<double> Theta_12_2 = Theta_12;
  for( int iter = 0; iter < 1000; iter++){
    Mat<double> Theta_12_1 = Theta_12;
    for(int i = 0; i<p ; i++){
      double b = n*(as_scalar(Theta_11_inv.row(i)*Theta_12)*W_22- \
                        (Theta_11_inv(i,i)*Theta_12(i))*W_22+S_12(i));
      double a = n*Theta_11_inv(i,i)*W_22;
      double lambda = P_12(i)/v_1 + (1-P_12(i))/v_0;
      Theta_12(i) = SoftThreshold(a,b,lambda);
    }
    
    if((abs(Theta_12_1 - Theta_12)).max() < 0.001) return(Theta_12);
    if(abs(Theta_12).max() > 10) return(Theta_12_2);
  }
  return Theta_12_2;
}

// [[Rcpp::export]] 
// Main Function
List GemBag(
  List S_l,
  NumericVector n_l,
  double v_0,
  double v_1,
  int maxiter,
  double p_1,
  double p_2,
  double tau
){
  // Input:
  //   S_l: sample covariance matrices
  //   n_l: sample sizes of the classes
  //   v_0: spike prior parameter
  //   v_1: slab prior parameter
  //   maxiter: maximum of number of iterations
  //   p_1: parameter of Bernoulli prior on binary latent indicator $\gamma_{ij}$
  //   tau: parameter of exponential prior on diagonal entries of precision matrix
  // 
  //  Return: 
  //    Theta_l: estimated precision matrices
  //    P_l: estimated posterior inclusion probabilities
  //    W_l: estimated covariance matrices
  
  double p = as<arma::mat>(S_l[0]).n_rows;
  int k = S_l.size();
  List W_l(k);
  List Theta_l(k);
  for(int i = 0; i < k; i++) {
    W_l[i] = eye<mat>(p, p);
    Theta_l[i] = eye<mat>(p, p);
  } 
  
  for(int iter =0; iter<maxiter; iter++){
    // E-step
    // Calculate P_l
    mat tmp = ones(p, p);
    for(int i = 0; i < k; i++) {
      tmp = tmp % (p_2*v_0/v_1*exp(-abs(as<arma::mat>(Theta_l[i]))*(1.0/v_1-1.0/v_0)) + (1-p_2));
    }
    mat shrk = 1.0 / ( 1.0 + (1.0-p_1)/p_1 / tmp);
    
    List P_l(k);
    for(int i = 0; i < k; i++) {
      P_l[i] = shrk % (1.0/(1.0 + v_1/v_0*exp(-abs(as<arma::mat>(Theta_l[i]))/v_0 + abs(as<arma::mat>(Theta_l[i]))/v_1)*(1-p_2)/p_2));
    }
    
    for(int i = 0; i < k; i++) { 
      mat Theta = as<arma::mat>(Theta_l[i]);
      mat W = as<arma::mat>(W_l[i]);
      mat P = as<arma::mat>(P_l[i]);
      mat S = as<arma::mat>(S_l[i]);
      int n = n_l[i];
      
      for(int o=0; o<5; o++){
        // M-step
        for(int i = 0; i<p; i++){
          
          uvec idx(p);
          std::iota(std::begin(idx),std::end(idx),0);
          idx.shed_row(i);
          
          uvec indice(1);
          indice << i;
          Mat<double> W_11 = W.submat(idx, idx);
          
          Mat<double> W_12 = W.submat(idx, indice);
          double W_22 = W(i,i);
  
          
          Mat<double> Theta_12 = Theta.submat(idx, indice);
          Mat<double> Theta_11 = Theta(idx,idx);
          double Theta_22 = Theta(i,i);
          
          Mat<double> P_12 = P.submat(idx, indice);
          Mat<double> S_12 = S.submat(idx, indice);
          
          Mat<double> Theta_11_inv = W_11 - W_12*W_12.t()/W_22;
          W_22 = S(i,i) + 2/n*tau;
          Theta_12 = Cor_desc(Theta_11_inv, W_22, Theta_12, S_12, v_0, v_1, P_12, n);
          Theta_22 = 1.0/W_22 + as_scalar(Theta_12.t()*Theta_11_inv*Theta_12);
  
          Theta(indice, idx) = Theta_12.t();
          Theta(idx, indice) = Theta_12;
          Theta(i,i) = Theta_22;
          
          double tmp = as_scalar(Theta_22 - Theta_12.t()*Theta_11_inv*Theta_12);
          Mat<double> tmp1 = Theta_11_inv *  Theta_12;
          W_11 = Theta_11_inv + tmp1 * tmp1.t()/tmp;
          W_12 = -tmp1/tmp;
          
          W(idx,idx) = W_11;
          W(indice, idx) = W_12.t();
          W(idx, indice) = W_12;
          W(i,i)  = W_22;
        }
      }
      Theta_l[i] = Theta;
      W_l[i] = W;
      P_l[i] = P;
    }
  }
  // Calculate P_l
  mat tmp = ones(p, p);
  for(int i = 0; i < k; i++) {
    tmp = tmp % (p_2*v_0/v_1*exp(-abs(as<arma::mat>(Theta_l[i]))*(1.0/v_1-1.0/v_0)) + (1-p_2));
  }
  mat shrk = 1.0 / ( 1.0 + (1.0-p_1)/p_1 / tmp);
  
  List P_l(k);
  for(int i = 0; i < k; i++) {
    P_l[i] = shrk % (1.0/(1.0 + v_1/v_0*exp(-abs(as<arma::mat>(Theta_l[i]))/v_0 + abs(as<arma::mat>(Theta_l[i]))/v_1)*(1-p_2)/p_2));
  }
  
  List returnlist(3);
  returnlist(0) = Theta_l;
  returnlist(1) = P_l;
  returnlist(2) = W_l;

  return returnlist;
  
}

