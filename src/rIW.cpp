#include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>     

#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]

Eigen::MatrixXd rIW (const double& df, const Eigen::MatrixXd& S, const int& K)
{
  Eigen::MatrixXd I(K, K); I.setIdentity();
  Eigen::MatrixXd S_1 = S.llt().solve(I);
  
  Eigen::MatrixXd Z(K, K); Z.setZero();
  for (int hs = 0; hs < K; ++hs)
  {
    Z(hs, hs) = sqrt(R::rchisq(df - hs));
    // for (int h = 0; h < hs; ++h)
    //   Z(h, hs) = R::rnorm(0, 1);
  }
  for (int hs = 0; hs < K; ++hs)
  {
    for (int h = 0; h < hs; ++h)
      Z(h, hs) = R::rnorm(0, 1);
  }
  Eigen::MatrixXd L(S_1.llt().matrixU());
  Z = Z.triangularView<Eigen::Upper>() * L;
  return (Z.transpose() * Z).inverse();
}
