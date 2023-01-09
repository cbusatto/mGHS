#include <RcppEigen.h>

#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppEigen)]]

long long factorial (const int& n)
{
  long long f = 1;
  for (int i = 1; i <= n; ++i)
    f *= i;
  return f;
}