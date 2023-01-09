#include <RcppEigen.h>
#include "Gamma2.h"
#include "VVLA.h"

#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE
#define pi 3.141592653589793

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppEigen)]]


double DVLA (const double& va, const double& x) 
{
  int i = 1;
  
  const double eps = 1.0e-12;
  double ep = std::exp(-0.25 * x * x);
  double a0 = std::pow(std::abs(x), va) * ep;
  double r = 1.0, pd = 1.0;
  
  while ((i <= 16) && (std::abs(r / pd) > eps))
  {
    r = -0.5 * r * (2.0 * i - va - 1.0) * (2.0 * i - va - 2.0) / (i * x * x);
    pd += r;
    i += 1;
  }
  
  pd *= a0;
  if (x < 0.0) 
  {
    double x1 = -x;
    double vl = VVLA(va, x1);
    double gl = gamma2(-va);
    
    pd = pi * vl / gl + std::cos(pi * va) * pd;
  }
  
  return pd;
}

