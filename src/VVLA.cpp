#include <RcppEigen.h>
#include "Gamma2.h"
#include "DVLA.h"

#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE
#define pi 3.141592653589793

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppEigen)]]


double VVLA (const double& va, const double& x) 
{
  int i = 1;
  
  const double eps = 1.0e-12;
  double qe = std::exp(0.25*x*x);
  double a0 = std::pow(std::abs(x), -va-1.0) * std::sqrt(2.0 / pi) * qe;
  double r = 1.0, pv = 1.0;

  while ((i <= 18) && (std::abs(r / pv) > eps)) 
  {
    r = 0.5 * r * (2.0 * i + va - 1.0) * (2.0 * i + va) / (i * x * x);
    pv += r;
    i += 1;
  }
  
  pv *= a0;
  if (x < 0.0) 
  {
    double x1 = -x;
    double pdl = DVLA(va, x1);
    double gl = gamma2(-va);
    
    double dsl = std::sin(pi * va) * std::sin(pi * va);
    pv = dsl * gl / pi * pdl - std::cos(pi * va) * pv;
  }
  
  return pv;
}

