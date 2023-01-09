#include <RcppEigen.h>
#include "Gamma2.h"

#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE
#define pi 3.141592653589793

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppEigen)]]


double DVSA (const double& va, const double& x) 
{
  const double eps = 1.0e-15;
  double ep = std::exp(-0.25 * x * x);
  double pd;
  
  if (va == 0.0) 
  {
    pd = ep;
  } 
  else 
  {
    double sq2 = std::sqrt(2.0);
    double va0 = 0.5 * (1.0 - va);  
    
    if (x == 0.0) 
    {
      if ((va0 <= 0.0) && (va0 == int(va0))) 
      {
        pd = 0.0;
      } 
      else 
      {
        double ga0 = gamma2(va0);
        pd = std::sqrt(pi) / (ga0 * std::pow(2.0, -0.5 * va));
      }
    } 
    else 
    {
      int i = 1;
      double g1 = gamma2(-va);
      double vt = -0.5 * va;
      double a0 = ep / g1 * std::pow(2.0, vt - 1.0);
      double g0 = gamma2(vt);
      pd = g0;
      double r = 1.0, r1 = 1.0, vm = 0.0, gm = 0.0;
      
      while ((i <= 250) && (std::abs(r1) > std::abs(pd)*eps))
      {
        vm = 0.5 * (i-va);
        gm = gamma2(vm);
        r = -r * sq2 * x / i;
        r1 = gm * r;
        pd += r1;
        
        i += 1;
      }
      
      pd *= a0;
    }
  }
  
  return pd;
}