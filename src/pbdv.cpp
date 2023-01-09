#include <RcppEigen.h>
#include "Gamma2.h"
#include "DVLA.h"
#include "DVSA.h"
#include "VVLA.h"

#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE
#define pi 3.141592653589793

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppEigen)]]


template <typename T> int sign (T val) 
{
  return (T(0) < val) - (val < T(0));
}

// [[Rcpp::export]]
double pbdv (const double& va, const double& x) 
{
  // Function to compute Parabolic Cylinder function D_va(x)
  
  int i = 0;
  
  double xa = std::abs(x);
  double v = va + sign(va);
  int nv = (int)v;
  double v0 = v - (double)nv, v1;
  
  int na = std::abs(nv);
  double ep = std::exp(-0.25 * x * x);
  
  double pd0, pd1;
  double f = 0.0, f0 = 1.0e-30, f1 = 0.0;
  
  int ja = 0;
  if (na >= 1)
    ja = 1;
  
  Eigen::VectorXd dv(na+1); dv.setZero();
  double pdf;
  
  if (v >= 0.0) 
  {
    if (v0 == 0.0) 
    {
      pd0 = ep;
      pd1 = x * ep;  
    } 
    else 
    {
      v1 = v0;
      if (xa <= 5.8)
        pd1 = DVSA(v1, x);
      if (xa > 5.8)
        pd1 = DVLA(v1, x);
      pd0 = pd1;
      
      if (ja == 1) 
      {
        v1 = v0 + 1;
        if (xa <= 5.8)
          pd1 = DVSA(v1, x);
        if (xa > 5.8)
          pd1 = DVLA(v1, x);
      }
    }
    
    dv(0) = pd0;
    dv(1) = pd1;
    
    for (i = 2; i <= na; i++) 
    {
      pdf = x * pd1 - (i + v0 - 1.0) * pd0;
      dv(i) = pdf;
      pd0 = pd1;
      pd1 = pdf;
    }
  } 
  else 
  {
    if (x <= 0.0) 
    {
      if (xa <= 5.8) 
      {
        pd0 = DVSA(v0, x);
        v1 = v0 - 1.0;
        pd1 = DVSA(v1, x);
      } 
      else 
      {
        pd0 = DVLA(v0, x);
        v1 = v0 - 1.0;
        pd1 = DVLA(v1, x);
      }
      
      dv(0) = pd0;
      dv(1) = pd1;
      
      for (i = 2; i <= na; i++) 
      {
        pdf = (-x * pd1 + pd0) / (i - 1.0 - v0);
        dv(i) = pdf;
        pd0 = pd1;
        pd1 = pdf;
      }
    } 
    else 
    {
      if (x <= 2.0) 
      {
        double v2 = nv + v0;
        
        if (nv == 0.0)
          v2 = v2 - 1.0;
        
        int nk = int(-v2);
        f1 = DVSA(v2, x);
        v1 = v2 + 1.0;
        f0 = DVSA(v1, x);
        
        dv(nk) = f1;
        dv(nk-1) = f0;
        
        for (i = nk-2; i >= 0; i--) 
        {
          f = x * f0 + (i - v0 + 1.0) * f1;
          dv(i) = f;
          f1 = f0;
          f0 = f;
        }
      } 
      else 
      {
        if (xa <= 5.8) 
          pd0 = DVSA(v0, x);
        if (xa > 5.8)
          pd0 = DVLA(v0, x);
        
        dv(0) = pd0;
        int m = 100 + na;
        
        for (i = m; i >= 0; i--) 
        {
          f = x * f0 + (i - v0 + 1.0) * f1;
          if (i <= na)
            dv(i) = f;
          f1 = f0;
          f0 = f;
        }
        
        double s0 = pd0 / f;
        dv *= s0;
      }
    }
  }
  
  return dv(na-1);
}

