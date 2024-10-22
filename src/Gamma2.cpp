#include <RcppEigen.h>
#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE
#define pi 3.141592653589793

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppEigen)]]


double gamma2 (const double& x) 
{
  // Gamma function (accept negative values)
  
  int i;
  double ga = 0.0;
  
  if (x == (int)x) 
  {
    if (x > 0) 
    {
      ga = 1.0;
      double m1 = x - 1.0;
      for (i = 2; i <= m1; i++)
        ga *= i;
    } 
    else 
    {
      ga = 1.0e+300;
    }
  } 
  else 
  {
    double r = 1.0;
    double z = 0.0;
    if (std::abs(x) > 1.0) 
    {
      z = std::abs(x);
      int m = (int)z;
      for (i = 1; i <= m; i++)
        r *= (z - i);
      z -= m;
    } 
    else 
    {
      z = x;
    }
    
    Eigen::VectorXd g(26);
    g << 1.0, 0.5772156649015329, -0.6558780715202538, -0.0420026350340952, 0.1665386113822915, -0.0421977345555443, -0.96219715278770e-2, 0.72189432466630e-2, -0.11651675918591e-2, -0.2152416741149e-3, 0.1280502823882e-3, -0.201348547807e-4, -0.12504934821e-5, 0.11330272320e-5, -0.2056338417e-6, 0.61160950e-8, 0.50020075e-8, -0.11812746e-8, 0.1043427e-9, 0.77823e-11, -0.36968e-11, 0.51e-12, -0.206e-13, -0.54e-14, 0.14e-14, 0.1e-15;
    
    double gr = g(25);
    for (i = 24; i >= 0; i--)
      gr = gr * z + g(i);
    
    ga = 1.0 / (gr * z);
    
    if (std::abs(x) > 1.0) 
    {
      ga *= r;
      if (x < 0.0)
        ga = -pi / (x * ga * std::sin(pi * x));
    }
  }
  
  return ga;
}


