#include <Rcpp.h>
#include <RcppEigen.h>
#include <RcppParallel.h>
#include <math.h>     
// #include <Rcpp/Benchmark/Timer.h>

#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]] 
// [[Rcpp::depends(RcppParallel)]]

#include "boost/multi_array.hpp"
#include "utils_rgamma3p.h"

typedef boost::multi_array<double, 1> vector_type;
typedef boost::multi_array<double, 2> matrix_type;
typedef boost::multi_array<double, 3> array_type;
typedef array_type::index_range range;

#define s2 1.414213562373095
#define s2pi 2.506628274631


Eigen::VectorXd inverse (Eigen::VectorXd x, const int& n)
{
  for (int j = 0; j < n; ++j)
    x(j) = 1.0 / x(j);
  return x;
}


Eigen::VectorXd sinverse (Eigen::VectorXd x, const int& n)
{
  for (int j = 0; j < n; ++j)
    x(j) = 1.0 / sqrt(x(j));
  
  return x;
}

Eigen::VectorXd delete_elem_from_boost_vec (const vector_type& x, const int& j, const int& nm1)
{
  // delete an element j from a (nm1+1)-dimensional boos::vector_type x, return Eigen::VectorXd x[-j]
  
  Eigen::VectorXd out(nm1);
  int is = 0;
  for (int js = 0; js < nm1+1; ++js) 
  {
    if(js != j) 
    {
      out(is) = x[js];
      is += 1;
    }
  }
  
  return out;
}


Eigen::VectorXd delete_elem_from_eigen_vec (const Eigen::VectorXd& x, const int& j, const int& nm1)
{
  // delete an element j from a (nm1+1)-dimensional Eigen::VectorXd x, return Eigen::VectorXd x[-j]
  
  if (j == 0) 
  {
    return x.tail(nm1);
  } 
  else if (j == nm1) 
  {
    return x.head(nm1);
  } 
  else 
  {
    Eigen::VectorXd a(nm1);
    a.head(j) = x.head(j);
    a.tail(nm1-j) = x.tail(nm1-j);
    
    return a;
  }
}


Eigen::MatrixXd delete_rc_from_eigen (const Eigen::MatrixXd& X, const int& j, const int& pm1) 
{
  // delete row j and col j from a (pm1+1) x (pm1+1) Eigen::MatrixXd X, return Eigen::MatrixXd X[-j, -j]
  
  if (j == 0)
  {
    return X.bottomRightCorner(pm1, pm1);
  } 
  else if (j == pm1)
  {
    return X.topLeftCorner(pm1, pm1);
  } 
  else 
  {
    Eigen::MatrixXd Xs(pm1, pm1);
    
    Xs.topLeftCorner(j, j) = X.topLeftCorner(j, j);
    Xs.bottomRightCorner(pm1-j, pm1-j) = X.bottomRightCorner(pm1-j, pm1-j);
    Xs.topRightCorner(j, pm1-j) = X.topRightCorner(j, pm1-j);
    Xs.bottomLeftCorner(pm1-j, j) = X.bottomLeftCorner(pm1-j, j);
    
    return Xs;
  }
}


Eigen::MatrixXd delete_rc_from_boost1 (const matrix_type& X, const int& j, const int& pm1) 
{
  // delete row j and col j from a (pm1+1) x (pm1+1) boost::matrix_type X, return Eigen::MatrixXd X[-j, -j]
  
  const int p = pm1 + 1;
  
  matrix_type Xs_b(boost::extents[pm1][pm1]);
  if (j == 0)
  {
    Xs_b = X[boost::indices[range(1, p)][range(1, p)]];
  } 
  else if (j == pm1)
  {
    Xs_b = X[boost::indices[range(0, pm1)][range(0, pm1)]];
  } 
  else 
  {
    Xs_b[boost::indices[range(0, j)][range(0, j)]] = X[boost::indices[range(0, j)][range(0, j)]];
    Xs_b[boost::indices[range(0, j)][range(j, pm1)]] = X[boost::indices[range(0, j)][range(j+1, p)]];
    Xs_b[boost::indices[range(j, pm1)][range(0, j)]] = X[boost::indices[range(j+1, p)][range(0, j)]];
    Xs_b[boost::indices[range(j, pm1)][range(j, pm1)]] = X[boost::indices[range(j+1, p)][range(j+1, p)]];
  }
  
  Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > Xs_e(Xs_b.data(), pm1, pm1);
  return Xs_e;  
}


Eigen::MatrixXd delete_rc_from_boost2 (matrix_type& X, const int& j, const int& pm1) 
{
  // delete row j and col j from a (pm1+1) x (pm1+1) boost::matrix_type X, return Eigen::MatrixXd X[-j, -j]
  const int p = pm1 + 1;
  Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > Xs_e(X.origin(), p, p);
  
  if (j == 0)
  {
    return Xs_e.bottomRightCorner(pm1, pm1);
  } 
  else if (j == pm1)
  {
    return Xs_e.topLeftCorner(pm1, pm1);
  } 
  else 
  {
    Eigen::MatrixXd Xs(pm1, pm1);
    
    Xs.topLeftCorner(j, j) = Xs_e.topLeftCorner(j, j);
    Xs.bottomRightCorner(pm1-j, pm1-j) = Xs_e.bottomRightCorner(pm1-j, pm1-j);
    Xs.topRightCorner(j, pm1-j) = Xs_e.topRightCorner(j, pm1-j);
    Xs.bottomLeftCorner(pm1-j, j) = Xs_e.bottomLeftCorner(pm1-j, j);
    
    return Xs;
  }
}


Eigen::MatrixXd get_matrix (const int& j, const matrix_type& X, const int& p, const int& K) 
{
  // get Eigen::MatrixXd X[-j, j, ] from a p x p x K boost::array_type
  
  int i, js, h;
  
  Eigen::MatrixXd out(p-1, K); out.setZero();
  for (h = 0; h < K; ++h) 
  {
    i = 0;
    for (js = 0; js < p; ++js) 
    {
      if (js != j) 
      {
        out(i, h) = X[js][h];
        i += 1;
      }
    }
  }
  
  return out;
}


Eigen::MatrixXd get_Sigma_h (matrix_type& Sigma, const int& p) 
{
  // get Eigen::MatrixXd X[-j, j, ] from a p x p x K boost::array_type
  
  Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > Sigma_h(Sigma.origin(), p, p);
  
  return Sigma_h;
}


Eigen::MatrixXd get_delta (const Eigen::VectorXd& tau, const Eigen::MatrixXd& lambda_12, const int& pm1, const int& K)
{
  // get Eigen::VectorXd e_ij = tau * lambda (element-wise product)
  Eigen::MatrixXd out(pm1, K);
  for (int h = 0; h < K; ++h)
    out.col(h) = lambda_12.col(h) * tau(h);
  
  return out;
}


Eigen::VectorXd get_s2 (const vector_type& s_h, const int& p)
{
  Eigen::VectorXd out(p);
  std::copy(s_h.origin(), s_h.origin() + p, out.data());
  return out;
}


matrix_type fill_S (const Eigen::MatrixXd& S, const int& p)
{
  matrix_type X_b(boost::extents[p][p]);
  std::copy(S.data(), S.data() + p * p, X_b.origin());
  
  return X_b;
}


matrix_type fill_Xt (const Eigen::MatrixXd& Xt, const int& n, const int& j_train)
{
  matrix_type X_b(boost::extents[n][j_train]);
  
  for (int i = 0; i < n; i++)
    for (int j = 0; j < j_train; j++)
      X_b[i][j] = Xt(i, j);
  
  return X_b;
}


matrix_type get_Sigma (boost::multi_array_ref<double, 3>::array_view<2>::type Sigma, const Eigen::MatrixXd& inv_Omega_11, const int& j, const int& pm1, const int& p)
{
  matrix_type X_b(boost::extents[pm1][pm1]);
  std::copy(inv_Omega_11.data(), inv_Omega_11.data() + pm1 * pm1, X_b.origin());
  
  if (j == 0)
  {
    Sigma[boost::indices[range(1, p)][range(1, p)]] = X_b;
  } 
  else if (j == pm1)
  {
    Sigma[boost::indices[range(0, pm1)][range(0, pm1)]] = X_b;
  } 
  else 
  {
    Sigma[boost::indices[range(0, j)][range(0, j)]] = X_b[boost::indices[range(0, j)][range(0, j)]];
    Sigma[boost::indices[range(0, j)][range(j+1, p)]] = X_b[boost::indices[range(0, j)][range(j, pm1)]];
    Sigma[boost::indices[range(j+1, p)][range(0, j)]] = X_b[boost::indices[range(j, pm1)][range(0, j)]];
    Sigma[boost::indices[range(j+1, p)][range(j+1, p)]] = X_b[boost::indices[range(j, pm1)][range(j, pm1)]];
  }
  
  return Sigma;
}


Eigen::MatrixXd get_Sigma2(Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >& Sigma, const Eigen::MatrixXd& inv_Omega_11, const Eigen::VectorXd& v, const double& v22, const int& j, const int& pm1)
{
  if (j == 0)
  {
    Sigma(0, 0) = v22;
    Sigma.bottomRightCorner(pm1, pm1) = inv_Omega_11;
    Sigma.topRightCorner(1, pm1) = v.transpose();
    Sigma.bottomLeftCorner(pm1, 1) = v;
  } 
  else if (j == pm1)
  {
    Sigma(pm1, pm1) = v22;
    Sigma.topLeftCorner(pm1, pm1) = inv_Omega_11;
    Sigma.topRightCorner(pm1, 1) = v;
    Sigma.bottomLeftCorner(1, pm1) = v.transpose();
  } 
  else 
  {
    
    Sigma(j, j) = v22;
    Sigma.topLeftCorner(j, j) = inv_Omega_11.topLeftCorner(j, j);
    Sigma.bottomRightCorner(pm1-j, pm1-j) = inv_Omega_11.bottomRightCorner(pm1-j, pm1-j);
    Sigma.topRightCorner(j, pm1-j) = inv_Omega_11.topRightCorner(j, pm1-j);
    Sigma.bottomLeftCorner(pm1-j, j) = inv_Omega_11.bottomLeftCorner(pm1-j, j);
    
    Sigma.block(0, j, j, 1) = v.head(j);
    Sigma.block(j, 0, 1, j) = v.head(j).transpose();
    Sigma.block(j+1, j, pm1-j, 1) = v.tail(pm1-j);
    Sigma.block(j, j+1, 1, pm1-j) = v.tail(pm1-j).transpose();
  }
  
  return Sigma;
}


Eigen::VectorXd get_eij (const Eigen::VectorXd& tau, const vector_type& lambda, const Eigen::VectorXd& diag_V, const vector_type& omega, const int& K)
{
  // get Eigen::VectorXd e_ij = tau * lambda (element-wise product)
  
  Eigen::VectorXd eij(K);
  for (int i = 0; i < K; ++i)
    eij(i) = diag_V(i) * omega[i] / sqrt(tau(i) * lambda[i]);
  
  return eij;
}


Eigen::VectorXd get_diag_vec_h (const Eigen::VectorXd& tau, const vector_type& lambda, const int& h, const int& Km1)
{
  // get Eigen::VectorXd x = tau * lambda without the h-th element (element-wise product)
  
  Eigen::VectorXd out(Km1);
  int is = 0;
  for (int i = 0; i < Km1+1; ++i)
  {
    if (i != h)
    {
      out(is) = 1.0 / sqrt(tau(i) * lambda[i]);
      is += 1;
    }
  }
  
  return out;
}


Eigen::VectorXd rmvnorm (Eigen::VectorXd& v, const int& j, const double& s_22, const Eigen::VectorXd& s_12, const Eigen::MatrixXd& inv_Omega_11, const Eigen::VectorXd& d_inv, const Eigen::VectorXd& temp_2, const int& pm1)
{
  // sample from multivariate normal distribution
  
  Eigen::MatrixXd Dev = s_22 * inv_Omega_11;
  Dev.diagonal() += d_inv;
  Eigen::VectorXd tmp = d_inv.asDiagonal() * temp_2 - s_12;
  
  const Eigen::MatrixXd L(Dev.llt().matrixL());
  
  tmp = L.triangularView<Eigen::Lower>().solve(tmp);
  v = L.transpose().triangularView<Eigen::Upper>().solve(v + tmp);
  
  return v;
}


// [[Rcpp::export]]

double pgamma3p1 (const double& x, const double& a, const double& b) 
{

  const double b2a = b / (2.0 * a);
  const double b2a2 = b2a * b2a;
  const double b2as = s2 * b2a;
  
  const double u = a * x - b2a;
  const double p1 = R::pnorm(-b2as, 0.0, 1.0, 1, 0);
  const double e1 = std::exp(-b2a2);
  
  const double tmp1 = 1.77245 * (R::pnorm(s2 * u, 0.0, 1.0, 1, 0) - p1) * b2a + 0.5 * (e1 - std::exp(-u * u));
  
  const double d1 = e1 - s2pi * b2as * (p1 - 1.0);
  const double out = (2.0 * tmp1) / d1;
  
  return out;
}


struct lambda_parallel_step : public RcppParallel::Worker {
  
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::
  // inputs
  
  const double mth, mth2;
  
  const Eigen::VectorXd v, temp_1, eta_12;
  
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::
  // output
  
  RcppParallel::RVector<double> lambda_12_par, eta_12_par;
  
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::
  // initialization
  
  lambda_parallel_step (const Eigen::VectorXd& v, const Eigen::VectorXd& temp_1, const Eigen::VectorXd& eta_12, const double& mth, const double& mth2, Rcpp::NumericVector& lambda_12_par, Rcpp::NumericVector& eta_12_par) 
    : v(v), temp_1(temp_1), eta_12(eta_12), mth(mth), mth2(mth2), lambda_12_par(lambda_12_par), eta_12_par(eta_12_par) {}
  
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::
  // main
  
  void operator()(std::size_t k_0, std::size_t k_1) 
  {
    std::size_t j;
    
    // :::::::::::::::::::::::::::::::::::::::::::::::::::::
    // parallel cycle
    
    for (j = k_0; j < k_1; j++) 
    {
      double tmp1 = v(j);
      double alpha_l = std::sqrt(0.5 * tmp1 * tmp1 / mth + 1.0 / eta_12(j));
      double beta_l = temp_1(j) * tmp1 / mth2;
      
      tmp1 = rg3p_c1(alpha_l, beta_l);
      tmp1 = 1.0 / (tmp1 * tmp1);
      
      lambda_12_par[j] = tmp1;
      eta_12_par[j] = 1.0 / R::rgamma(1.0, 1.0 / (1.0 + 1.0 / tmp1));
    }
  }
};


std::tuple<Eigen::VectorXd, Eigen::VectorXd> lambda_pstep (const Eigen::VectorXd& v, const Eigen::VectorXd& temp_1, const Eigen::VectorXd& eta_12, const double& mth, const double& mth2)
{
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::
  // output
  
  int pm1 = v.size();
  Rcpp::NumericVector lambda_12_par(pm1), eta_12_par(pm1);
  
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::
  // parallel step
  
  lambda_parallel_step Lambda_parallel_step(v, temp_1, eta_12, mth, mth2, lambda_12_par, eta_12_par);
  RcppParallel::parallelFor (0, pm1, Lambda_parallel_step);
  
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::
  // return output
  
  Eigen::Map<Eigen::VectorXd> l12(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(lambda_12_par));
  Eigen::Map<Eigen::VectorXd> e12(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(eta_12_par));
  
  return std::make_tuple(l12, e12);
}

// end file

