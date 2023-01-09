#include <RcppEigen.h>
#include <RcppNumerical.h>
#include <RcppParallel.h>

#include <Eigen/Core>

#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppNumerical)]]

#include "utils_InfHS.h"
#include "solve_quartic.h"


// :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// IG log-expectation
// :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

double mean_ig_logx (const double& a, const double& b)
{
  return std::log(b) - R::digamma(a);
}


// :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// Solve covariance matrix 
// :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Eigen::VectorXd inverse (Eigen::VectorXd& x)
{
  for (int j = 0; j < x.size(); ++j)
    x(j) = 1.0 / x(j);
  return x;
}


// Eigen::VectorXd repelem (const double& x, const int& n) 
// {
//   Eigen::VectorXd out(n);
//   for (int i = 0; i < n; i++) 
//     out(i) = x;
//   
//   return out;
// }


Eigen::MatrixXd solve_largep (const Eigen::MatrixXd& X, Eigen::VectorXd& dinv) 
{
  const int n = X.rows();
  const Eigen::VectorXd d = inverse(dinv);
  const Eigen::MatrixXd XD = X * d.asDiagonal();
  
  Eigen::MatrixXd D1 = XD * X.transpose();
  D1.diagonal() += repelem(1.0, n);
  
  Eigen::MatrixXd R = D1.llt().matrixL();
  D1 = R.triangularView<Eigen::Lower>().solve(XD);
  
  D1 = -D1.transpose() * D1;
  D1.diagonal() += d;
  
  return D1;
}


std::tuple<double, Eigen::MatrixXd> solve_largep_fastdet (const Eigen::MatrixXd& X, Eigen::VectorXd& dinv) 
{
  const int n = X.rows();
  const Eigen::VectorXd d = inverse(dinv);
  const Eigen::MatrixXd XD = X * d.asDiagonal();
  
  Eigen::MatrixXd D1 = XD * X.transpose();
  D1.diagonal() += repelem(1.0, n);
  
  // log-determinant
  double ldet = 0.0;
  for (int i = 0; i < d.size(); i++)
    ldet += std::log(d(i));
  ldet -= std::log(D1.determinant());
  
  Eigen::MatrixXd R = D1.llt().matrixL();
  D1 = R.triangularView<Eigen::Lower>().solve(XD);
  
  D1 = -D1.transpose() * D1;
  D1.diagonal() += d;
  
  return std::make_tuple(ldet, D1);
}


Eigen::MatrixXd solve_smallp (const Eigen::MatrixXd& X) 
{
  const int p = X.cols();
  Eigen::MatrixXd I(p, p); I.setIdentity();
  Eigen::MatrixXd R = X.llt().matrixL();
  Eigen::MatrixXd L = R.triangularView<Eigen::Lower>().solve(I);
  
  return L.transpose().triangularView<Eigen::Upper>() * L;
}



// :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// compute S * X^t * y, diag(S) and tr(XSX^t) in linear regression 
// :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

std::tuple<Eigen::VectorXd, Eigen::VectorXd, double, double> step_VB_det (Eigen::VectorXd& dinv, const Eigen::MatrixXd& X, const Eigen::VectorXd& y)
{
  const int n = X.rows();
  const Eigen::VectorXd d = inverse(dinv);
  const Eigen::MatrixXd DX = d.asDiagonal() * X.transpose();
  Eigen::MatrixXd XDX = X * DX;
  Eigen::MatrixXd Sn = XDX;
  Sn.diagonal() += repelem(1.0, n);

  Eigen::MatrixXd R = Sn.llt().matrixL();
  double ldet = 0.0;
  for (int i = 0; i < n; i++)
    ldet -= 2.0 * std::log(R(i, i));

  Sn = DX * solve_smallp(Sn);

  Eigen::VectorXd dSigma(d.size());
  for (int j = 0; j < d.size(); j++)
  {
    ldet += std::log(d(j));

    dSigma(j) = d(j);
    for (int i = 0; i < n; i++)
      dSigma(j) -= Sn(j, i) * DX(j, i);
  }

  Eigen::MatrixXd SX = DX - Sn * XDX;
  Eigen::VectorXd mu = SX * y;

  double xSx = 0.0;
  for (int i = 0; i < n; i++)
    xSx += X.row(i) * SX.col(i);

  return std::make_tuple(mu, dSigma, ldet, xSx);
}


std::tuple<Eigen::VectorXd, Eigen::VectorXd, double> step_VB_det0 (Eigen::VectorXd& dinv, const Eigen::MatrixXd& X, const Eigen::VectorXd& y)
{
  const int n = X.rows();
  const Eigen::VectorXd d = inverse(dinv);
  const Eigen::MatrixXd DX = d.asDiagonal() * X.transpose();
  Eigen::MatrixXd Sn = X * DX;
  Sn.diagonal() += repelem(1.0, n);

  double ldet = 0.0;
  Eigen::MatrixXd R = Sn.llt().matrixL();
  for (int i = 0; i < n; i++)
    ldet -= 2.0 * std::log(R(i, i));
  Sn = DX * solve_smallp(Sn);

  Eigen::VectorXd dSigma(d.size());
  for (int j = 0; j < d.size(); j++)
  {
    ldet += std::log(d(j));

    dSigma(j) = d(j);
    for (int i = 0; i < n; i++)
      dSigma(j) -= Sn(j, i) * DX(j, i);
  }

  Eigen::VectorXd mu = X.transpose() * y;
  mu = d.asDiagonal() * mu - Sn * (DX.transpose() * mu);

  return std::make_tuple(mu, dSigma, ldet);
}


// :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// compute S * X^t * y and diag(S) in probit regression 
// :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

std::tuple<Eigen::VectorXd, Eigen::VectorXd, double> step_VB_pr_det (Eigen::VectorXd& dinv, const Eigen::MatrixXd& X, const Eigen::VectorXd& ys) 
{
  const int n = X.rows();
  const Eigen::VectorXd d = inverse(dinv);
  const Eigen::MatrixXd DX = d.asDiagonal() * X.transpose();
  Eigen::MatrixXd Sn = X * DX;
  Sn.diagonal() += repelem(1.0, n);
  
  double ldet = 0.0;
  Eigen::MatrixXd R = Sn.llt().matrixL();
  for (int i = 0; i < n; i++)
    ldet -= 2.0 * std::log(R(i, i));  
  Sn = DX * solve_smallp(Sn);
  
  Eigen::VectorXd dSigma(d.size());
  for (int j = 0; j < d.size(); j++)
  {
    ldet += std::log(d(j));
    
    dSigma(j) = d(j);
    for (int i = 0; i < n; i++)
      dSigma(j) -= Sn(j, i) * DX(j, i);
  }
  
  Eigen::VectorXd mu = X.transpose() * ys;
  mu = d.asDiagonal() * mu - Sn * (DX.transpose() * mu);
  
  return std::make_tuple(mu, dSigma, ldet);
}


// :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// log-normalizing constant and moments of lambda 
// :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

class f0: public Numer::Func {
private:
  double lmax, d, a, b;
public:
  f0(double lmax, double d, double a, double b) : lmax(lmax), d(d), a(a), b(b) {}
  
  double operator()(const double& x) const {
    return exp(-std::log(x) - d / (x * x) - a * a * x * x + b * x - lmax);
  }
};


// [[Rcpp::export]]

double log_kf (const double& d, const double& a, const double& b)
{
  double a2 = a * a;
  Eigen::VectorXd coeff(5);
  coeff(0) = 2.0 * d;
  coeff(1) = 0.0;
  coeff(2) = -1.0;
  coeff(3) = b;
  coeff(4) = -2.0 * a2;
  
  double fmax = solve_quartic(coeff);
  double fmax2 = fmax * fmax;
  double lmax = -std::log(fmax) - d / (fmax2) - a2 * fmax2 + b * fmax;
  
  double err_est;
  int err_code;
  f0 f(lmax, d, a, b);
  double i0 = integrate(f, 0.0, std::numeric_limits<double>::infinity(), err_est, err_code);
  if (i0 < 1e-04)
    i0 = integrate(f, 0.0, fmax + 1000.0, err_est, err_code);
  if (i0 < 1e-04)
    i0 = integrate(f, 0.0, fmax + 500.0, err_est, err_code);
  if (i0 < 1e-04)
    i0 = integrate(f, 0.0, fmax + 250.0, err_est, err_code);
  if (i0 < 1e-04)
    i0 = integrate(f, 0.0, fmax + 100.0, err_est, err_code);
  
  return std::log(i0) + lmax;
}



class f1: public Numer::Func {
private:
  double lmax, d, a, b, lkf;
public:
  f1(double lmax, double d, double a, double b, double lkf) : lmax(lmax), d(d), a(a), b(b), lkf(lkf) {}
  
  double operator()(const double& x) const {
    return exp(-d / (x * x) - a * a * x * x + b * x - lkf - lmax);
  }
};


// [[Rcpp::export]]

double m_e1 (const double& d, const double& a, const double& b, const double& lkf)
{
  double a2 = a * a;
  Eigen::VectorXd coeff(5);
  coeff(0) = 2.0 * d;
  coeff(1) = 0.0;
  coeff(2) = 0.0;
  coeff(3) = b;
  coeff(4) = -2.0 * a2;
  
  double fmax = solve_quartic(coeff);
  double fmax2 = fmax * fmax;
  double lmax = -d / (fmax2) - a2 * fmax2 + b * fmax - lkf;
  
  double err_est;
  int err_code;
  f1 f(lmax, d, a, b, lkf);
  double i0 = integrate(f, 0.0, std::numeric_limits<double>::infinity(), err_est, err_code);
  if (i0 < 1e-04)
    i0 = integrate(f, 0.0, fmax + 1000.0, err_est, err_code);
  if (i0 < 1e-04)
    i0 = integrate(f, 0.0, fmax + 500.0, err_est, err_code);
  if (i0 < 1e-04)
    i0 = integrate(f, 0.0, fmax + 250.0, err_est, err_code);
  if (i0 < 1e-04)
    i0 = integrate(f, 0.0, fmax + 100.0, err_est, err_code);
  
  return std::exp(std::log(i0) + lmax);
}



class f2: public Numer::Func {
private:
  double lmax, d, a, b, lkf;
public:
  f2(double lmax, double d, double a, double b, double lkf) : lmax(lmax), d(d), a(a), b(b), lkf(lkf) {}
  
  double operator()(const double& x) const {
    return exp(std::log(x) - d / (x * x) - a * a * x * x + b * x - lkf - lmax);
  }
};

// [[Rcpp::export]]

double m_e2 (const double& d, const double& a, const double& b, const double& lkf)
{
  double a2 = a * a;
  Eigen::VectorXd coeff(5);
  coeff(0) = 2.0 * d;
  coeff(1) = 0.0;
  coeff(2) = 1.0;
  coeff(3) = b;
  coeff(4) = -2.0 * a2;
  
  double fmax = solve_quartic(coeff);
  double fmax2 = fmax * fmax;
  double lmax = std::log(fmax) - d / (fmax2) - a2 * fmax2 + b * fmax - lkf;
  
  double err_est;
  int err_code;
  f2 f(lmax, d, a, b, lkf);
  double i0 = integrate(f, 0.0, std::numeric_limits<double>::infinity(), err_est, err_code);
  if (i0 < 1e-04)
    i0 = integrate(f, 0.0, fmax + 1000.0, err_est, err_code);
  if (i0 < 1e-04)
    i0 = integrate(f, 0.0, fmax + 500.0, err_est, err_code);
  if (i0 < 1e-04)
    i0 = integrate(f, 0.0, fmax + 250.0, err_est, err_code);
  if (i0 < 1e-04)
    i0 = integrate(f, 0.0, fmax + 100.0, err_est, err_code);
  
  return std::exp(std::log(i0) + lmax);
}



class f3: public Numer::Func {
private:
  double lmax, d, a, b, lkf;
public:
  f3(double lmax, double d, double a, double b, double lkf) : lmax(lmax), d(d), a(a), b(b), lkf(lkf) {}
  
  double operator()(const double& x) const {
    return exp(-3.0 * std::log(x) - d / (x * x) - a * a * x * x + b * x - lkf - lmax);
  }
};

// [[Rcpp::export]]

double m_e3 (const double& d, const double& a, const double& b, const double& lkf)
{
  double a2 = a * a;
  Eigen::VectorXd coeff(5);
  coeff(0) = 2.0 * d;
  coeff(1) = 0.0;
  coeff(2) = -3.0;
  coeff(3) = b;
  coeff(4) = -2.0 * a2;
  
  double fmax = solve_quartic(coeff);
  double fmax2 = fmax * fmax;
  double lmax = -3.0 * std::log(fmax) - d / (fmax2) - a2 * fmax2 + b * fmax - lkf;
  
  double err_est;
  int err_code;
  f3 f(lmax, d, a, b, lkf);
  double i0 = integrate(f, 0.0, std::numeric_limits<double>::infinity(), err_est, err_code);
  if (i0 < 1e-04)
    i0 = integrate(f, 0.0, fmax + 1000.0, err_est, err_code);
  if (i0 < 1e-04)
    i0 = integrate(f, 0.0, fmax + 500.0, err_est, err_code);
  if (i0 < 1e-04)
    i0 = integrate(f, 0.0, fmax + 250.0, err_est, err_code);
  if (i0 < 1e-04)
    i0 = integrate(f, 0.0, fmax + 100.0, err_est, err_code);
  
  return std::exp(std::log(i0) + lmax);
}


// ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// PARALLEL iHS STEP

struct ihs_parallel_step0 : public RcppParallel::Worker {
  // function that return the updated inverse matrix after updating m columns
  // Q, R starting QR decomposition
  // k number of trials
  // X data matrix
  // l_det log-determinant
  // v1 slab standard deviation
  // gamma current model
  // ind variables in the current model
  // hyp hyper-parameters of sigma2
  // phi current phi
  // change variables to remove in the k models
  // q vector with number of variables to be modified for each proposal
  
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::
  // inputs
  
  const double s0, s02inv;
  
  const Eigen::VectorXd mu_beta, dSigma;
  
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::
  // output
  
  RcppParallel::RVector<double> a_lambda, b_lambda, c_lambda, d_phi, mu_phi, mu_lambda, lambda_phi, bl_sum, tmp_j;
  
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::
  // initialization
  
  ihs_parallel_step0 (const Eigen::VectorXd& mu_beta, const Eigen::VectorXd& dSigma, const double& s0, const double& s02inv, Rcpp::NumericVector& a_lambda, Rcpp::NumericVector& b_lambda, Rcpp::NumericVector& c_lambda, Rcpp::NumericVector& d_phi, Rcpp::NumericVector& mu_lambda, Rcpp::NumericVector& mu_phi, Rcpp::NumericVector& lambda_phi, Rcpp::NumericVector& bl_sum, Rcpp::NumericVector& tmp_j) 
    : mu_beta(mu_beta), dSigma(dSigma), s0(s0), s02inv(s02inv), a_lambda(a_lambda), b_lambda(b_lambda), c_lambda(c_lambda), d_phi(d_phi), mu_lambda(mu_lambda), mu_phi(mu_phi), lambda_phi(lambda_phi), bl_sum(bl_sum), tmp_j(tmp_j) {}
  
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::
  // main
  
  void operator()(std::size_t k_0, std::size_t k_1) 
  {
    std::size_t j;
    
    // :::::::::::::::::::::::::::::::::::::::::::::::::::::
    // parallel cycle
    
    for (j = k_0; j < k_1; j++) 
    {
      double mu_beta2 = std::pow(mu_beta(j+1), 2.0) + dSigma(j+1);
      
      a_lambda[j] = 0.5 * mu_beta2;
      b_lambda[j] = s02inv;
      c_lambda[j] = 0.0;
      
      double b_lambda_sq = std::sqrt(b_lambda[j]);
      double lkf = log_kf(a_lambda[j], b_lambda_sq, c_lambda[j]);
      double e1 = m_e1(a_lambda[j], b_lambda_sq, c_lambda[j], lkf);
      double e2 = m_e2(a_lambda[j], b_lambda_sq, c_lambda[j], lkf);
      double e3 = m_e3(a_lambda[j], b_lambda_sq, c_lambda[j], lkf);
      
      bl_sum[j] = mu_beta2 * e3;
      
      d_phi[j] = 0.5 + s02inv * e2;
      
      mu_lambda[j+1] = e3;
      mu_phi[j] = 1.0 / d_phi[j];
      lambda_phi[j] = mu_phi[j] * e1;
      
      tmp_j[j] = -std::log(0.5) + lkf + a_lambda[j] * e3 + b_lambda[j] * e2 - std::log(d_phi[j]);
    }
  }
};


std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> ihs_pstep0 (const Eigen::VectorXd& mu_beta, const Eigen::VectorXd& dSigma, const double& s0, const double& s02inv) 
{
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::
  // output
  
  int p = mu_beta.size();
  Rcpp::NumericVector a_lambda(p-1), b_lambda(p-1), c_lambda(p-1), d_phi(p-1), mu_lambda(p), mu_phi(p-1), lambda_phi(p-1), bl_sum(p-1), tmp_j(p-1);
  
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::
  // parallel step
  
  ihs_parallel_step0 Ihs_parallel_step0(mu_beta, dSigma, s0, s02inv, a_lambda, b_lambda, c_lambda, d_phi, mu_lambda, mu_phi, lambda_phi, bl_sum, tmp_j);
  RcppParallel::parallelFor (0, p-1, Ihs_parallel_step0);
  
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::
  // return output
  
  Eigen::Map<Eigen::VectorXd> al(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(a_lambda));
  Eigen::Map<Eigen::VectorXd> bl(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(b_lambda));
  Eigen::Map<Eigen::VectorXd> cl(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(c_lambda));
  Eigen::Map<Eigen::VectorXd> d(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(d_phi));
  Eigen::Map<Eigen::VectorXd> lam(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(mu_lambda));
  Eigen::Map<Eigen::VectorXd> phi(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(mu_phi));
  Eigen::Map<Eigen::VectorXd> lamphi(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(lambda_phi));
  Eigen::Map<Eigen::VectorXd> bsum(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(bl_sum));
  Eigen::Map<Eigen::VectorXd> tj(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(tmp_j));
  
  return std::make_tuple(al, bl, cl, d, lam, phi, lamphi, bsum, tj);
}


struct ihs_parallel_step : public RcppParallel::Worker {
  // function that return the updated inverse matrix after updating m columns
  // Q, R starting QR decomposition
  // k number of trials
  // X data matrix
  // l_det log-determinant
  // v1 slab standard deviation
  // gamma current model
  // ind variables in the current model
  // hyp hyper-parameters of sigma2
  // phi current phi
  // change variables to remove in the k models
  // q vector with number of variables to be modified for each proposal
  
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::
  // inputs
  
  const double mu_tau2, mu_sigma, mu_sigma2, s0, s02inv;
  
  const Eigen::VectorXd mu_beta, mu_gamma, dSigma, dphi, muphi;
  const Eigen::MatrixXd Z, Sigma_gamma;
  
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::
  // output
  
  RcppParallel::RVector<double> a_lambda, b_lambda, c_lambda, d_phi, mu_phi, mu_lambda, lambda_phi, bl_sum, tmp_j;
  
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::
  // initialization
  
  ihs_parallel_step (const Eigen::MatrixXd& Z, const Eigen::VectorXd& mu_beta, const Eigen::VectorXd& dSigma, const Eigen::VectorXd& mu_gamma, const Eigen::MatrixXd& Sigma_gamma, const Eigen::VectorXd& dphi, const Eigen::VectorXd& muphi, const double& mu_tau2, const double& mu_sigma, const double& mu_sigma2, const double& s0, const double& s02inv, Rcpp::NumericVector& a_lambda, Rcpp::NumericVector& b_lambda, Rcpp::NumericVector& c_lambda, Rcpp::NumericVector& d_phi, Rcpp::NumericVector& mu_lambda, Rcpp::NumericVector& mu_phi, Rcpp::NumericVector& lambda_phi, Rcpp::NumericVector& bl_sum, Rcpp::NumericVector& tmp_j) 
    : Z(Z), mu_beta(mu_beta), dSigma(dSigma), mu_gamma(mu_gamma), Sigma_gamma(Sigma_gamma), dphi(dphi), muphi(muphi), mu_tau2(mu_tau2), mu_sigma(mu_sigma), mu_sigma2(mu_sigma2), s0(s0), s02inv(s02inv), a_lambda(a_lambda), b_lambda(b_lambda), c_lambda(c_lambda), d_phi(d_phi), mu_lambda(mu_lambda), mu_phi(mu_phi), lambda_phi(lambda_phi), bl_sum(bl_sum), tmp_j(tmp_j) {}
  
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::
  // main
  
  void operator()(std::size_t k_0, std::size_t k_1) 
  {
    std::size_t j;
    
    // :::::::::::::::::::::::::::::::::::::::::::::::::::::
    // parallel cycle
    
    for (j = k_0; j < k_1; j++) 
    {
      Eigen::VectorXd zj = Z.row(j).transpose();
      
      double mu_beta2 = std::pow(mu_beta(j+1), 2.0) + dSigma(j+1) * mu_sigma;
      double mu_j = zj.transpose() * mu_gamma;
      
      a_lambda[j] = 0.5 * mu_beta2 * mu_tau2 * mu_sigma2;
      b_lambda[j] = s02inv * muphi(j);
      c_lambda[j] = mu_j * muphi(j) / s0;
      
      double b_lambda_sq = std::sqrt(b_lambda[j]);
      double lkf_p = std::log(1.0 - R::pnorm(0, mu_j, std::sqrt(2.0 * s0 * dphi(j)), true, false));
      double lkf = log_kf(a_lambda[j], b_lambda_sq, c_lambda[j]);
      double e1 = m_e1(a_lambda[j], b_lambda_sq, c_lambda[j], lkf);
      double e2 = m_e2(a_lambda[j], b_lambda_sq, c_lambda[j], lkf);
      double e3 = m_e3(a_lambda[j], b_lambda_sq, c_lambda[j], lkf);
      
      bl_sum[j] = mu_beta2 * e3;
      
      d_phi[j] = 0.5 + s02inv * (e2 - 2.0 * e1 * mu_j + mu_j * mu_j + s0 * (zj.transpose() * Sigma_gamma * zj).value());
      
      mu_lambda[j+1] = e3;
      mu_phi[j] = 1.0 / d_phi[j];
      lambda_phi[j] = mu_phi[j] * e1;
      
      tmp_j[j] = -lkf_p + lkf + a_lambda[j] * e3 + b_lambda[j] * e2 - c_lambda[j] * e1 - std::log(d_phi[j]);
    }
  }
};


std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> ihs_pstep (const Eigen::MatrixXd& Z, const Eigen::VectorXd& mu_beta, const Eigen::VectorXd& dSigma, const Eigen::VectorXd& mu_gamma, const Eigen::MatrixXd& Sigma_gamma, const Eigen::VectorXd& dphi, const Eigen::VectorXd& muphi, const double& mu_tau2, const double& mu_sigma, const double& mu_sigma2, const double& s0, const double& s02inv) 
{
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::
  // output
  
  int p = mu_beta.size();
  Rcpp::NumericVector a_lambda(p-1), b_lambda(p-1), c_lambda(p-1), d_phi(p-1), mu_lambda(p), mu_phi(p-1), lambda_phi(p-1), bl_sum(p-1), tmp_j(p-1);
  
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::
  // parallel step
  
  ihs_parallel_step Ihs_parallel_step(Z, mu_beta, dSigma, mu_gamma, Sigma_gamma, dphi, muphi, mu_tau2, mu_sigma, mu_sigma2, s0, s02inv, a_lambda, b_lambda, c_lambda, d_phi, mu_lambda, mu_phi, lambda_phi, bl_sum, tmp_j);
  RcppParallel::parallelFor (0, p-1, Ihs_parallel_step);
  
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::
  // return output
  
  Eigen::Map<Eigen::VectorXd> al(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(a_lambda));
  Eigen::Map<Eigen::VectorXd> bl(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(b_lambda));
  Eigen::Map<Eigen::VectorXd> cl(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(c_lambda));
  Eigen::Map<Eigen::VectorXd> d(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(d_phi));
  Eigen::Map<Eigen::VectorXd> lam(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(mu_lambda));
  Eigen::Map<Eigen::VectorXd> phi(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(mu_phi));
  Eigen::Map<Eigen::VectorXd> lamphi(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(lambda_phi));
  Eigen::Map<Eigen::VectorXd> bsum(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(bl_sum));
  Eigen::Map<Eigen::VectorXd> tj(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(tmp_j));
  
  return std::make_tuple(al, bl, cl, d, lam, phi, lamphi, bsum, tj);
}



struct ihs_pr_parallel_step : public RcppParallel::Worker {
  // function that return the updated inverse matrix after updating m columns
  // Q, R starting QR decomposition
  // k number of trials
  // X data matrix
  // l_det log-determinant
  // v1 slab standard deviation
  // gamma current model
  // ind variables in the current model
  // hyp hyper-parameters of sigma2
  // phi current phi
  // change variables to remove in the k models
  // q vector with number of variables to be modified for each proposal
  
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::
  // inputs
  
  const double mu_tau2, s0, s02inv;
  
  const Eigen::VectorXd mu_beta, mu_gamma, dSigma, dphi, muphi;
  const Eigen::MatrixXd Z, Sigma_gamma;
  
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::
  // output
  
  RcppParallel::RVector<double> a_lambda, b_lambda, c_lambda, d_phi, mu_phi, mu_lambda, lambda_phi, bl_sum, tmp_j;
  
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::
  // initialization
  
  ihs_pr_parallel_step (const Eigen::MatrixXd& Z, const Eigen::VectorXd& mu_beta, const Eigen::VectorXd& dSigma, const Eigen::VectorXd& mu_gamma, const Eigen::MatrixXd& Sigma_gamma, const Eigen::VectorXd& dphi, const Eigen::VectorXd& muphi, const double& mu_tau2, const double& s0, const double& s02inv, Rcpp::NumericVector& a_lambda, Rcpp::NumericVector& b_lambda, Rcpp::NumericVector& c_lambda, Rcpp::NumericVector& d_phi, Rcpp::NumericVector& mu_lambda, Rcpp::NumericVector& mu_phi, Rcpp::NumericVector& lambda_phi, Rcpp::NumericVector& bl_sum, Rcpp::NumericVector& tmp_j) 
    : Z(Z), mu_beta(mu_beta), dSigma(dSigma), mu_gamma(mu_gamma), Sigma_gamma(Sigma_gamma), dphi(dphi), muphi(muphi), mu_tau2(mu_tau2), s0(s0), s02inv(s02inv), a_lambda(a_lambda), b_lambda(b_lambda), c_lambda(c_lambda), d_phi(d_phi), mu_lambda(mu_lambda), mu_phi(mu_phi), lambda_phi(lambda_phi), bl_sum(bl_sum), tmp_j(tmp_j) {}
  
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::
  // main
  
  void operator()(std::size_t k_0, std::size_t k_1) 
  {
    std::size_t j;
    
    // :::::::::::::::::::::::::::::::::::::::::::::::::::::
    // parallel cycle
    
    for (j = k_0; j < k_1; j++) 
    {
      Eigen::VectorXd zj = Z.row(j).transpose();
      
      double mu_beta2 = std::pow(mu_beta(j+1), 2.0) + dSigma(j+1);
      double mu_j = zj.transpose() * mu_gamma;
      
      a_lambda[j] = 0.5 * mu_beta2 * mu_tau2;
      b_lambda[j] = s02inv * muphi(j);
      c_lambda[j] = mu_j * muphi(j) / s0;
      
      double b_lambda_sq = std::sqrt(b_lambda[j]);
      double lkf_p = std::log(1.0 - R::pnorm(0, mu_j, std::sqrt(2.0 * s0 * dphi(j)), true, false));
      double lkf = log_kf(a_lambda[j], b_lambda_sq, c_lambda[j]);
      double e1 = m_e1(a_lambda[j], b_lambda_sq, c_lambda[j], lkf);
      double e2 = m_e2(a_lambda[j], b_lambda_sq, c_lambda[j], lkf);
      double e3 = m_e3(a_lambda[j], b_lambda_sq, c_lambda[j], lkf);
      
      bl_sum[j] = mu_beta2 * e3;
      
      d_phi[j] = 0.5 + s02inv * (e2 - 2.0 * e1 * mu_j + mu_j * mu_j + s0 * (zj.transpose() * Sigma_gamma * zj).value());
      
      mu_lambda[j+1] = e3;
      mu_phi[j] = 1.0 / d_phi[j];
      lambda_phi[j] = mu_phi[j] * e1;
      
      tmp_j[j] = -lkf_p + lkf + a_lambda[j] * e3 + b_lambda[j] * e2 - c_lambda[j] * e1 - std::log(d_phi[j]);
    }
  }
};


std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> ihs_pstep_pr (const Eigen::MatrixXd& Z, const Eigen::VectorXd& mu_beta, const Eigen::VectorXd& dSigma, const Eigen::VectorXd& mu_gamma, const Eigen::MatrixXd& Sigma_gamma, const Eigen::VectorXd& dphi, const Eigen::VectorXd& muphi, const double& mu_tau2, const double& s0, const double& s02inv) 
{
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::
  // output
  
  int p = mu_beta.size();
  Rcpp::NumericVector a_lambda(p-1), b_lambda(p-1), c_lambda(p-1), d_phi(p-1), mu_lambda(p), mu_phi(p-1), lambda_phi(p-1), bl_sum(p-1), tmp_j(p-1);
  
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::
  // parallel step
  
  ihs_pr_parallel_step Ihs_pr_parallel_step(Z, mu_beta, dSigma, mu_gamma, Sigma_gamma, dphi, muphi, mu_tau2, s0, s02inv, a_lambda, b_lambda, c_lambda, d_phi, mu_lambda, mu_phi, lambda_phi, bl_sum, tmp_j);
  RcppParallel::parallelFor (0, p-1, Ihs_pr_parallel_step);
  
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::
  // return output
  
  Eigen::Map<Eigen::VectorXd> al(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(a_lambda));
  Eigen::Map<Eigen::VectorXd> bl(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(b_lambda));
  Eigen::Map<Eigen::VectorXd> cl(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(c_lambda));
  Eigen::Map<Eigen::VectorXd> d(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(d_phi));
  Eigen::Map<Eigen::VectorXd> lam(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(mu_lambda));
  Eigen::Map<Eigen::VectorXd> phi(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(mu_phi));
  Eigen::Map<Eigen::VectorXd> lamphi(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(lambda_phi));
  Eigen::Map<Eigen::VectorXd> bsum(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(bl_sum));
  Eigen::Map<Eigen::VectorXd> tj(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(tmp_j));
  
  return std::make_tuple(al, bl, cl, d, lam, phi, lamphi, bsum, tj);
}
