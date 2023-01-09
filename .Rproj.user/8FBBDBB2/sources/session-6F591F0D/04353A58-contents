#include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>     

#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]] 

#include "boost/multi_array.hpp"

#include "utils_rgamma3p.h"
#include "utils_mGHS.h"
#include "utils_rIW.h"

typedef boost::multi_array<double, 1> vector_type;
typedef boost::multi_array<double, 2> matrix_type;
typedef boost::multi_array<double, 3> array_type;


// [[Rcpp::export]]

Rcpp::List mGHS (const int& B, const int& bn, const Eigen::VectorXd& n, const Rcpp::List& S, const int& p, const Eigen::VectorXd hyp_ta, const int& ping = 2500) 
{
  // B number of iteration
  // bn burn-in
  // n vector K-dim
  // S array p x p x K
  
  typedef array_type::index_range range;
  
  // ::::::::::::::::::::::::::::::::::::::::::::
  // initialization
  
  // const int p = S.size(), pm1 = p-1;
  const int pm1 = p-1;
  const int K = n.size(), Km1 = K-1;
  
  int h = 0, hs = 0, b = 0, j = 0, js = 0, c = 0, cs = 0;
  
  // working objects
  const double IW_df = 0.5 * p * pm1;
  double gamma = 0.0;
  double alpha_l = 0.0, beta_l = 0.0;
  double alpha_t = 0.0, beta_t = 0.0;
  double s_22 = 0.0;
  double mth = 0.0, mth2 = 0.0;
  double tmp1 = 0.0, tmp2 = 0.0, sv = 0.0;
  double a_MH = 0.0;
  
  Eigen::VectorXd delta_12_inv(Km1); delta_12_inv.setZero();
  Eigen::VectorXd diag_V(K); diag_V.setZero();
  Eigen::VectorXd mu(K); mu.setOnes();
  Eigen::VectorXd e_ij(K); e_ij.setZero();
  Eigen::VectorXd v(pm1); v.setZero();
  Eigen::VectorXd temp_1(pm1); temp_1.setZero();
  Eigen::VectorXd temp_2(pm1); temp_2.setZero();
  Eigen::VectorXd temp_3(pm1); temp_3.setZero();
  Eigen::VectorXd s_12(pm1); s_12.setZero();
  Eigen::VectorXd sigma_12(pm1); sigma_12.setZero();
  // Eigen::VectorXd eta_12(pm1); eta_12.setZero();
  Eigen::VectorXd delta_12(pm1); delta_12.setZero();
  
  Eigen::MatrixXd Rinv(Km1, Km1); Rinv.setZero();
  Eigen::MatrixXd rtR(Km1, K); rtR.setZero();
  Eigen::MatrixXd IW_scale(K, K); IW_scale.setZero();
  Eigen::MatrixXd Psi(K, K); Psi.setZero();
  Eigen::MatrixXd R_s(K, K); R_s.setZero();
  Eigen::MatrixXd inv_Omega_11(pm1, pm1); inv_Omega_11.setZero();
  
  // save values
  Eigen::VectorXd tau(K); tau.setOnes();
  Eigen::VectorXd zeta(K); zeta.setOnes();
  
  array_type Omega(boost::extents[p][p][K]);
  array_type Sigma(boost::extents[p][p][K]);
  for (h = 0; h < K; h++)
  {
    for (j = 0; j < p; j++)
    {
      Omega[j][j][h] = 1.0;
      Sigma[j][j][h] = 1.0;
    }
  }
  
  array_type Lambda(boost::extents[p][p][K]);
  std::fill(Lambda.origin(), Lambda.origin() + K*p*p, 1.0);
  
  array_type Eta(boost::extents[p][p][K]);
  std::fill(Eta.origin(), Eta.origin() + K*p*p, 1.0);
  
  Eigen::MatrixXd R(K, K); R.setIdentity();
  array_type S1(boost::extents[p][p][K]);
  for (h = 0; h < K; h++)
    S1[boost::indices[range()][range()][h]] = fill_S(S[h], p);
  
  // output
  
  Eigen::VectorXd tau_out(K); tau_out.setZero();
  Eigen::MatrixXd Omega_jj_out(p, K); Omega_jj_out.setZero();
  Eigen::MatrixXd Omega_ij_out((int)IW_df, K); Omega_ij_out.setZero();
  Eigen::MatrixXd lambda_out((int)IW_df, K); lambda_out.setZero();
  Eigen::MatrixXd lambda_prob((int)IW_df, K); lambda_prob.setZero();
  Eigen::MatrixXd eta_out((int)IW_df, K); eta_out.setZero();
  Eigen::MatrixXd R_out(K, K); R_out.setZero();
  
  // Eigen::MatrixXd pi_chain(B-bn, (int)IW_df * K); pi_chain.setZero();
  
  double ta_s = 0.0, l2 = 0.0;
  
  Eigen::VectorXd l1(K); l1.setZero();
  Eigen::VectorXd qta2((int)IW_df); qta2.setZero();
  Eigen::VectorXd ta(K);

  for (h = 0; h < K; h++) 
    ta(h) = 0.5;
  
  Eigen::MatrixXd Z((int)IW_df, K); Z.setZero();
  Eigen::MatrixXd qta1((int)IW_df, K); qta1.setZero();
  
  for (h = 0; h < K; h++)
  {
    l1(h) = IW_df * std::log(0.5);
    for (j = 0; j < IW_df; j++)
    {
      qta1(j, h) = 0.5;
      Z(j, h) = R::rbinom(1, 0.5);
    }
  } 
  
  Eigen::VectorXd ta_out(K); ta_out.setZero();
  Eigen::MatrixXd l1_out(B, K); l1_out.setZero();
  
  // ::::::::::::::::::::::::::::::::::::::::::::
  // start algorithm
  // ::::::::::::::::::::::::::::::::::::::::::::
  
  for (b = 0; b < B; b++)
  {
    for (h = 0; h < K; h++)
    {
      for (j = 0; j < p; j++)
      {
        // ::::::::::::::::::::::::::::::::::::::::::::
        // sample gamma and v
        
        s_22 = S1[j][j][h];
        gamma = R::rgamma(0.5 * n(h) + 1.0, 2.0 / s_22);
        
        cs = 0;
        for (js = 0; js < p; js++)
        {
          if (js != j) 
          {
            c = 0;
            for (hs = 0; hs < K; hs++)
            {
              if (hs != h)
              {
                delta_12_inv(c) = Omega[js][j][hs] / std::sqrt(tau(hs) * Lambda[js][j][hs]);
                c += 1;
              }
            }
            
            tmp1 = tau(h) * Lambda[js][j][h];
            delta_12(cs) = tmp1;
            
            v(cs) = R::rnorm(0, 1);
            temp_1(cs) = (rtR.col(h).transpose() * delta_12_inv).value();
            temp_2(cs) = std::sqrt(tmp1) * temp_1(cs);
            
            s_12(cs) = S1[js][j][h];
            sigma_12(cs) = Sigma[js][j][h];
            
            cs += 1;
          }
        }
        
        delta_12 = inverse(delta_12, pm1) / mu(h);
        inv_Omega_11 = delete_rc_from_boost1(Sigma[boost::indices[range(0, p)][range(0, p)][h]], j, pm1);
        inv_Omega_11 -= sigma_12 * (sigma_12.transpose() / Sigma[j][j][h]);
        
        v = rmvnorm(v, j, s_22, s_12, inv_Omega_11, delta_12, temp_2, pm1);
        
        // ::::::::::::::::::::::::::::::::::::::::::::
        // sample lambda and eta
        
        mth = mu(h) * tau(h);
        mth2 = std::sqrt(tau(h)) * mu(h);
        
        temp_2 = inv_Omega_11 * v;
        temp_3 = temp_2 / gamma;
        
        c = 0;
        for (js = 0; js < p; js++)
        {
          if (js != j) 
          {
            tmp1 = v(c);
            Omega[js][j][h] = tmp1;
            Omega[j][js][h] = tmp1;
            
            alpha_l = std::sqrt(0.5 * tmp1 * tmp1 / mth + 1.0 / Eta[js][j][h]);
            beta_l = temp_1(c) * tmp1 / mth2;
            tmp1 = rg3p_c1(alpha_l, beta_l);
            
            tmp1 = 1.0 / (tmp1 * tmp1);
            Lambda[js][j][h] = tmp1;  
            Lambda[j][js][h] = tmp1;  
            
            tmp1 = 1.0 / R::rgamma(1, 1.0 / (1.0 + 1.0 / tmp1));
            Eta[js][j][h] = tmp1;
            Eta[j][js][h] = tmp1;
            
            tmp1 = -temp_3(c);
            Sigma[js][j][h] = tmp1;
            Sigma[j][js][h] = tmp1;
            
            c += 1;
          }
        }
        
        inv_Omega_11 += temp_2 * temp_3.transpose();
        
        Omega[j][j][h] = gamma + v.transpose() * temp_2;
        Sigma[j][j][h] = 1.0 / gamma;
        
        Sigma[boost::indices[range(0, p)][range(0, p)][h]] = get_Sigma(Sigma[boost::indices[range(0, p)][range(0, p)][h]], inv_Omega_11, j, pm1, p);
      }
    }
    
    // ::::::::::::::::::::::::::::::::::::::::::::
    // sample tau and zeta
    
    diag_V.setZero();
    
    for (h = 0; h < K; h++)
    {
      alpha_t = 0.0;
      beta_t = 0.0;
      
      // sv = std::pow(Omega[0][0][h], 2);
      sv = 0.0;
      
      for (j = 1; j < p; j++)
      {
        // sv += std::pow(Omega[j][j][h], 2);
        
        for (js = 0; js < j; js++)
        {
          tmp2 = Omega[js][j][h];
          sv += tmp2 * tmp2;   
          
          tmp1 = tmp2 / std::sqrt(Lambda[js][j][h]);
          alpha_t += tmp1 * tmp1;
          
          c = 0;
          for (hs = 0; hs < K; hs++)
          {
            if (hs != h)
            {
              delta_12_inv(c) = Omega[js][j][hs] / std::sqrt(tau(hs) * Lambda[js][j][hs]);
              c += 1;
            }
          }
          
          beta_t += tmp1 * (rtR.col(h).transpose() * delta_12_inv).value();
        }
      }
      
      diag_V(h) = 1.0 / std::sqrt(sv);
      alpha_t = std::sqrt(0.5 / mu(h) * alpha_t + 1.0 / zeta(h));
      
      tau(h) = 1.0 / std::pow(rg3p_approx(alpha_t, beta_t / mu(h), (int)IW_df), 2);
      // tau(h) = 1.0 / std::pow(rg3p(alpha_t, beta_t / mu(h), (int)IW_df), 2);
      zeta(h) = 1.0 / R::rgamma(1.0, 1.0 / (1.0 + 1.0 / tau(h)));
    }
    
    // ::::::::::::::::::::::::::::::::::::::::::::
    // sample R
    
    IW_scale.setZero();
    for (j = 1; j < p; j++)
    {
      for (js = 0; js < j; js++)
      {
        for (h = 0; h < K; h++)
          e_ij(h) = diag_V(h) * Omega[js][j][h] / std::sqrt(tau(h) * Lambda[js][j][h]);
        IW_scale += e_ij * e_ij.transpose();
      }
    }
    
    Psi = rIW(IW_df, IW_scale, K);
    e_ij = Psi.diagonal();
    e_ij = sinverse(e_ij, K);
    R_s = e_ij.asDiagonal() * Psi * e_ij.asDiagonal();
    
    a_MH = exp((K+1) * 0.5 * (std::log(R_s.determinant()) - std::log(R.determinant())));
    
    if (R::runif(0, 1) < a_MH)
    {
      R = R_s;
      
      for (h = 0; h < K; h++)
      {
        Rinv = delete_rc_from_eigen(R, h, Km1);
        Rinv = Rinv.inverse();
        rtR.col(h) = delete_elem_from_eigen_vec(R.col(h), h, Km1);
        rtR.col(h) = rtR.col(h).transpose() * Rinv;
        
        mu(h) = 1.0 - rtR.col(h).transpose() * delete_elem_from_eigen_vec(R.col(h), h, Km1);
      }
    }
    
    for (h = 0; h < K; h++)
    {
      ta_s = R::rbeta(hyp_ta(0), hyp_ta(1));
      l2 = 0.0;

      mth = mu(h) * tau(h);
      mth2 = std::sqrt(tau(h)) * mu(h);

      cs = 0;
      for (j = 0; j < p-1; j++)
      {
        for (js = j+1; js < p; js++)
        {
          c = 0;
          for (hs = 0; hs < K; hs++)
          {
            if (hs != h)
            {
              delta_12_inv(c) = Omega[js][j][hs] / std::sqrt(tau(hs) * Lambda[js][j][hs]);
              c += 1;
            }
          }

          tmp1 = Omega[js][j][h];
          alpha_l = std::sqrt(0.5 * tmp1 * tmp1 / mth + 1.0 / Eta[js][j][h]);
          beta_l = (rtR.col(h).transpose() * delta_12_inv).value() * tmp1 / mth2;

          qta2(cs) = 1 - pgamma3p1(std::sqrt((1.0 - ta_s) / ta_s), alpha_l, beta_l);
          if (Z(cs, h) == 1)
          {
            l2 += std::log(qta2(cs));
          } else 
          {
            l2 += std::log(1.0 - qta2(cs));
          }
          cs += 1;
          
          a_MH = l2 - l1(h);
          if (std::log(R::runif(0, 1)) < a_MH)
          {
            ta(h) = ta_s;
            l1(h) = l2;
            qta1.col(h) = qta2;
          }
        }
      }

      for (j = 0; j < IW_df; j++) 
      {
        Z(j, h) = R::rbinom(1, qta1(j, h));
      }
    }
    
    l1_out.row(b) = l1;
    
    if (b >= bn)
    {
      tau_out += tau;
      R_out += R;
      ta_out += ta;
      
      cs = 0;
      for (h = 0; h < K; h++)
      {
        c = 0;
        for (j = 0; j < p; j++)
        {
          Omega_jj_out(j, h) += Omega[j][j][h];
          for (js = j+1; js < p; js++)
          {
            Omega_ij_out(c, h) += Omega[js][j][h];
            lambda_out(c, h) += Lambda[js][j][h];
            lambda_prob(c, h) += 1.0 / (1.0 + Lambda[js][j][h]);
            eta_out(c, h) += Eta[js][j][h];
            c += 1;
            cs += 1;
          }
        }
      }
    }
    
    if ((ping != 0) && (b % ping == 0))
      Rcpp::Rcout << "\n iter: " << b;
  }
  
  // ::::::::::::::::::::::::::::::::::::::::::::
  // output
  // ::::::::::::::::::::::::::::::::::::::::::::
  
  return Rcpp::List::create(Rcpp::Named("Omega_ij") = Omega_ij_out / (B - bn), Rcpp::Named("Omega_jj") = Omega_jj_out / (B - bn), Rcpp::Named("ta") = ta_out / (B - bn), Rcpp::Named("l1_chain") = l1_out, Rcpp::Named("lambda") = lambda_out / (B - bn), Rcpp::Named("lambda_prob") = lambda_prob / (B - bn), Rcpp::Named("eta") = eta_out / (B - bn), Rcpp::Named("R") = R_out / (B - bn), Rcpp::Named("tau") = tau_out / (B - bn));
}


// [[Rcpp::export]]

Rcpp::List mGHS2 (const int& B, const int& bn, const Eigen::VectorXd& n, const Rcpp::List& S, const int& p, const Eigen::VectorXd hyp_ta, const int& ping = 2500) 
{
  // B number of iteration
  // bn burn-in
  // n vector K-dim
  // S array p x p x K
  
  typedef array_type::index_range range;
  
  // ::::::::::::::::::::::::::::::::::::::::::::
  // initialization
  
  // const int p = S.size(), pm1 = p-1;
  const int pm1 = p-1;
  const int K = n.size(), Km1 = K-1;
  
  int h = 0, hs = 0, b = 0, j = 0, js = 0, c = 0, cs = 0;
  
  // working objects
  const double IW_df = 0.5 * p * pm1;
  double gamma = 0.0;
  double alpha_l = 0.0, beta_l = 0.0;
  double alpha_t = 0.0, beta_t = 0.0;
  double s_22 = 0.0;
  double mth = 0.0, mth2 = 0.0;
  double tmp1 = 0.0, tmp2 = 0.0, sv = 0.0;
  double a_MH = 0.0;
  
  Eigen::VectorXd delta_12_inv(Km1); delta_12_inv.setZero();
  Eigen::VectorXd diag_V(K); diag_V.setZero();
  Eigen::VectorXd mu(K); mu.setOnes();
  Eigen::VectorXd e_ij(K); e_ij.setZero();
  Eigen::VectorXd v(pm1); v.setZero();
  Eigen::VectorXd temp_1(pm1); temp_1.setZero();
  Eigen::VectorXd temp_2(pm1); temp_2.setZero();
  Eigen::VectorXd temp_3(pm1); temp_3.setZero();
  Eigen::VectorXd s_12(pm1); s_12.setZero();
  Eigen::VectorXd sigma_12(pm1); sigma_12.setZero();
  Eigen::VectorXd eta_12(pm1); eta_12.setZero();
  Eigen::VectorXd delta_12(pm1); delta_12.setZero();
  Eigen::VectorXd lambda_12_par(pm1); lambda_12_par.setZero();
  Eigen::VectorXd eta_12_par(pm1); eta_12_par.setZero();
  
  Eigen::MatrixXd Rinv(Km1, Km1); Rinv.setZero();
  Eigen::MatrixXd rtR(Km1, K); rtR.setZero();
  Eigen::MatrixXd IW_scale(K, K); IW_scale.setZero();
  Eigen::MatrixXd Psi(K, K); Psi.setZero();
  Eigen::MatrixXd R_s(K, K); R_s.setZero();
  Eigen::MatrixXd inv_Omega_11(pm1, pm1); inv_Omega_11.setZero();
  
  // save values
  Eigen::VectorXd tau(K); tau.setOnes();
  Eigen::VectorXd zeta(K); zeta.setOnes();
  
  array_type Omega(boost::extents[p][p][K]);
  array_type Sigma(boost::extents[p][p][K]);
  for (h = 0; h < K; h++)
  {
    for (j = 0; j < p; j++)
    {
      Omega[j][j][h] = 1.0;
      Sigma[j][j][h] = 1.0;
    }
  }
  
  array_type Lambda(boost::extents[p][p][K]);
  std::fill(Lambda.origin(), Lambda.origin() + K*p*p, 1.0);
  
  array_type Eta(boost::extents[p][p][K]);
  std::fill(Eta.origin(), Eta.origin() + K*p*p, 1.0);
  
  Eigen::MatrixXd R(K, K); R.setIdentity();
  array_type S1(boost::extents[p][p][K]);
  for (h = 0; h < K; h++)
    S1[boost::indices[range()][range()][h]] = fill_S(S[h], p);
  
  // output
  
  Eigen::VectorXd tau_out(K); tau_out.setZero();
  Eigen::MatrixXd Omega_jj_out(p, K); Omega_jj_out.setZero();
  Eigen::MatrixXd Omega_ij_out((int)IW_df, K); Omega_ij_out.setZero();
  Eigen::MatrixXd lambda_out((int)IW_df, K); lambda_out.setZero();
  Eigen::MatrixXd lambda_prob((int)IW_df, K); lambda_prob.setZero();
  Eigen::MatrixXd eta_out((int)IW_df, K); eta_out.setZero();
  Eigen::MatrixXd R_out(K, K); R_out.setZero();
  
  // Eigen::MatrixXd pi_chain(B-bn, (int)IW_df * K); pi_chain.setZero();
  
  double ta_s = 0.0, l2 = 0.0;
  
  Eigen::VectorXd l1(K); l1.setZero();
  Eigen::VectorXd qta2((int)IW_df); qta2.setZero();
  Eigen::VectorXd ta(K);
  
  for (h = 0; h < K; h++) 
    ta(h) = 0.5;
  
  Eigen::MatrixXd Z((int)IW_df, K); Z.setZero();
  Eigen::MatrixXd qta1((int)IW_df, K); qta1.setZero();
  
  for (h = 0; h < K; h++)
  {
    l1(h) = IW_df * std::log(0.5);
    for (j = 0; j < IW_df; j++)
    {
      qta1(j, h) = 0.5;
      Z(j, h) = R::rbinom(1, 0.5);
    }
  } 
  
  Eigen::VectorXd ta_out(K); ta_out.setZero();
  Eigen::MatrixXd l1_out(B, K); l1_out.setZero();
  
  // ::::::::::::::::::::::::::::::::::::::::::::
  // start algorithm
  // ::::::::::::::::::::::::::::::::::::::::::::
  
  for (b = 0; b < B; b++)
  {
    for (h = 0; h < K; h++)
    {
      for (j = 0; j < p; j++)
      {
        // ::::::::::::::::::::::::::::::::::::::::::::
        // sample gamma and v
        
        s_22 = S1[j][j][h];
        gamma = R::rgamma(0.5 * n(h) + 1.0, 2.0 / s_22);
        
        cs = 0;
        for (js = 0; js < p; js++)
        {
          if (js != j) 
          {
            c = 0;
            for (hs = 0; hs < K; hs++)
            {
              if (hs != h)
              {
                delta_12_inv(c) = Omega[js][j][hs] / std::sqrt(tau(hs) * Lambda[js][j][hs]);
                c += 1;
              }
            }
            
            tmp1 = tau(h) * Lambda[js][j][h];
            delta_12(cs) = tmp1;
            
            v(cs) = R::rnorm(0, 1);
            temp_1(cs) = (rtR.col(h).transpose() * delta_12_inv).value();
            temp_2(cs) = std::sqrt(tmp1) * temp_1(cs);
            
            s_12(cs) = S1[js][j][h];
            sigma_12(cs) = Sigma[js][j][h];
            
            cs += 1;
          }
        }
        
        delta_12 = inverse(delta_12, pm1) / mu(h);
        inv_Omega_11 = delete_rc_from_boost1(Sigma[boost::indices[range(0, p)][range(0, p)][h]], j, pm1);
        inv_Omega_11 -= sigma_12 * (sigma_12.transpose() / Sigma[j][j][h]);
        
        v = rmvnorm(v, j, s_22, s_12, inv_Omega_11, delta_12, temp_2, pm1);
        
        // ::::::::::::::::::::::::::::::::::::::::::::
        // sample lambda and eta
        
        mth = mu(h) * tau(h);
        mth2 = std::sqrt(tau(h)) * mu(h);
        
        temp_2 = inv_Omega_11 * v;
        temp_3 = temp_2 / gamma;
        
        c = 0;
        for (js = 0; js < p; js++)
        {
          if (js != j) 
          {
            eta_12(c) = Eta[js][j][h];
            
            tmp1 = v(c);
            Omega[js][j][h] = tmp1;
            Omega[j][js][h] = tmp1;
            
            tmp1 = -temp_3(c);
            Sigma[js][j][h] = tmp1;
            Sigma[j][js][h] = tmp1;
            
            c += 1;
          }
        }
        
        std::tie(lambda_12_par, eta_12_par) = lambda_pstep(v, temp_1, eta_12, mth, mth2);
        
        c = 0;
        for (js = 0; js < p; js++)
        {
          if (js != j) 
          {
            tmp1 = lambda_12_par(c);
            // tmp1 = 1.0 / (tmp1 * tmp1);
            Lambda[js][j][h] = tmp1;  
            Lambda[j][js][h] = tmp1;  
            
            tmp1 = eta_12_par(c);
            Eta[js][j][h] = tmp1;
            Eta[j][js][h] = tmp1;
            
            c += 1;
          }
        }
        
        inv_Omega_11 += temp_2 * temp_3.transpose();
        
        Omega[j][j][h] = gamma + v.transpose() * temp_2;
        Sigma[j][j][h] = 1.0 / gamma;
        
        Sigma[boost::indices[range(0, p)][range(0, p)][h]] = get_Sigma(Sigma[boost::indices[range(0, p)][range(0, p)][h]], inv_Omega_11, j, pm1, p);
      }
    }
    
    // ::::::::::::::::::::::::::::::::::::::::::::
    // sample tau and zeta
    
    diag_V.setZero();
    
    for (h = 0; h < K; h++)
    {
      alpha_t = 0.0;
      beta_t = 0.0;
      
      // sv = std::pow(Omega[0][0][h], 2);
      sv = 0.0;
      
      for (j = 1; j < p; j++)
      {
        // sv += std::pow(Omega[j][j][h], 2);
        
        for (js = 0; js < j; js++)
        {
          tmp2 = Omega[js][j][h];
          sv += tmp2 * tmp2;   
          
          tmp1 = tmp2 / std::sqrt(Lambda[js][j][h]);
          alpha_t += tmp1 * tmp1;
          
          c = 0;
          for (hs = 0; hs < K; hs++)
          {
            if (hs != h)
            {
              delta_12_inv(c) = Omega[js][j][hs] / std::sqrt(tau(hs) * Lambda[js][j][hs]);
              c += 1;
            }
          }
          
          beta_t += tmp1 * (rtR.col(h).transpose() * delta_12_inv).value();
        }
      }
      
      diag_V(h) = 1.0 / std::sqrt(sv);
      alpha_t = std::sqrt(0.5 / mu(h) * alpha_t + 1.0 / zeta(h));
      
      tau(h) = 1.0 / std::pow(rg3p_approx(alpha_t, beta_t / mu(h), (int)IW_df), 2);
      // tau(h) = 1.0 / std::pow(rg3p(alpha_t, beta_t / mu(h), (int)IW_df), 2);
      zeta(h) = 1.0 / R::rgamma(1.0, 1.0 / (1.0 + 1.0 / tau(h)));
    }
    
    // ::::::::::::::::::::::::::::::::::::::::::::
    // sample R
    
    IW_scale.setZero();
    for (j = 1; j < p; j++)
    {
      for (js = 0; js < j; js++)
      {
        for (h = 0; h < K; h++)
          e_ij(h) = diag_V(h) * Omega[js][j][h] / std::sqrt(tau(h) * Lambda[js][j][h]);
        IW_scale += e_ij * e_ij.transpose();
      }
    }
    
    Psi = rIW(IW_df, IW_scale, K);
    e_ij = Psi.diagonal();
    e_ij = sinverse(e_ij, K);
    R_s = e_ij.asDiagonal() * Psi * e_ij.asDiagonal();
    
    a_MH = exp((K+1) * 0.5 * (std::log(R_s.determinant()) - std::log(R.determinant())));
    
    if (R::runif(0, 1) < a_MH)
    {
      R = R_s;
      
      for (h = 0; h < K; h++)
      {
        Rinv = delete_rc_from_eigen(R, h, Km1);
        Rinv = Rinv.inverse();
        rtR.col(h) = delete_elem_from_eigen_vec(R.col(h), h, Km1);
        rtR.col(h) = rtR.col(h).transpose() * Rinv;
        
        mu(h) = 1.0 - rtR.col(h).transpose() * delete_elem_from_eigen_vec(R.col(h), h, Km1);
      }
    }
    
    for (h = 0; h < K; h++)
    {
      ta_s = R::rbeta(hyp_ta(0), hyp_ta(1));
      l2 = 0.0;
      
      mth = mu(h) * tau(h);
      mth2 = std::sqrt(tau(h)) * mu(h);
      
      cs = 0;
      for (j = 0; j < p-1; j++)
      {
        for (js = j+1; js < p; js++)
        {
          c = 0;
          for (hs = 0; hs < K; hs++)
          {
            if (hs != h)
            {
              delta_12_inv(c) = Omega[js][j][hs] / std::sqrt(tau(hs) * Lambda[js][j][hs]);
              c += 1;
            }
          }
          
          tmp1 = Omega[js][j][h];
          alpha_l = std::sqrt(0.5 * tmp1 * tmp1 / mth + 1.0 / Eta[js][j][h]);
          beta_l = (rtR.col(h).transpose() * delta_12_inv).value() * tmp1 / mth2;
          
          qta2(cs) = 1 - pgamma3p1(std::sqrt((1.0 - ta_s) / ta_s), alpha_l, beta_l);
          if (Z(cs, h) == 1)
          {
            l2 += std::log(qta2(cs));
          } else 
          {
            l2 += std::log(1.0 - qta2(cs));
          }
          cs += 1;
          
          a_MH = l2 - l1(h);
          if (std::log(R::runif(0, 1)) < a_MH)
          {
            ta(h) = ta_s;
            l1(h) = l2;
            qta1.col(h) = qta2;
          }
        }
      }
      
      for (j = 0; j < IW_df; j++) 
      {
        Z(j, h) = R::rbinom(1, qta1(j, h));
      }
    }
    
    l1_out.row(b) = l1;
    
    if (b >= bn)
    {
      tau_out += tau;
      R_out += R;
      ta_out += ta;
      
      cs = 0;
      for (h = 0; h < K; h++)
      {
        c = 0;
        for (j = 0; j < p; j++)
        {
          Omega_jj_out(j, h) += Omega[j][j][h];
          for (js = j+1; js < p; js++)
          {
            Omega_ij_out(c, h) += Omega[js][j][h];
            lambda_out(c, h) += Lambda[js][j][h];
            lambda_prob(c, h) += 1.0 / (1.0 + Lambda[js][j][h]);
            eta_out(c, h) += Eta[js][j][h];
            c += 1;
            cs += 1;
          }
        }
      }
    }
    
    if ((ping != 0) && (b % ping == 0))
      Rcpp::Rcout << "\n iter: " << b;
  }
  
  // ::::::::::::::::::::::::::::::::::::::::::::
  // output
  // ::::::::::::::::::::::::::::::::::::::::::::
  
  return Rcpp::List::create(Rcpp::Named("Omega_ij") = Omega_ij_out / (B - bn), Rcpp::Named("Omega_jj") = Omega_jj_out / (B - bn), Rcpp::Named("ta") = ta_out / (B - bn), Rcpp::Named("l1_chain") = l1_out, Rcpp::Named("lambda") = lambda_out / (B - bn), Rcpp::Named("lambda_prob") = lambda_prob / (B - bn), Rcpp::Named("eta") = eta_out / (B - bn), Rcpp::Named("R") = R_out / (B - bn), Rcpp::Named("tau") = tau_out / (B - bn));
}

// end file

