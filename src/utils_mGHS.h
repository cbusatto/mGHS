/* File: utils_mGHS.h */

#ifndef utils_mGHS_H
#define utils_mGHS_H

double pgamma3p1 (const double& x, const double& a, const double& b);
  
Eigen::VectorXd inverse (Eigen::VectorXd& x, const int& n);
Eigen::VectorXd sinverse (Eigen::VectorXd& x, const int& n);
Eigen::VectorXd delete_elem_from_boost_vec (const boost::multi_array<double, 1>& x, const int& j, const int& nm1);
Eigen::VectorXd delete_elem_from_eigen_vec (const Eigen::VectorXd& x, const int& j, const int& nm1);
Eigen::VectorXd get_s2 (const boost::multi_array<double, 1>& s_h, const int& p);
Eigen::VectorXd get_diag_vec_h (const Eigen::VectorXd& tau, const boost::multi_array<double, 1>& lambda, const int& h, const int& Km1);
Eigen::VectorXd get_eij (const Eigen::VectorXd& tau, const boost::multi_array<double, 1>& lambda, const Eigen::VectorXd& diag_V, const boost::multi_array<double, 1>& omega, const int& K);
Eigen::VectorXd rmvnorm (Eigen::VectorXd& v, const int& j, const double& s_22, const Eigen::VectorXd& s_12, const Eigen::MatrixXd& inv_Omega_11, const Eigen::VectorXd& d_inv, const Eigen::VectorXd& temp_2, const int& pm1);

Eigen::MatrixXd delete_rc_from_eigen (const Eigen::MatrixXd& X, const int& j, const int& pm1);
Eigen::MatrixXd delete_rc_from_boost1 (const boost::multi_array<double, 2>& X, const int& j, const int& pm1);  
Eigen::MatrixXd delete_rc_from_boost2 (boost::multi_array<double, 2>& X, const int& j, const int& pm1);  
Eigen::MatrixXd get_matrix (const int& j, const boost::multi_array<double, 2>& X, const int& p, const int& K);
Eigen::MatrixXd get_Sigma_h (boost::multi_array<double, 2>& X, const int& p);
Eigen::MatrixXd get_delta (const Eigen::VectorXd& tau, const Eigen::MatrixXd& lambda_12, const int& pm1, const int& K);
Eigen::MatrixXd bb (const boost::multi_array<double, 3>& A, const int& p, const int& h) ;
  
Eigen::MatrixXd get_Sigma2 (Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >& Sigma, const Eigen::MatrixXd& inv_Omega_11, const Eigen::VectorXd& v, const double& v22, const int& j, const int& pm1);
    
boost::multi_array<double, 2> get_Sigma (boost::multi_array_ref<double, 3>::array_view<2>::type Sigma, const Eigen::MatrixXd& inv_Omega_11, const int& j, const int& pm1, const int& p);
boost::multi_array<double, 2> fill_S (const Eigen::MatrixXd& S, const int& p);
boost::multi_array<double, 2> fill_Xt (const Eigen::MatrixXd& S, const int& n, const int& p);

std::tuple<Eigen::VectorXd, Eigen::VectorXd> lambda_pstep (const Eigen::VectorXd& v, const Eigen::VectorXd& temp_1, const Eigen::VectorXd& eta_12, const double& mth, const double& mth2);
  
#endif

