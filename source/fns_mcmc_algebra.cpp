// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]


// Beta update function
//
// [[Rcpp::export]]
Rcpp::List penalized_beta_update(
    const Eigen::VectorXd & y,
    const Eigen::MatrixXd & X,
    const Eigen::MatrixXd & lambda_S,
    const double & sigmasq
) {

  int N = X.rows();
  int K = X.cols();

  Eigen::MatrixXd xtx_plus_pen(K, K);
  xtx_plus_pen = X.transpose() * X + lambda_S;

  Eigen::MatrixXd beta_helper(K, K);
  beta_helper = xtx_plus_pen.inverse();

  Eigen::VectorXd beta_expect(K);
  beta_expect = beta_helper * X.transpose() * y;

  Eigen::MatrixXd beta_var(K, K);
  beta_var = sigmasq * beta_helper;

  return Rcpp::List::create(Rcpp::Named("coeff_mean") = beta_expect,
                            Rcpp::Named("coeff_var") = beta_var,
                            Rcpp::Named("xtx_pen") = xtx_plus_pen);
}

// Beta update In Place
//
// [[Rcpp::export]]
void penalized_beta_inplace(
    const Eigen::VectorXd & y,
    const Eigen::MatrixXd & X,
    const Eigen::MatrixXd & lambda_S,
    const double & sigmasq,
    Eigen::Map<Eigen::VectorXd> & beta_expect,
    Eigen::Map<Eigen::MatrixXd> & beta_var,
    Eigen::Map<Eigen::MatrixXd> & xtx_plus_pen
) {

  int N = X.rows();
  int K = X.cols();

  xtx_plus_pen = X.transpose() * X + lambda_S;

  Eigen::MatrixXd beta_helper(K, K);
  beta_helper.block(0, 0, K, K) = xtx_plus_pen.inverse();

  beta_expect = beta_helper * X.transpose() * y;

  beta_var.block(0, 0, K, K) = sigmasq * beta_helper;
}

// Multivariate Normal random draw In Place
//
// [[Rcpp::export]]
void normal_beta_draw(
    const Eigen::VectorXd & expect,
    const Eigen::MatrixXd & covar,
    Eigen::VectorXd & beta,
    Eigen::MatrixXd & beta_var_L
){
  int K = expect.size();

  beta_var_L = covar.llt().matrixL();

  Eigen::VectorXd z_sample(K);
  for (auto i=0; i<K; i++)
    z_sample(i) = R::rnorm(0, 1);

  beta = beta_var_L * z_sample + expect;
}

// Sigma squared update function In Place
//
// [[Rcpp::export]]
void penalized_sigma_sq_inplace(
    const Eigen::VectorXd & y,
    const Eigen::MatrixXd & X,
    const Eigen::VectorXd & beta,
    const Eigen::MatrixXd & xtx_plus_pen,
    const double & b_pri,
    double & b_post
) {
  b_post = b_pri -
    y.transpose() * X * beta +
    0.5 * y.transpose() * y +
    0.5 * beta.transpose() * xtx_plus_pen * beta;
}
