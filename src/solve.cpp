#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::VectorXd chol_solve2(
  const Eigen::MatrixXd& M,
  const Eigen::VectorXd& b
){

  const Eigen::LLT<Eigen::MatrixXd> LLt(M);
  return LLt.solve(b);
}
