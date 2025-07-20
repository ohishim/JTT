#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

template<typename EigenType>
std::vector<EigenType> list2vec(const Rcpp::List& L, const int m){

  std::vector<EigenType> out(m);

  for (int j = 0; j < m; ++j) {
    out[j] = Rcpp::as<EigenType>(L[j]);
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector JTT_score(
    const int q,
    const int m,
    const Eigen::MatrixXi& E,
    const Rcpp::List& Xy_,
    const Eigen::VectorXd& yPy,
    const Rcpp::List& M_,
    const double D0,
    const double alp
){

  const std::vector<Eigen::VectorXd> Xy = list2vec<Eigen::VectorXd>(Xy_, m);
  const std::vector<Eigen::MatrixXd> M = list2vec<Eigen::MatrixXd>(M_, m);

  Rcpp::NumericVector out(q);

  for(int j=0; j < q; ++j){

    const Eigen::VectorXi& Ej = E.row(j);
    const int k = Ej[0]; const int l = Ej[1];

    const Eigen::VectorXd& Xy_k = Xy[k];
    const Eigen::VectorXd& Xy_l = Xy[l];
    const Eigen::VectorXd Xy_kl = Xy_k + Xy_l;

    const Eigen::MatrixXd& M_k = M[k];
    const Eigen::MatrixXd& M_l = M[l];
    const Eigen::LLT<Eigen::MatrixXd> LLTD(M_k + M_l);

    const double D1 = yPy(Ej).sum() - LLTD.solve(Xy_kl).dot(Xy_kl);

    out[j] = D0*D1 - alp;
  }

  return(out);
}
