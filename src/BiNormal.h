#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

#include <RcppNumerical.h>
using namespace Numer;

// P(a1 < X1 < b1, a2 < X2 < b2), (X1, X2) ~ N([0], [1   rho])
//                                            ([0], [rho   1])
class BiNormal: public MFunc
{
private:
  const arma::mat prec;
  const arma::colvec mu;
  double constant;  // sqrt(det(precision)) * 0.5 / M_PI
public:
  BiNormal(const arma::colvec& mu_, const arma::mat precision_)
  : mu(mu_), prec(precision_){
    double det_prec = prec(0, 0) * prec(1, 1) - prec(1, 0) * prec(0, 1);
    constant = sqrt(prec) * 0.5 / M_PI;
  }
  // PDF of bivariate normal
  double operator()(Constvec& x)
  {
    double z = (x - mu);
    return const2 * std::exp(-z / const1);
  }
};