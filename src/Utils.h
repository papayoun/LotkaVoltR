// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(mvtnorm)]]
#include <mvtnormAPI.h>
#ifndef DebugMethods_H
#define DebugMethods_H

class Utils{
private:
  static arma::vec triangl(const arma::mat& X){
    int n = X.n_cols;
    arma::vec res(n * (n-1) / 2);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < i; ++j) {
        res(j + i * (i-1) / 2) = X(i, j);
      }
    }
    return res;
  }
public:
  static arma::mat kalmanFiltVar(const arma::mat& cov1, const arma::mat& cov2){
    arma::mat prec1 = arma::inv_sympd(cov1);
    arma::mat prec2 = arma::inv_sympd(cov2);
    arma::mat out = arma::inv_sympd(prec1 + prec2);
    return out;
  }
  static arma::colvec kalmanFiltMean(const arma::colvec& mu1, const arma::colvec& mu2,
                              const arma::mat& prec1, const arma::mat& prec2,
                              const arma::mat& kalCov){
    return kalCov * (prec1 * mu1 + prec2 * mu2);
  }
  static void checkSDP2d(const arma::mat& mat, std::string name = ""){
    double det = mat(0, 0) * mat(1, 1) - mat(0, 1) * mat(1, 0);
    double diffDiag = std::abs(mat(0, 1) - mat(1, 0));
    bool out = ((mat(0, 0) > 0) and (mat(1, 1) > 0) and (det > 0) and
                  (diffDiag < pow(10, -10)));
    if(not out){
      std::cout << "Matrix " << name << ": " << std::endl;
      DebugMethods::debugprint(mat, "la matrice");
      Rcpp::stop("The matrix is not semi positive definite");
    }
  }
  static arma::mat outerprod(const arma::colvec & x) {
    arma::mat m = x * x.t();
    return m;
  }
  static double innerprod(const arma::colvec & x, const arma::colvec & y) {
    return arma::sum(x % y);
  }
  static double frobprod(const arma::mat & x, const arma::mat & y){
    return arma::accu(x % y);
  }
  
  static arma::colvec rmvnorm(const arma::colvec& mu, 
                              const arma::mat& sigma) {
    arma::colvec Y(size(mu)); Y.randn();
    return mu + arma::chol(sigma, "lower") * Y;
  }
  
  static arma::mat rmvnorm(unsigned int n, const arma::colvec& mu, 
                           const arma::mat& sigma) {
    arma::mat out(mu.size(), n);
    for(unsigned int i = 0; i < n; i ++){
      out.col(i) = rmvnorm(mu, sigma);
    }
    return out.t();// returns a n * 2 matrix
  }
  static double dmvnorm(const arma::colvec& x,  
                        const arma::colvec& mean,  
                        const arma::mat& sigma, 
                        bool logd = false) { 
    int xdim = x.n_rows;
    arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
    double rootisum = arma::sum(log(rooti.diag()));
    double constants = - (xdim / 2.0) * std::log(2.0 * M_PI);
    arma::colvec z = rooti * (x - mean);    
    double out = constants - 0.5 * innerprod(z, z) + rootisum;     
    if (!logd) {
      out = exp(out);
    }
    return(out);
  }
  static bool checkAllPos(const arma::colvec& hiddenState){
    bool cond = true;
    unsigned int n = hiddenState.size();
    for(unsigned int i = 0; i < n; i ++){
      if(hiddenState(i) <= 0){
        cond = false;
        break;
      }
    }
    return cond;// Returns true if all the coordinates are strictly positive
  };
  static bool checkAllPos(const Rcpp::NumericVector& vec){
    bool cond = true;
    unsigned int n = vec.size();
    for(unsigned int i = 0; i < n; i ++){
      if(vec(i) <= 0){
        cond = false;
        break;
      }
    }
    return cond;// Returns true if all the coordinates are strictly positive
  };
  static bool checkAllPosOr0(const Rcpp::NumericVector& vec){
    bool cond = true;
    unsigned int n = vec.size();
    for(unsigned int i = 0; i < n; i ++){
      if(vec(i) < 0){
        cond = false;
        break;
      }
    }
    return cond;// Returns true iff all the coordinates are strictly positive
  };
  static unsigned int countNeg(Rcpp::NumericVector vec){
    unsigned int out = 0;
    unsigned int n = vec.length();
    for(unsigned int i = 0; i < n; i++){
      if(vec(i) < 0)
        out += 1;
    }
    return out;
  }
  static arma::colvec vec_log(const arma::colvec& x){
    unsigned int n = x.size();
    arma::colvec output(n);
    for(unsigned int i = 0; i < n; i ++){
      if(x (i) <= 0){
        Rcpp::stop("You're trying to take the log of a non positive number");
      }
      output(i) = log(x(i));
    }
    return output;
  }
  static arma::colvec rmvlognorm(const arma::colvec& mu, 
                           const arma::mat& sigma){
    arma::colvec logOutput = rmvnorm(mu, sigma);
    return exp(logOutput);
  }
  static arma::mat rmvlognorm(unsigned int n,
                              const arma::colvec& mu, 
                              const arma::mat& sigma){
    arma::mat out(mu.size(), n);
    for(unsigned int i = 0; i < n; i ++){
      out.col(i) = rmvlognorm(mu, sigma);
    }
    return out.t();
  }
  static double dmvlognorm(const arma::colvec& y, const arma::colvec& mu,
                           const arma::mat& sigma,
                           bool logd = false){
    arma::colvec logY = vec_log(y);
    double sum_log_inv = -arma::sum(vec_log(y));
    double dnorm_log = dmvnorm(logY, mu, sigma, true);
    double out = sum_log_inv + dnorm_log;
    if (!logd) {
      out = exp(out);
    }
    return(out);
  }
  static double pmvnorm_cpp(const arma::colvec& lower_bound,
                            const arma::colvec& upper_bound,
                            const arma::colvec& mean,
                            const arma::mat& cov,
                            double abseps = 1e-3){
    int n = lower_bound.n_elem;
    for(unsigned int i = 0; i < n; i ++){
      if(lower_bound(i) >= upper_bound(i)){
        Rcpp::stop("One element of lower_bound is >= than upper_bound");
      }
    }
    arma::vec sdVec(2);
    sdVec(0) = pow(cov(0, 0), -0.5); sdVec(1) = pow(cov(1, 1), -0.5);
    arma::mat sdMat = arma::diagmat(sdVec);
    arma::mat corMat = sdMat * cov * sdMat;
    arma::vec lowertrivec = triangl(corMat);
    int nu = 0;
    int maxpts = 25000;     // default in mvtnorm: 25000
    double releps = 0;      // default in mvtnorm: 0
    int rnd = 1;            // Get/PutRNGstate
    arma::vec new_lb =  sdMat * (lower_bound - mean); 
    arma::vec new_ub =  sdMat * (upper_bound - mean); 
    double* lower_bound_ = new_lb.memptr();
    double* upper_bound_ = new_ub.memptr();
    double* correlationMatrix = lowertrivec.memptr();
    int* infin = new int[n];
    double* delta = new double[n];
    for (int i = 0; i < n; ++i) {
      infin[i] = 2; // as in mvtnorm:::mvt
      if(std::isinf(lower_bound(i))){
        infin[i] = 0;
        if(std::isinf(upper_bound(i))){
          infin[i] = -1;
        }
      }
      else if(std::isinf(upper_bound(i))){
        infin[i] = 1;
      }
      delta[i] = 0.0;
    }
    // return values
    double error;
    int inform;
    double value;
    
    mvtnorm_C_mvtdst(&n, &nu, lower_bound_, upper_bound_,
                     infin, correlationMatrix, delta,
                     &maxpts, &abseps, &releps,
                     &error, &value, &inform, &rnd);
                     delete[] (infin);
                     delete[] (delta);
                     
    return value;
  }
  static double qmvnorm(const arma::colvec& lower_bound,
                        const arma::colvec& mean,
                        const arma::mat& cov,
                        double abseps = 1e-3){
    arma::colvec up(2); up.fill(R_PosInf);
    return pmvnorm_cpp(lower_bound, up, mean, cov, abseps);
  }
  static arma::colvec rmvnorm_low_trunc(const arma::colvec& mu, 
                                     const arma::mat& sigma,
                                     const arma::colvec& lower_bound) {
    bool cond_arret = false;
    arma::colvec out(mu.size());
    while(!cond_arret){
      out = rmvnorm(mu, sigma);
      cond_arret = checkAllPos(out - lower_bound);
    }
    return out;
  }
  static arma::mat rmvnorm_low_trunc(unsigned int n,
                                     const arma::colvec& mu, 
                           const arma::mat& sigma,
                           const arma::colvec& lower_bound) {
    arma::mat out(mu.size(), n);
    for(unsigned int i = 0; i < n; i ++){
      out.col(i) = rmvnorm_low_trunc(mu, sigma, lower_bound);
    }
    return out.t();// returns a n * 2 matrix
  }
  static double dmvnorm_low_trunc(const arma::colvec& x,
                               const arma::colvec& mu, 
                               const arma::mat& sigma,
                               const arma::colvec& lower_bound) {
    double constant = qmvnorm(lower_bound, mu, sigma);
    double out = dmvnorm(x, mu, sigma) / constant;
    return out;
  }
  // Sampling with replacement
  static Rcpp::IntegerVector findInterval(Rcpp::NumericVector& x,
                                          Rcpp::NumericVector& breaks) {
    // equivalent of findInterval function
    Rcpp::IntegerVector out(x.size());
    Rcpp::NumericVector::iterator it, pos;
    Rcpp::IntegerVector::iterator out_it;
    for(it = x.begin(), out_it = out.begin(); it != x.end();++it, ++out_it) {
      pos = std::upper_bound(breaks.begin(), breaks.end(), *it);
      *out_it = std::distance(breaks.begin(), pos);
    }
    return out;
  }
  static Rcpp::NumericVector sampleReplace(const Rcpp::NumericVector& x,
                                           const int& size,
                                           const Rcpp::NumericVector& probs){
    int nx=x.size();int np=probs.size();
    if(nx!=np){
      Rcpp::stop("Error, x and probs must have same length");
    }
    Rcpp::NumericVector cumprob(nx + 1);
    cumprob.fill(0);
    for(int i = 1; i < (nx + 1); i++){
      cumprob(i) = cumprob(i - 1) + probs(i - 1);
    }
    Rcpp::NumericVector us = Rcpp::runif(size, 0, cumprob(nx));
    Rcpp::IntegerVector inds = findInterval(us, cumprob) - 1;
    return x[inds];
  }
  static Rcpp::IntegerVector sampleReplace(const Rcpp::IntegerVector& x,
                                           const int& size,
                                           const Rcpp::NumericVector& probs){
    int nx=x.size();int np=probs.size();
    if(nx!=np){
      Rcpp::stop("Error, x and probs must have same length");
    }
    Rcpp::NumericVector cumprob(nx+1);
    cumprob.fill(0);
    for(int i=1;i<(nx+1);i++){
      cumprob(i) = cumprob(i-1) + probs(i-1);
    }
    Rcpp::NumericVector us = Rcpp::runif(size, 0, cumprob(nx));
    Rcpp::IntegerVector inds = findInterval(us, cumprob) - 1;
    return x[inds];
  }
};

#endif