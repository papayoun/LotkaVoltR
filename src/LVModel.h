#include <RcppArmadillo.h>      // declarations for both Rcpp and RcppArmadillo offering Armadillo classes
// // [[Rcpp::depends(RcppArmadillo)]]
#ifndef Utils_H
#define Utils_H

//' @title Class SDE model
//' @name LV_simple
//' @description bivariate SDE LV model
//' @export LVModel
class LVModel{  
private:
  ////////////// ATTRIBUTES /////////////////////////////////
  // Model parameters
  arma::colvec preyParam;
  arma::colvec predParam;
  arma::mat gamma;
  arma::mat param;
  arma::colvec mu0;
  arma::mat sigma0;
  // Density computation tuning parameters
  arma::colvec lb = arma::zeros(2);
  arma::mat twoAndOnes;
  double intens_delta = 5;
  double alpha_pow = 0.28;
  void createParam(){
    arma::mat par_tmp(2, 3);
    par_tmp(0, 0) = preyParam(0); 
    par_tmp(0, 1) = -preyParam(1); 
    par_tmp(0, 2) = -preyParam(2); 
    par_tmp(1, 0) = -predParam(0); 
    par_tmp(1, 1) =  predParam(1); 
    par_tmp(1, 2) = -predParam(2); 
    param = par_tmp;
  }
  arma::colvec getDrift(const arma::colvec& currentState) const{
    arma::colvec output(2);
    output(0) = currentState(0) * (param(0, 0) + param(0, 1) * currentState(0) + 
      param(0, 2) * currentState(1));
    output(1) = currentState(1) * (param(1, 0) + param(1, 1) * currentState(0) + 
      param(1, 2) * currentState(1));
    return output;
  }
  arma::mat getDiffusion(const arma::colvec& currentState) const{
    arma::mat output(2, 2);
    output(0, 0) = currentState(0) * gamma(0, 0);
    output(0, 1) = currentState(0) * gamma(0, 1);
    output(1, 0) = currentState(1) * gamma(1, 0);
    output(1, 1) = currentState(1) * gamma(1, 1);
    return output;
  }
  arma::mat getBrownianVariance(const arma::colvec& currentState) const{
    arma::mat sigma = getDiffusion(currentState);
    return sigma * sigma.t();
  }
  double getB_1(const arma::colvec& currentState) const{
    arma::colvec output(2);
    output(0) = (param(0, 0) + 2 * param(0, 1) * currentState(0) + 
      param(0, 2) * currentState(1));
    output(1) = (param(1, 0) + param(1, 1) * currentState(0) + 
      2 * param(1, 2) * currentState(1));
    return arma::sum(output);
  };
public:
  // Construction method
  LVModel(){};
  LVModel(const arma::colvec& a1, const arma::colvec& a2, 
          const arma::colvec& m0, const arma::mat& s0, const arma::mat& gam)
    : preyParam(a1), predParam(a2), mu0(m0), sigma0(s0), gamma(gam){
    Utils::checkSDP2d(sigma0, "sigma0");
    arma::mat tmp(2, 2); tmp.fill(1);
    tmp(0, 0) = tmp(1, 1) = 2;
    twoAndOnes = tmp;
    for(unsigned int i = 0; i < 3; i ++){
      bool cond1 = preyParam(i) > 0;
      bool cond2 = predParam(i) > 0;
      if(!(cond1 & cond2)){
        Rcpp::stop("All parameters of the pred-prey dynamics must be > 0");
      };
    }
    createParam();
  };
  arma::colvec getEulerMean(const arma::colvec& start, 
                            const double time_lag) const{
    return start + getDrift(start) * time_lag;
  }
  arma::colvec rLV_euler(const arma::colvec& start, const double time_lag) const{
    arma::colvec mean = getEulerMean(start, time_lag);
    arma::mat covar = getBrownianVariance(start) * time_lag;
    return Utils::rmvnorm(mean, covar); 
  }
  double dLV_euler(const arma::colvec& start, const arma::colvec& end,
                         const double time_lag) const{
    arma::colvec mean = getEulerMean(start, time_lag);
    arma::mat covar = getBrownianVariance(start) * time_lag;
    return Utils::dmvnorm(end, mean, covar); 
  }
  arma::mat simulateTrajectory(const arma::colvec& simulationTimes) const{
    arma::colvec x0 = Utils::rmvnorm_low_trunc(mu0, sigma0, lb);
    unsigned int simulationSize = simulationTimes.size();
    arma::mat output(2, simulationSize); output.col(0) = x0;
    for(unsigned int i = 1; i < simulationSize; i++){
      double delta = simulationTimes(i) - simulationTimes(i - 1);
      output.col(i) = rLV_euler(output.col(i - 1), delta);  
    }
    return output.t();
  };
  arma::mat getParam() const{
    return param;
  };
  double getIntensity(const double u) const{
    return intens_delta  * pow(u, alpha_pow - 1);
  }
  double rhoWeight(const arma::colvec& x, const arma::colvec y,
                   const double u, const double intens) const{
    arma::mat sigmaX = getBrownianVariance(x);
    // std::cout << "u vaut " << u << std::endl;
    // std::cout << "intens vaut " << intens << std::endl;
    arma::mat sigmaY = getBrownianVariance(y);
    arma::mat sigmaInvX = arma::inv_sympd(sigmaX);
    arma::colvec bx = getDrift(x);
    arma::colvec eul_inc = y - x - u * bx;
    arma::colvec lambda = - sigmaInvX * eul_inc / u;// Lambda_k
    arma::mat k_theta = Utils::outerprod(lambda) - sigmaInvX / u;
    arma::colvec vecOnes(2); vecOnes.fill(1);
    arma::mat matY(2, 2); matY.col(0) = y; matY.col(1) = y;
    arma::mat sigmaY_1 = (twoAndOnes % matY) % gamma * gamma.t();
    double t1 = 0.5 * Utils::frobprod(sigmaY - sigmaX, k_theta);
    // std::cout << "t1 vaut " << t1 << std::endl;
    double t2 = 0.5 * Utils::frobprod(twoAndOnes, gamma * gamma.t());
    // std::cout << "t2 vaut " << t2 << std::endl;
    arma::colvec term3a = sigmaY_1 * vecOnes - getDrift(y) + bx;
    double t3 = Utils::innerprod(term3a, lambda);
    // std::cout << "t3 vaut " << t3 << std::endl;
    double t4 = - getB_1(y);
    // std::cout << "t4 vaut " << t4 << std::endl;
    // std::cout << "c(" << t1 << ", "<< t2 << ", " << t3 << ", " << t4 << ") / " << intens <<std::endl;
    return 1 + (t1 + t2 + t3 + t4) / intens;
  }
  double getDensX0(arma::colvec x0){
    return Utils::dmvnorm_low_trunc(x0, mu0, sigma0, lb);
  }
  double getUnbiasSignDensEstim(const arma::colvec& x0, const arma::colvec& xT,
                                const double time_lag) const{
    double current_time = 0;
    double weight = 1; // Initial weight
    arma::colvec current_start(2); // Current value of the skeleton
    current_start(0) = x0(0); current_start(1) = x0(1);
    int compt = 0;
    while(true){
      double unif = Rcpp::runif(1, 0, 1)(0);
      double delta_tau = pow(-alpha_pow / intens_delta * log(unif), 1 / alpha_pow);
      if(current_time + delta_tau > time_lag){
        break;
      }
      else{
        compt += 1;
        arma::colvec current_end = rLV_euler(current_start, delta_tau);
        double current_intens = getIntensity(delta_tau);
        double rho_weight = rhoWeight(current_start, current_end, delta_tau, 
                                      current_intens);
        weight *= rho_weight;
        // std::cout << "weight at compt " << compt << ": " << weight << std::endl;
        current_time += delta_tau;
        current_start = current_end;
      }
    }
    double delta_tau = time_lag - current_time;
    // std::cout << "Il y a eu tant de temps:" << compt << std::endl;
    // std::cout << "current_start <- c(" << current_start(0) << ", " << current_start(1) << ")" << std::endl;
    // std::cout << "delta_tau <-  " <<  delta_tau << std::endl;
    double last_q = dLV_euler(current_start, xT, delta_tau);
    // std::cout << "last_q <-  " <<  last_q << std::endl;
    return weight * dLV_euler(current_start, xT, delta_tau);
  }
  double unbiasedDensityEstimate(const arma::colvec& x0, 
                                 const arma::colvec& xT,
                                 const double time_lag,
                                 const unsigned int sampleSize) const{
    bool allPos = Utils::checkAllPos(xT);
    if(not allPos){
      return(0); // Simulated vector must be positive
    }
    Rcpp::NumericVector samples(sampleSize);
    for(unsigned int i = 0; i < sampleSize; i++){
      samples(i) = getUnbiasSignDensEstim(x0, xT, time_lag);
    }
    return Rcpp::mean(samples);
  }
  // GETTERS
  arma::mat getSigma0() const{
    return sigma0;
  }
  arma::colvec getMu0() const{
    return mu0;
  }
};

// Exposes the class to Rcpp
RCPP_MODULE(LV_Module) {
  using namespace Rcpp;
  class_<LVModel>("LVModel")
    .constructor<arma::colvec, arma::colvec, arma::colvec, arma::mat, arma::mat>("constructor")
    .method("simulate", &LVModel::simulateTrajectory)
    .method("signed_density", &LVModel::getUnbiasSignDensEstim)
    .method("euler_density", &LVModel::dLV_euler)
  ;
}

#endif