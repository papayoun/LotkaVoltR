#include <RcppArmadillo.h>      // declarations for both Rcpp and RcppArmadillo offering Armadillo classes
// [[Rcpp::depends(RcppArmadillo)]]
#ifndef LVModel_H
#define LVModel_H

//' @title Class of partially diffusion process ruled by a LV model
//' @name LVPOD
//' @description POD with the bivariate SDE LV model
//' and the observation model Y_t = X_t * exp(N(0, Sigma_obs)) 
//' @export LVPOD
class LVPOD{
private:
  LVModel model;
  arma::mat observationCov;
  arma::colvec meanNoise;
  arma::colvec detectVec;
  arma::mat detectMat;
  void computeMeanNoise(){
    meanNoise = - 0.5 * arma::diagvec(observationCov);
  }
  void computeDetectMat(){
    unsigned int p = detectVec.size();
    arma::mat tmp(p, p, arma::fill::zeros);
    for(unsigned int i = 0; i < p; i ++){
      tmp(i, i) = detectVec(i);
    }
    detectMat = tmp;
  }
  arma::colvec getMean(const arma::colvec& hiddenState) const{
    return Utils::vec_log(detectMat * hiddenState) + meanNoise;
  }
public:
  LVPOD()
    : model(){
    arma::mat tmp(2, 2, arma::fill::eye);
    observationCov = tmp;
    Utils::checkSDP2d(observationCov, "obsCov");
    arma::colvec tmp1(2, arma::fill::randu);
    detectVec = tmp1;
    computeMeanNoise();
    computeDetectMat();
  }
  LVPOD(const arma::colvec& a1, const arma::colvec& a2, 
        const arma::colvec& mu0, const arma::mat& sigma0,
        const arma::mat& gam,
        const arma::mat& cov, const arma::colvec& qs)
    : model(a1, a2, mu0, sigma0, gam), observationCov(cov), detectVec(qs){
    Utils::checkSDP2d(observationCov, "obsCov");
    computeMeanNoise();
    computeDetectMat();
  }
  Rcpp::List simulateHSandObs(const Rcpp::NumericVector& simulationTimes){
    arma::mat hiddenStates = model.simulateTrajectory(simulationTimes);
    arma::mat noise = Utils::rmvnorm(simulationTimes.size(),
                                               meanNoise,
                                               observationCov);
    arma::mat observations = (hiddenStates * detectMat) % exp(noise);
    return  Rcpp::List::create(Rcpp::Named("observations") = observations,
                               Rcpp::Named("hiddenStates") = hiddenStates,
                               Rcpp::Named("simulationTimes") = simulationTimes);
  };
  void setObsCov(const double newSigma){observationCov = newSigma;};
  double unbiasedDensity(const arma::colvec& x0, const arma::colvec& xF, 
                         double time_lag, 
                         const unsigned int sampleSize) const{
    return model.unbiasedDensityEstimate(x0, xF, time_lag, sampleSize);
  }
  double observationDensityUnit(const arma::colvec& hiddenState, 
                                const arma::colvec& observation) const{
    bool allPos = Utils::checkAllPos(hiddenState);
    if(not allPos){
      return 0;
    }
    arma::colvec mu = getMean(hiddenState);
    return Utils::dmvlognorm(observation, mu, observationCov, false);
  }
  Rcpp::NumericVector observationDensity(const arma::mat& hiddenStates,
                                         const arma::colvec& observation) const{
    unsigned int particleSize = hiddenStates.n_rows; // a generation aof particle
    // is a n * 2 matrix
    Rcpp::NumericVector output(particleSize);
    for(unsigned int i = 0; i < particleSize; i++){
      arma::colvec colPart = hiddenStates.row(i).t();
      output(i) = observationDensityUnit(colPart, observation);
    }
    return output;
  }
  arma::mat getSigmaObs() const{
    return observationCov;
  }
  double getDensX0(arma::colvec x0){
    return model.getDensX0(x0);
  }
  arma::colvec getEulerMean(const arma::colvec& start, const double time_lag) const{
    return model.getEulerMean(start, time_lag);
  }
  arma::colvec getScaledObs(const arma::colvec& obs) const{
    arma::colvec out(2);
    out(0) = obs(0) / detectVec(0);
    out(1) = obs(1) / detectVec(1);
    return out;
  }
  double getEulerDensity(const arma::colvec& start, const arma::colvec& end,
                         const double time_lag){
    return model.dLV_euler(start, end, time_lag);
  }
  // GETTERS
  LVModel getModel() const{
    return model;
  }
};

// Exposes the class to Rcpp
RCPP_MODULE(LVPOD_Module) {
  using namespace Rcpp;
  class_<LVPOD>("LVPOD")
    .constructor<arma::colvec, arma::colvec, arma::colvec, arma::mat, arma::mat, arma::mat, arma::colvec>("constructor")
    .method("simulate", &LVPOD::simulateHSandObs)
    .method("density", &LVPOD::unbiasedDensity)
    .method("obs_density", &LVPOD::observationDensityUnit)
    .method("euler_density", &LVPOD::getEulerDensity)
  ;
}

#endif