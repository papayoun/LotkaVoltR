#include <RcppArmadillo.h>      // declarations for both Rcpp and RcppArmadillo offering Armadillo classes
#ifndef LVPOD_H
#define LVPOD_H

//' @title ProposalModel for LVPOD
//' @name ProposalLVModel
//' @description Proposoal model with optimal filter
//' POD with the bivariate SDE LV model
//' and the observation model Y_t = X_t * exp(N(0, Sigma_obs))
//' @export ProposalLVModel
class ProposalLVModel{
private:
  LVPOD trueModel;
  arma::mat randomWalkCov ;
  arma::mat randomWalkPrec;
  arma::mat aroundObsCov  ;
  arma::mat aroundObsPrec ;
  arma::colvec lb = arma::zeros(2);
  double weightDyn;
  double weightObs;
  void computeAroundObsCov(){
    arma::mat id(2, 2, arma::fill::eye);
    aroundObsCov = trueModel.getSigmaObs();
    computeAOPrec();// recomputing precision
  }
  void computeRWPrec(){
    randomWalkPrec = arma::inv_sympd(randomWalkCov);
  }
  void computeAOPrec(){
    aroundObsPrec = arma::inv_sympd(aroundObsCov);
  }
  Rcpp::List computePropMoments(const arma::vec& oldX, 
                                const arma::vec& newY,
                                const double timeLag,
                                const double wDyn = 1,
                                const double wObs = 1) const{
    arma::mat diagPartMat(2, 2, arma::fill::zeros);
    diagPartMat.diag() = oldX;
    arma::mat dynamicsCov = diagPartMat * randomWalkCov * diagPartMat * timeLag;
    arma::mat dynamicsPrec = arma::inv_sympd(dynamicsCov) * wDyn;
    arma::mat obsPrec = aroundObsPrec * wObs;
    arma::mat outCov = arma::inv_sympd(dynamicsPrec + obsPrec);
    arma::colvec mDyn = trueModel.getEulerMean(oldX, timeLag);
    arma::mat outMean = outCov * (dynamicsPrec * mDyn + obsPrec * newY);
    return Rcpp::List::create(Rcpp::Named("mean") = outMean,
                              Rcpp::Named("cov")  = outCov);
  }
public:
  ProposalLVModel(Rcpp::List L)
    : trueModel(), aroundObsCov(arma::eye(2, 2)), randomWalkCov(arma::eye(2, 2)){
    arma::colvec a1(3, arma::fill::zeros);
    a1 = Rcpp::as<arma::vec>(L["a1"]);
    arma::vec a2(3);
    a2 = Rcpp::as<arma::vec>(L["a2"]);
    arma::mat gam(2, 2); gam = Rcpp::as<arma::mat>(L["gam"]);
    arma::vec mu0(2);
    mu0 = Rcpp::as<arma::vec>(L["mu0"]);
    arma::mat sigma0(2, 2); sigma0 = Rcpp::as<arma::mat>(L["sigma0"]);
    arma::mat cov(2, 2); cov = Rcpp::as<arma::mat>(L["cov"]);
    arma::vec qs(2);
    qs = Rcpp::as<arma::vec>(L["qs"]);
    arma::mat RWC(2, 2); RWC = Rcpp::as<arma::mat>(L["RWC"]);
    double wD; double wO;
    wD = L["wD"]; wO = L["wO"];
    weightDyn = wD; weightObs = wO;
    LVPOD tmp(a1, a2, mu0, sigma0, gam, cov, qs);
    trueModel = tmp;
    randomWalkCov = RWC;
    Utils::checkSDP2d(randomWalkCov, "RWC");
    computeAroundObsCov();
    computeRWPrec();
  };
  Rcpp::NumericVector evalInitialModelDensity(const arma::mat& particles) const{
    // DebugMethods::debugprint(aroundObsCov, "Covariance");
    unsigned int particleSize = particles.n_rows;
    LVModel hiddenModel = trueModel.getModel();
    arma::mat sig0 = hiddenModel.getSigma0();
    arma::colvec m0 = hiddenModel.getMu0();
    Rcpp::NumericVector out(particleSize); out.fill(0);
    for(unsigned int i = 0; i < particleSize; i++){
      arma::colvec colParticle = particles.row(i).t();
      out(i) = Utils::dmvnorm_low_trunc(colParticle, m0, sig0, lb);
    }
    return out;
  };
  arma::mat simFirstPart(unsigned int particleSize, 
                         const arma::colvec& observation) const{
    
    arma::colvec scaledObs = trueModel.getScaledObs(observation);
    LVModel hiddenModel = trueModel.getModel();
    arma::mat sig0 = hiddenModel.getSigma0();
    arma::mat prec0 = arma::inv_sympd(sig0);
    arma::colvec m0 = hiddenModel.getMu0();
    arma::mat covar = Utils::kalmanFiltVar(sig0, 
                                           aroundObsCov);
    arma::colvec mean = Utils::kalmanFiltMean(m0, scaledObs, prec0,
                                              aroundObsPrec, covar);
    // Using a log normal
    arma::colvec logNormMean = log(mean) - covar.diag() * 0.5;
    arma::mat out = Utils::rmvlognorm(particleSize, logNormMean, covar);
    // arma::mat out = Utils::rmvnorm_low_trunc(particleSize, mean, covar, lb);
    return out;
  };
  Rcpp::NumericVector evalInitialPropDensity(const arma::mat& particles, 
                                         const arma::colvec& observation) const{
    // DebugMethods::debugprint(aroundObsCov, "Covariance");
    unsigned int particleSize = particles.n_rows;
    arma::colvec scaledObs = trueModel.getScaledObs(observation);
    LVModel hiddenModel = trueModel.getModel();
    arma::mat sig0 = hiddenModel.getSigma0();
    arma::mat prec0 = arma::inv_sympd(sig0);
    arma::colvec m0 = hiddenModel.getMu0();
    arma::mat covar = Utils::kalmanFiltVar(sig0, 
                                           aroundObsCov);
    arma::colvec mean = Utils::kalmanFiltMean(m0, scaledObs, prec0,
                                              aroundObsPrec, covar);
    arma::colvec logNormMean = log(mean) - covar.diag() * 0.5;
    Rcpp::NumericVector out(particleSize); out.fill(0);
    for(unsigned int i = 0; i < particleSize; i++){
      arma::colvec colParticle = particles.row(i).t();
      // out(i) = Utils::dmvnorm_low_trunc(colParticle, mean, covar, lb);
      out(i) = Utils::dmvlognorm(colParticle, logNormMean, covar);
    }
    return out;
  };
  arma::mat simNextPart(const arma::mat& oldParticles, 
                        const arma::colvec& newObs, 
                        const double time_lag) const{
    unsigned int particleSize = oldParticles.n_rows;
    arma::mat output(oldParticles.n_cols, particleSize);
    arma::colvec scaledFutObs = trueModel.getScaledObs(newObs);
    for(unsigned int i = 0; i < particleSize; i++){
      arma::colvec oldPart       = oldParticles.row(i).t();
      Rcpp::List propMom = computePropMoments(oldPart, scaledFutObs, time_lag,
                                              weightDyn, weightObs);
      arma::colvec propMean = propMom["mean"];
      arma::mat propCov = propMom["cov"];
      arma::colvec logNormMean = log(propMean) - propCov.diag() * 0.5;
      output.col(i) = Utils::rmvlognorm(logNormMean, propCov);
      // output.col(i) = Utils::rmvnorm_low_trunc(propMean, propCov, lb);
    }  
    return output.t();
  };
  Rcpp::NumericVector densityNextParticle(const arma::mat& oldParticles,
                                          const arma::mat& newParticles,
                                          const arma::colvec newObs,
                                          double time_lag) const{
    //using euler method
    unsigned int particleSize = oldParticles.n_rows;
    Rcpp::NumericVector output(particleSize); output.fill(-1);
    arma::colvec scaledFutObs = trueModel.getScaledObs(newObs);
    for(unsigned int i = 0; i < particleSize; i++){
      arma::colvec oldPart       = oldParticles.row(i).t();
      arma::colvec newPart = newParticles.row(i).t();
      Rcpp::List propMom = computePropMoments(oldPart, scaledFutObs, time_lag,
                                              weightDyn, weightObs);
      arma::colvec propMean = propMom["mean"];
      arma::mat propCov = propMom["cov"];
      arma::colvec logNormMean = log(propMean) - propCov.diag() * 0.5;
      output(i) = Utils::dmvlognorm(newPart, logNormMean, propCov);
      // output(i) = Utils::dmvnorm_low_trunc(newPart, propMean, propCov, lb);
    }
    return output;
  }
  double evalTransitionDensityUnit(const arma::colvec oldParticle, 
                                   const arma::colvec newParticle,
                                   double time_lag,
                                   unsigned int sampleSize) const{
    return trueModel.unbiasedDensity(oldParticle, newParticle, 
                                     time_lag, sampleSize);
  };
  Rcpp::NumericVector evalTransitionDensity(const arma::mat& oldParticles,
                                            const arma::mat& newParticles,
                                            double time_lag, 
                                            unsigned int sampleSize, 
                                            unsigned int max_try = 500){
    unsigned int particleSize = oldParticles.n_rows;
    Rcpp::NumericVector output(particleSize); output.fill(0);
    bool stop_cond = false;
    unsigned int n_try = 0;
    DebugMethods db;
    while((not stop_cond) and (n_try < max_try)){
      db.here();
      n_try += 1;
      Rcpp::NumericVector cur_est(particleSize);
      for(unsigned int i = 0; i < particleSize; i++){
        arma::colvec ancester = oldParticles.row(i).t();
        arma::colvec child = newParticles.row(i).t();
        cur_est(i) = evalTransitionDensityUnit(ancester, child, 
                time_lag, sampleSize);
      }
      output = ((n_try - 1) * output + cur_est) / n_try;
      // double tol = pow(10, -10);
      // for(unsigned int i = 0; i < particleSize; i++){
      //   if(std::abs(output(i)) < tol){
      //     output(i) = 0;
      //   }
      // }
      // std::cout << "Le minimum est " << Rcpp::min(output) << std::endl;
      stop_cond = Utils::checkAllPosOr0(output);
    }
    if(not stop_cond){
      unsigned int n_neg = Utils::countNeg(output);
      std::cout << "There are " << n_neg << " negative weights" << std::endl;
      Rcpp::stop("Not all estimated transition densities are positive. Increase max_try");
    }
    return output;
  };
  Rcpp::NumericVector evalDensX0(const arma::mat& particles0){
    Rcpp::NumericVector out(particles0.n_rows);
    for(unsigned int i = 0; i < particles0.n_rows; i++){
      out(i) = trueModel.getDensX0(particles0.row(i).t());
    }
    return out;
  }
  Rcpp::NumericVector evalObsDensity(const arma::mat& newParticles,
                                     const arma::colvec& observation) const{
    return trueModel.observationDensity(newParticles, observation);
  };
  // Rcpp::NumericVector getParams() const{
  //   return trueModel.getParams();
  // }
  // void setParams(double newTheta, double newSigma){
  //   trueModel.setParams(newTheta, newSigma);
  //   updatearoundObsCov();
  // }
  // void setParams(const Rcpp::NumericVector& newParams){
  //   trueModel.setParams(newParams);
  //   updatearoundObsCov();
  // }
};

// Exposes the class to Rcpp
RCPP_MODULE(ProposalLVModel_Module) {
  using namespace Rcpp;
  class_<ProposalLVModel>("ProposalLVModel")
    // .constructor<arma::mat, arma::colvec, arma::colvec, arma::mat, arma::mat, arma::mat>("constructor")
    .constructor<Rcpp::List>("constructor")
    .method("simX0", & ProposalLVModel::simFirstPart)
    .method("simNextX", & ProposalLVModel::simNextPart)
    .method("densX0", &ProposalLVModel::evalDensX0)
    .method("transDens", &ProposalLVModel::evalTransitionDensity)
    .method("obsDens", &ProposalLVModel::evalObsDensity)
    .method("modelDens0", &ProposalLVModel::evalInitialModelDensity)
    .method("propDens0", &ProposalLVModel::evalInitialPropDensity)
    .method("propDens", &ProposalLVModel::densityNextParticle)
  ;
}
#endif