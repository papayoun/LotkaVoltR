#include <RcppArmadillo.h>      // declarations for both Rcpp and RcppArmadillo offering Armadillo classes
#ifndef ProposalLVModel_H
#define ProposalLVModel_H

class LVParticleFilter{
private:
  // // // Attributes // // //
  // Observations objects
  const arma::mat observations; //2 * n
  const Rcpp::NumericVector observationTimes;
  //Model objects
  ProposalLVModel propModel;
  // Filtering objects
  arma::cube particleSet;
  Rcpp::NumericMatrix filteringWeights;
  const unsigned int particleSize;
  const Rcpp::IntegerVector particleInd; // vector 1:particleSize
  const unsigned int densitySampleSize;
  Rcpp::NumericVector obsDens;
  // // // Methods //
  Rcpp::NumericVector normWeights(const Rcpp::NumericVector& unNormedWeights) const{
    return unNormedWeights / sum(unNormedWeights);};
  void setInitalParticles(){
    arma::mat firstPart = propModel.simFirstPart(particleSize, 
                                                  observations.col(0));
    particleSet.slice(0) = firstPart;
    Rcpp::NumericVector obsDens = propModel.evalObsDensity(firstPart,
                                                           observations.col(0));
    Rcpp::NumericVector modDens = propModel.evalInitialModelDensity(firstPart);
    Rcpp::NumericVector propDens = propModel.evalInitialPropDensity(firstPart,
                                                                    observations.col(0));
    Rcpp::NumericVector unNormedWeights = obsDens * modDens / propDens;
    filteringWeights(Rcpp::_, 0) = normWeights(unNormedWeights);
  };
  // Propagation method
  void propagateParticles(const unsigned int& ancestorIndex){
    double time_lag = (observationTimes[ancestorIndex + 1] -
                        observationTimes[ancestorIndex]);
    Rcpp::NumericVector currentWeights = filteringWeights(Rcpp::_, ancestorIndex);
    Rcpp::IntegerVector selectedInd = Utils::sampleReplace(particleInd, 
                                                           particleSize,
                                                           currentWeights);
    arma::mat selectedParticles(particleSet.slice(ancestorIndex));
    for(unsigned int i =0; i < particleSize; i++){
      selectedParticles.row(i) = particleSet.slice(ancestorIndex).row(selectedInd(i));
    }
    arma::colvec futureObs = observations.col(ancestorIndex + 1);
    arma::mat newParts = propModel.simNextPart(selectedParticles,
                                               futureObs, time_lag);
    particleSet.slice(ancestorIndex + 1) = newParts;
    Rcpp::NumericVector propDens = propModel.densityNextParticle(selectedParticles, 
                                                                 newParts, 
                                                                 futureObs, 
                                                                 time_lag);
    Rcpp::NumericVector transDens = propModel.evalTransitionDensity(selectedParticles, 
                                                                    newParts,
                                                                    time_lag, 
                                                                    densitySampleSize);// the true is to use GPE2
    obsDens = propModel.evalObsDensity(newParts, futureObs);
    Rcpp::NumericVector unNormedWeights = transDens * obsDens / propDens;
    filteringWeights(Rcpp::_, ancestorIndex + 1) = normWeights(unNormedWeights);
  };
public:
  // Constructor
  LVParticleFilter(const arma::mat& obs_, 
                   const Rcpp::NumericVector& obsTimes_,//obs
                   const Rcpp::List& params, // for the LV Model
                   unsigned int n_part, unsigned int n_dens_samp)
    : observations(obs_), observationTimes(obsTimes_),
      propModel(params), particleSize(n_part),
      particleSet(n_part, 2, obsTimes_.size()),
      filteringWeights(n_part, obsTimes_.size()),
      particleInd(Rcpp::seq_len(n_part)  - 1),
      obsDens(n_part), densitySampleSize(n_dens_samp){};
  // Getters
  arma::cube getParticles() const{return particleSet;};
  Rcpp::NumericMatrix getWeights() const{return filteringWeights;};
  void runPF(){
    setInitalParticles();
    for(unsigned int k = 0; k < (observationTimes.size() - 1); k++){
      propagateParticles(k);
    }
  }
};

#endif