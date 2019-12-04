#include <RcppArmadillo.h>      // declarations for both Rcpp and RcppArmadillo offering Armadillo classes
#ifndef ProposalLVModel_H
#define ProposalLVModel_H

//' @title Class ParticleFilter
//' @name LV_PF
//' @description Class to run a particle filter
//' @export LVParticleFilter
class LVParticleFilter{
private:
  // // // Attributes // // //
  // Observations objects
  const arma::mat observations; //2 * n
  const Rcpp::NumericVector observationTimes; // size n
  //Model objects
  ProposalLVModel propModel;
  // Filtering objects
  const unsigned int particleSize; // Number of particles
  std::vector<arma::mat> particleSet;
  Rcpp::NumericMatrix filteringWeights; // particleSize * n matrix
  const Rcpp::IntegerVector particleInd; // vector 1:particleSize
  const unsigned int densitySampleSize;
  // backward sampling objects
  unsigned int backwardSamplingCounter; // To keep track of how many tries were done
  Rcpp::IntegerVector  backwardIndexCandidates; // simulate several at a time
  // Inference objects
  const Rcpp::IntegerVector particleIndexes;
  const unsigned int backwardSampleSize = 30;
  std::vector<arma::mat> tau_EX_old; // cube of dim particleSize * 2 * observationSize
  std::vector<arma::mat> tau_EX; // cube of dim particleSize * 2 * observationSize
  std::vector<arma::mat> tau_EStep;// numberModels matrices of dim particleSize * observationSize
  // // // Methods //
  Rcpp::NumericVector normWeights(const Rcpp::NumericVector& unNormedWeights) const{
    return unNormedWeights / sum(unNormedWeights);};
  // Initialization methods
  void initializeTauEStep(const unsigned int numberModels){
    std::vector<arma::mat>  initial(numberModels);
    arma::mat zeroMatrix(particleSize, 2);
    zeroMatrix.fill(0);
    for(unsigned int i = 0; i < numberModels; i++){
      initial[i] = zeroMatrix;
    }
    tau_EStep = initial;
  };
  void initializeTau_EX(bool both = false){
    arma::mat zeroMatrix(particleSize, 2);
    zeroMatrix.fill(0);
    for(unsigned int i = 0; i < observations.n_cols; i++){
      tau_EX[i] = zeroMatrix;
    }
    if(both){
      for(unsigned int i = 0; i < observations.n_cols; i++){
        tau_EX_old[i] = zeroMatrix;
      }
    }
  };
  void initializeParticleSet(){
    for(int i = 0; i < observations.n_cols; i++){
      arma::mat zeroMatrix(particleSize, 2); 
      zeroMatrix.fill(0.0);
      particleSet[i] = zeroMatrix;
    }
  };
  void setInitalParticles(){ // Simulating initial particles
    arma::mat firstPart = propModel.simFirstPart(particleSize, 
                                                  observations.col(0));
    particleSet[0] = firstPart;
    Rcpp::NumericVector obsDens = propModel.evalObsDensity(firstPart,
                                       observations.col(0));
    Rcpp::NumericVector modDens = propModel.evalInitialModelDensity(firstPart);
    Rcpp::NumericVector propDens = propModel.evalInitialPropDensity(firstPart,
                                                                    observations.col(0));
    Rcpp::NumericVector unNormedWeights = obsDens * modDens / propDens;
    filteringWeights(Rcpp::_, 0) = normWeights(unNormedWeights);
  };
  // Propagation method
  void updateTauTracking_IS(const unsigned int& timeStampIndex, 
                            const unsigned int& childParticleIndex,
                            const unsigned int& ancestorParticleIndex,
                            const double IS_weight,
                            const bool& tracked_X){
    arma::rowvec h_term(2);
    h_term.fill(0);
    if(tracked_X){
      h_term = particleSet[timeStampIndex].row(ancestorParticleIndex);
    }
    arma::rowvec old = tau_EX_old[timeStampIndex].row(ancestorParticleIndex);
    tau_EX[timeStampIndex].row(childParticleIndex) += IS_weight * (old + h_term);
  };
  void updateTauEStep_IS(const unsigned int& childIndex, 
                         const unsigned int& childParticleIndex,
                         const unsigned int& ancestorParticleIndex,
                         const double IS_weight,
                         const std::vector<ProposalLVModel>& testedModels){
    double time_lag = observationTimes(childIndex) - observationTimes(childIndex - 1);
    arma::rowvec old_particle = particleSet[childIndex - 1].row(ancestorParticleIndex);
    arma::rowvec new_particle = particleSet[childIndex].row(childParticleIndex);
    for(unsigned int m = 0; m < testedModels.size(); m++){
     ProposalLVModel model = testedModels[m];
      double sampledLogQ = model.getModel().unbiasedLogDensityEstimate(particleSet(ancestorParticleIndex, childIndex - 1),
                                          ,
                                          
                                          logDensitySampleSize,
                                          skeletonSimulationMaxTry);
      double transDens = model.evalTransitionDensity(selectedParticles, 
                                                         newParts,
                                                         time_lag, 
                                                         densitySampleSize)[0];
      double logtransDens = - pow(10, 9);
      if(transDens > 0){
        logtransDens = log(transDens);
      }
      double logObsDensityTerm =  log(model.observationDensity(particleSet(childParticleIndex, childIndex),
                                                               observations(childIndex)));
      tauEStep[m](childParticleIndex, childIndex) += IS_weight * (tauEStep[m](ancestorParticleIndex, childIndex - 1) +
        sampledLogQ + logObsDensityTerm);
    }
  };
  void propagateParticles(const unsigned int& ancestorIndex){// Particle propagation
    double time_lag = (observationTimes[ancestorIndex + 1] -
                        observationTimes[ancestorIndex]);
    Rcpp::NumericVector currentWeights = filteringWeights(Rcpp::_, ancestorIndex);
    Rcpp::IntegerVector selectedInd = Utils::sampleReplace(particleInd, 
                                                           particleSize,
                                                           currentWeights);
    arma::mat selectedParticles(particleSet[ancestorIndex]);
    for(unsigned int i =0; i < particleSize; i++){
      selectedParticles.row(i) = particleSet[ancestorIndex].row(selectedInd(i));
    }
    arma::colvec futureObs = observations.col(ancestorIndex + 1);
    arma::mat newParts = propModel.simNextPart(selectedParticles,
                                               futureObs, time_lag);
    particleSet[ancestorIndex + 1] = newParts;
    Rcpp::NumericVector propDens = propModel.densityNextParticle(selectedParticles, 
                                                                 newParts, 
                                                                 futureObs, 
                                                                 time_lag);
    Rcpp::NumericVector transDens = propModel.evalTransitionDensity(selectedParticles, 
                                                                    newParts,
                                                                    time_lag, 
                                                                    densitySampleSize);// the true is to use GPE2
    Rcpp::NumericVector obsDens = propModel.evalObsDensity(newParts, futureObs);
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
      filteringWeights(n_part, obsTimes_.size()),
      particleInd(Rcpp::seq_len(n_part)  - 1),
      tau_EX_old(obs_.n_cols),
      tau_EX(obs_.n_cols),
      particleSet(obs_.n_cols),
      particleIndexes(Rcpp::seq_len(n_part)  - 1),
      densitySampleSize(n_dens_samp){
    initializeParticleSet();
  };
  // Getters
  arma::cube getParticles() const{
    arma::cube output(particleSize, 2, observations.n_cols);
    for(int i = 0; i < observations.n_cols; i ++){
      output.slice(i) = particleSet[i];
    }
    return output;
  };
  arma::cube get_tau_EX() const{
    arma::cube output(particleSize, 2, observations.n_cols);
    for(int i = 0; i < observations.n_cols; i ++){
      output.slice(i) = tau_EX[i];
    }
    return output;
  };
  Rcpp::NumericMatrix getWeights() const{
    return filteringWeights;
  };
  // Main function
  void runPF(){
    // DebugMethods db;
    // db.here();
    setInitalParticles();
    // db.here();
    for(unsigned int k = 0; k < (observationTimes.size() - 1); k++){
      // db.here();
      propagateParticles(k);
    }
  };
  Rcpp::NumericMatrix eval_EX_smoothing(){
    unsigned int observationSize = observations.n_cols;
    initializeTau_EX(true);// Initialize both tau_EX and tau_EX_old to 0 arrays
    setInitalParticles();
    DebugMethods db;
    for(int k = 0; k < (observationSize - 1); k++){
      // DebugMethods db;
      std::cout << "observation " << k <<std::endl;
      propagateParticles(k);
      // initializeBackwardSampling(k);// Samples of ancestor index is made here
      Rcpp::NumericVector currentWeights = filteringWeights(Rcpp::_, k);
      for(unsigned int i = 0; i < particleSize; i++){// i indexes particles
        // db.here();
        // DebugMethods::debugprint(tau_EX[0], "EX");
        // DebugMethods::debugprint(tau_EX_old[0], "EX_old");
        // setDensityUpperBound(k + 1, i);// Density upperbound for particle xi_{k+1}^i
        double sum_IS_weights = 0;
        arma::mat curParticle = particleSet[k + 1].row(i);
        // Choosing ancestoir
        Rcpp::IntegerVector ancestInd = Utils::sampleReplace(particleIndexes,
                                                             backwardSampleSize,
                                                             currentWeights);
        Rcpp::NumericVector IS_weights(backwardSampleSize);
        for(unsigned int l = 0; l < backwardSampleSize; l++){
          arma::mat ancParticle = particleSet[k].row(ancestInd(l));
          IS_weights(l) = propModel.evalTransitionDensity(ancParticle,
                     curParticle,
                     observationTimes(k + 1) - observationTimes(k),
                     densitySampleSize, 200)(0);
          sum_IS_weights += IS_weights(l);
        }
        IS_weights = IS_weights / sum_IS_weights;
        for(unsigned int l = 0; l < backwardSampleSize; l++){
          for(int j = 0; j < k + 1; j++){
            bool get_h = (j == k); // In this cas, the h function is either x or 0
            updateTauTracking_IS(j, i, ancestInd(l), IS_weights(l), get_h);
          }
        }
      }
      tau_EX_old = tau_EX; // Updating the functionnal tau_{k-1} <- tau_k
      initializeTau_EX(false); // reputing tau_EX to 0
    }
    tau_EX_old[observationSize - 1] = particleSet[observationSize - 1]; 
    Rcpp::NumericVector lastWeights = filteringWeights(Rcpp::_, 
                                                       observationSize - 1);
    Rcpp::NumericMatrix output(observationSize, 2);
    output.fill(0);
    for(int k = 0; k < observationSize; k++){
      arma::mat currentTau = tau_EX_old[k];
      for(int i = 0; i < particleSize; i ++){
        output(k, 0) += lastWeights(i) * currentTau(i, 0);
        output(k, 1) += lastWeights(i) * currentTau(i, 1);
      }
    }
    return output;
}// End of eval_EX_smoothing
  Rcpp::NumericVector evalEStep_IS(const std::vector<ProposalLVModel>& testedModels){
    unsigned int numberModels = testedModels.size();
    Rcpp::NumericVector output(numberModels);
    initializeTauEStep(numberModels);// Initialize matrix of 0
    setInitalParticles();
    unsigned int observationSize = observations.n_cols;
    for(int k = 0; k < (observationSize - 1); k++){
      // DebugMethods db;
      std::cout << "observation " << k <<std::endl;
      propagateParticles(k);
      // initializeBackwardSampling(k);// Samples of ancestor index is made here
      Rcpp::NumericVector currentWeights = filteringWeights(Rcpp::_, k);
      for(unsigned int i = 0; i < particleSize; i++){// i indexes particles
        // db.here();
        // DebugMethods::debugprint(tau_EX[0], "EX");
        // DebugMethods::debugprint(tau_EX_old[0], "EX_old");
        // setDensityUpperBound(k + 1, i);// Density upperbound for particle xi_{k+1}^i
        double sum_IS_weights = 0;
        arma::mat curParticle = particleSet[k + 1].row(i);
        // Choosing ancestor
        Rcpp::IntegerVector ancestInd = Utils::sampleReplace(particleIndexes,
                                                             backwardSampleSize,
                                                             currentWeights);
        Rcpp::NumericVector IS_weights(backwardSampleSize);
        for(unsigned int l = 0; l < backwardSampleSize; l++){
          arma::mat ancParticle = particleSet[k].row(ancestInd(l));
          IS_weights(l) = propModel.evalTransitionDensity(ancParticle,
                     curParticle,
                     observationTimes(k + 1) - observationTimes(k),
                     densitySampleSize, 200)(0);
          sum_IS_weights += IS_weights(l);
        }
        IS_weights = IS_weights / sum_IS_weights;
        for(unsigned int l = 0; l < backwardSampleSize; l++){
          updateTauEStep_IS(k + 1, i, ancestInd(l), IS_weights(l), testedModels);
          //k + 1 is the time index from which the backward is done, 
          //i is the corresponding particle of this generation
        }
      }
    }
    Rcpp::NumericVector lastWeights = filteringWeights(Rcpp::_, observationSize - 1);
    for(int m = 0; m < numberModels; m++){
      output[m] = sum(lastWeights * tau_EStep[m](Rcpp::_, observationSize - 1));
    }
    return output;
  };// end of evalEstep method;
}; // End of the class

RCPP_MODULE(PF_Module) {
  using namespace Rcpp;
  class_<LVParticleFilter>("LVParticleFilter")
    .constructor<arma::mat, Rcpp::NumericVector, Rcpp::List, 
  unsigned int, unsigned int>("constructor")
    .method("runPF", &LVParticleFilter::runPF)
    .method("get_particles", &LVParticleFilter::getParticles)
    .method("get_weights", &LVParticleFilter::getWeights)
    .method("runSmoothing", &LVParticleFilter::eval_EX_smoothing)
    .method("get_tau_EX", &LVParticleFilter::get_tau_EX)
  ;
}


#endif