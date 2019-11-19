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
  std::vector<arma::mat> tau_EX_old; // cube of dim particleSize * 2 * observationSize
  std::vector<arma::mat> tau_EX; // cube of dim particleSize * 2 * observationSize
  std::vector<Rcpp::NumericMatrix> tau_EStep;// numberModels matrices of dim particleSize * observationSize
  // // // Methods //
  Rcpp::NumericVector normWeights(const Rcpp::NumericVector& unNormedWeights) const{
    return unNormedWeights / sum(unNormedWeights);};
  // Initialization methods
  void initializeTauEStep(const unsigned int numberModels){
    std::vector<Rcpp::NumericMatrix>  initial(numberModels);
    for(unsigned int i = 0; i < numberModels; i++){
      Rcpp::NumericMatrix zeroMatrix(particleSize, observations.n_cols); 
      zeroMatrix.fill(0);
      initial[i] = zeroMatrix;
    }
    tau_EStep = initial;
  };
  void initializeTau_EX(){
    arma::mat zeroMatrix(particleSize, 2);
    for(unsigned int i = 0; i < observations.n_cols; i++){
      tau_EX_old[i] = zeroMatrix;
      tau_EX[i] = zeroMatrix;
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
  void updateTauTracking_IS(const unsigned int& childIndex, 
                            const unsigned int& childParticleIndex,
                            const unsigned int& ancestorParticleIndex,
                            const double IS_weight,
                            const bool& tracked_X){
    arma::rowvec h_term(2);
    h_term.fill(0);
    if(tracked_X){
      h_term = particleSet[ancestorParticleIndex].row(childIndex - 1);
    }
    arma::rowvec old = tau_EX_old[childIndex - 1].row(ancestorParticleIndex);
    tau_EX[childIndex].row(childParticleIndex) += IS_weight * (old + h_term);
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
//   arma::mat eval_EX_smoothing(){
//     unsigned int observationSize = observations.n_cols;
//     arma::mat output(observationSize, 2);
//     initializeTau_EX();// Initialize matrix of 0
//     setInitalParticles();
//     for(int k = 0; k < (observationSize - 1);k++){
//       propagateParticles(k);
//       // initializeBackwardSampling(k);// Samples of ancestor index is made here
//       Rcpp::NumericVector currentWeights = filteringWeights(Rcpp::_, k);
//       for(unsigned int i = 0; i < particleSize; i++){// i indexes particles
//         // setDensityUpperBound(k + 1, i);// Density upperbound for particle xi_{k+1}^i
//         sum_IS_weights = 0;
//         double curParticle = particleSet(i, k + 1);
//         // Choosing ancestoir
//         Rcpp::IntegerVector ancestInd = GenericFunctions::sampleReplace(particleIndexes,
//                                                                         backwardSampleSize,
//                                                                         currentWeights);
//         Rcpp::NumericVector ancestPart(backwardSampleSize);
//         Rcpp::NumericVector IS_weights(backwardSampleSize);
//         for(unsigned int l = 0; l < backwardSampleSize; l++){
//           ancestPart(l) = particleSet(ancestInd(l), k);
//           IS_weights(l) = propModel.evalTransitionDensityUnit(ancestPart(l), 
//                      curParticle,
//                      observationTimes(k), 
//                      observationTimes(k + 1),
//                      densitySampleSize, 
//                      false);
//           sum_IS_weights += IS_weights(l);
//         }
//         IS_weights = IS_weights / sum_IS_weights;
//         for(unsigned int l = 0; l < backwardSampleSize; l++){
//           updateTauTracking_IS(k + 1, i, ancestInd(l), IS_weights(l), update_X);
//           //k + 1 is the time index from which the backward is done, 
//           //i is the corresponding particle of this generation
//         }
//         // std::cout << "sum of IS_Weights" << sum_IS_weights << std::endl;
//       }
//     }
}; // End of the class

RCPP_MODULE(PF_Module) {
  using namespace Rcpp;
  class_<LVParticleFilter>("LVParticleFilter")
    .constructor<arma::mat, Rcpp::NumericVector, Rcpp::List, 
  unsigned int, unsigned int>("constructor")
    .method("runPF", &LVParticleFilter::runPF)
    .method("get_particles", &LVParticleFilter::getParticles)
    .method("get_weights", &LVParticleFilter::getWeights)
  ;
}


#endif