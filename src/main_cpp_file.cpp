// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "DebugMethods.h"
#include "Utils.h"
#include "LVModel.h"
#include "LVPOD.h"
#include "ProposalLVModel.h"
#include "LVParticleFilter.h"
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// ProposalLVModel generateModel(const ProposalLVModel model_ref, 
//                               double sd, double thresh = 0.01){
//   Rcpp::NumericVector a1_ref = 
//   double newTheta = theta0 + Rcpp::rnorm(1, 0, sd)[0];
//   double newSigma2 = sigma20 + Rcpp::rnorm(1, 0, sd)[0];
//   if(newSigma2 < thresh)
//     newSigma2 = thresh;
//   SINE_POD newModel(newTheta, newSigma2);
//   return newModel;
// }

//'@title Run a Particle Filter on a POD based on a continuous LV model
//'@name get_E_step
//'@param a1 a2 Prey predatorsparameters
//'@param obs matrix of observations
//'@param obsTimes_ vector of ordered obs times
//'@param cov qs RWC some parameters
//'@param n_part number of particles
//'@param n_dens_samp number of sampling for density
//'@return a list with particles and weights
//'@export
// [[Rcpp::export]]
Rcpp::NumericVector get_E_step(arma::mat obs_,
              Rcpp::NumericVector obsTimes_,//obs
              Rcpp::List myParams, // for the LV Model
              Rcpp::List testedParams, // For the tested Models
              unsigned int n_part, unsigned int n_dens_samp) {
    LVParticleFilter myPF(obs_, obsTimes_, myParams, n_part, n_dens_samp);
    int n_models = testedParams.length();
    std::vector<ProposalLVModel> testedModels(n_models); 
    for(int m = 0; m < n_models; m ++){
      Rcpp::List paramsList = testedParams[m];
      ProposalLVModel propModel(paramsList);
      testedModels[m] = propModel;
    }
    Rcpp::NumericVector output = myPF.evalEStep_IS(testedModels);
    return output;
}