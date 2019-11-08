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


//'@title Run a Particle Filter on a POD based on a continuous LV model
//'@name runLV_PF
//'@param a1 a2 Prey predatorsparameters
//'@param obs matrix of observations
//'@param obsTimes_ vector of ordered obs times
//'@param cov qs RWC some parameters
//'@param n_part number of particles
//'@param n_dens_samp number of sampling for density
//'@return a list with particles and weights
//'@export
// [[Rcpp::export]]
void runLV_PF(arma::mat obs_, Rcpp::NumericVector obsTimes_,//obs
                       Rcpp::List myParams, // for the LV Model
                       unsigned int n_part, unsigned int n_dens_samp) {
    // LVParticleFilter myPF(obs_, obsTimes_, myParams, n_part, n_dens_samp);
    // myPF.runPF();
    // return Rcpp::List::create(Rcpp::Named("parts") = myPF.getParticles(),
    //                           Rcpp::Named("weights") = myPF.getWeights());
}