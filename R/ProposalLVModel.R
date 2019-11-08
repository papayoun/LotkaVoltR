loadModule("ProposalLVModel_Module", TRUE)

#' @title Create
#' @name PLV_create
#' @description  function that creates a LVPOD instance
#' @param a1 parameters for the prey dynamics
#' @param a2 parameters for the predator dynamics
#' @param gam Covariance parameters for the dynamics
#' @param cov Covariance parameters for the observation model
#' @param qs Capturability parameter
#' @param RWC random walk covariance
#' @return a Proposal LVPOD
#' @export
PLV_create <- function(L = list(a1 = c(0.5, 0.1, 0.2),
                                a2 = c(0.1, 0.4, 0.1),
                                gam = matrix(c(1, -0.3,
                                               -0.3, 1), nrow = 2),
                                mu0 = c(50, 20),
                                sigma0 = diag(1, 2),
                                cov = diag(1, 2),
                                qs = c(0.2, 0.3),
                                RWC = diag(10, 2))){
  model <- new(ProposalLVModel, L)
  model
}

#' @title Simulate initial particle
#' @name PLV_simX0
#' @description  function that creates a LVPOD instance
#' @param PLV PLV model
#' @param Obs First observation
#' @param N Number of particles
#' @return a Proposal LVPOD
#' @export

PLV_simX0 <- function(PLV, Obs, N){
  PLV$simX0(N, Obs)
}

#' @title GetWeights of a Particle
#' @name PLV_getW0
#' @description  function that creates a LVPOD instance
#' @param PLV PLV model
#' @param Obs First observation
#' @return a set of weights
#' @export

PLV_getW0 <- function(PLV, part0, Obs){
  oD <- PLV$obsDens(part0, Obs)
  mD <- PLV$modelDens0(part0)
  pD <- PLV$propDens0(part0, Obs)
  uW <- oD * mD / pD
  cbind(weight = uW / sum(uW), oD = oD, pD = pD, mD = mD)
}
