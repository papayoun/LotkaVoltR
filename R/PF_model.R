loadModule("PF_Module", TRUE, loadNow = F)

#' @title Create a particle filter
#' @name PF_create
#' @description  function that creates a LVPOD instance
#' @param L a list of parameters containing a1, a2, gam, cov, qs, RWC, wO and wD
#' @param observations a matrix of observations
#' @param obs_times a vector of strictly increasing times
#' @param n_particles number of particles
#' @param n_euler_skel number of skeleton (Monte Carlo effort) for CIS
#' @return a ParticleFilter
#' @export
PF_create <- function(L = list(a1 = c(0.5, 0.1, 0.2),
                               a2 = c(0.1, 0.4, 0.1),
                               gam = matrix(c(1, -0.3,
                                              -0.3, 1), nrow = 2),
                               mu0 = c(50, 20),
                               sigma0 = diag(1, 2),
                               cov = diag(1, 2),
                               qs = c(0.2, 0.3),
                               RWC = diag(10, 2),
                               wO = 1, wD = 1),
                      observations, obs_times, n_particles, n_euler_skel){
  model <- new(LVParticleFilter, observations, obs_times, L, n_particles,
               n_euler_skel)
  model
}