loadModule("LVPOD_Module", TRUE)

#' @title Create
#' @name LV_create
#' @description  function that creates a LVPOD instance
#' @param a1 parameters for the prey dynamics
#' @param a2 parameters for the predator dynamics
#' @param gam Covariance parameters for the dynamics
#' @param cov Covariance parameters for the observation model
#' @param qs Capturability parameter
#' @return a LVPOD
#' @export
LV_create <- function(a1 = c(0.5, 0.1, 0.2),
                      a2 = c(0.1, 0.4, 0.1),
                      mu0 = c(50, 20),
                      sigma0 = diag(1, 2),
                      gam = matrix(c(1, -0.3,
                                     -0.3, 1), nrow = 2),
                      cov = diag(1, 2),
                      qs = c(0.2, 0.3)){
  model <- new(LVPOD, a1, a2, mu0, sigma0, gam, cov, qs)
  model
}

#' @title simulate
#' @name LV_simulate
#' @description  function simulating from a Lotka Volterra model
#' @param LV_model an object of class LVModel
#' @param simTimes simulation times. Should be with low increment so that the scheme
#' @param selection selected points, default to all
#' is valid
#' @return a LVModel
#' @export
LV_simulate <- function(LV_model,
                        times = seq(0, 10, length.out = 2 * 10^4),
                        selection = T){
  all_simu <- LV_model$simulate(times)
  sel_sim <- do.call(cbind, 
                     lapply(all_simu,  function(x){
                       if(is.null(dim(x)))
                         return(x[selection])
                       else
                         return(x[selection,])
                     }))
  colnames(sel_sim) <- c("Y1", "Y2", "X1", "X2", "t")
  sel_sim
}

#' @title Unbiased density for a given Lotka Volterra between two points
#' @name LV_density
#' @description  function performing CIS
#' @param LV_model an object of class LVModel
#' @param x0 starting 2d point
#' @param xT ending 2d point
#' @param time_lag time_lag between observations
#' @return a (signed) estimate of the density
#' @export
LV_density <- function(LV_model, x0, xT, time_lag, n_samples = 50){
  LV_model$density(x0, xT, time_lag, n_samples)
}