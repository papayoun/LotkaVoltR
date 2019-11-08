loadModule("LV_Module", what = T, loadNow = F)

#' @title Create
#' @name LV_create
#' @description  function that creates a LVPOD instance
#' @param a1 parameters for the prey dynamics
#' @param a2 parameters for the predator dynamics
#' @param gam Covariance parameters for the dynamics
#' @param mu0 mean
#' @param sigma0 cov
#' @return a LVPOD
#' @export
LV_create <- function(a1 = c(12, 0.05, 1), a2 = c(2, 0.2, 0.1), 
                      mu0 = c(50, 20),
                      gam = matrix(c(0.2, -0.1, -0.1, 0.2), nrow = 2),
                      sigma0 = diag(1, 2)){
  model <- new(LVModel, a1, a2, mu0, sigma0, gam)
  model
}
