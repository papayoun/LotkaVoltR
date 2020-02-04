## ----setup, include=FALSE, message = FALSE-------------------------------
library(tidyverse)
knitr::opts_chunk$set(echo = TRUE)



initial_param <- list(a1 = a1, a2 = a2, mu0 = mu0,
                        sigma0 = Sigma0, RWC = diag(0.005, 2), 
                        qs = q_values, cov = Sigma_obs,
                        wD = 1, wO = 1,
                        gam = Gamma)
observations <- hares_lynx_data %>% 
  select(Hares, Lynx) %>% 
  as.matrix()
observation_times <- pull(hares_lynx_data, Year)
library(LotkaVoltR)
library(parallel)
mclapply(1:30, function(my_seed){
  set.seed(my_seed)
  particle_filter <- PF_create(initial_param, t(observations),
                               observation_times, 200,
                               n_euler_skel = 100)
  smoothing_exp <- particle_filter$runSmoothing()
  
  
  ## ------------------------------------------------------------------------
  particles <- particle_filter$get_particles()
  weights <- particle_filter$get_weights()
  sel_inds <- particle_filter$get_indexes() + 1
  genealogy <- matrix(NA, nrow = nrow(weights), ncol = ncol(weights))
  for(year in 1:ncol(sel_inds)){
    tmp_gen <- genealogy[sel_inds[, year], ]
    tmp_gen[,year] <- sel_inds[, year]
    genealogy <- tmp_gen
  }
  
  smoothing_traj <- array(NA, dim = dim(particles))
  for(i in 1:ncol(genealogy)){
    smoothing_traj[,, i] <- particles[genealogy[, i],, i]
  }
  
  save(smoothing_exp, particles, weights, sel_inds, smoothing_traj, genealogy,
       file = paste0("smoothing_from_th0_seed_", my_seed, ".RData"))
}, mc.cores = 10)

