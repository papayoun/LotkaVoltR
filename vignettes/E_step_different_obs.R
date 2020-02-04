## ----chosen_parameters, eval = T-----------------------------------------
rm(list = ls()) # Cleaning environment

# Dynamics parameters

a1 <- c(12, 0.05, 1) # Prey parameters
a2 <- c(2, 0.2, 0.1) # Predator parameters
Gamma <- matrix(c(0.5, 0.1, 0.1, 0.2), nrow = 2) # Diffusion parameters
mu0 <- c(20, 10) # Mean of the inital distribution
Sigma0 <- diag(1, 2) # Variance of the initial distribution

# Observation process parameters

Sigma_obs <- matrix(c(0.01, 0.005, 0.005, 0.01), ncol = 2) # Observation noise
q_values <- c(0.2, 0.3) # Known q values

# Model creation

library(LotkaVoltR) # Dedicated library
# Creation of a model instance in R
POLV_model <- POLV_create(a1 = a1, a2 = a2, gam = Gamma, mu0 = mu0, 
                          sigma0 = Sigma0, cov = Sigma_obs, qs = q_values)
parameters_list <- list(a1 = a1, a2 = a2, mu0 = mu0,
                        sigma0 = Sigma0, RWC = diag(0.005, 2), 
                        qs = q_values, cov = Sigma_obs,
                        wD = 1, wO = 1,
                        gam = Gamma)
# Simulation parameters

simulation_times <- seq(from = 0, by = 1e-6, # Simulation time step, small!!
                        length.out = 10^6 + 1)

# If the simulation must be done at small time steps, the output can be thinned
selection <- seq(1, length(simulation_times), length.out = 201)
my_inference_function <- function(my_seed, my_name){
  set.seed(my_seed)
  simulated_process <- POLV_simulate(POLV_model,
                                     times = simulation_times, 
                                     selection = selection)
  observation_times <- simulation_times[selection]
  observations <- simulated_process[, c("Y1", "Y2")]
  
  
  ## ------------------------------------------------------------------------
  true_param <- list(a1 = a1, a2 = a2, mu0 = mu0,
                     sigma0 = Sigma0, RWC = diag(0.005, 2), 
                     qs = q_values, cov = Sigma_obs,
                     wD = 1, wO = 1,
                     gam = Gamma)
  sd_list <- list(a1 = c(2, 0.05, 0.5),
                  a2 = c(1, 0.25, 0.05),
                  cov = 0,
                  gam = 0)
  if(my_name == "all"){
    sd_list$cov = c(0.001, 0.0005, 0.001)
    sd_list$gam = c(0.1, 0.05, 0.1)
    print("Go for all")
  }
  else if(my_name == "only_a"){
    print("only a is estimated")
  }
  else{
    stop("unknown name")
  }
  
  ## ----functions_par_candidate_generation----------------------------------
  def_pos <- function(mat){
    cond1 <- all(diag(mat) > 0)
    cond2 <- det(mat) > 0
    cond1 & cond2
  }
  generate_candidate_list <- function(par_list,
                                      standard_dev_list){
    out <- par_list
    out$a1 <- abs(par_list$a1 + standard_dev_list$a1 * rnorm(3))
    out$a2 <- abs(par_list$a2 + standard_dev_list$a2 * rnorm(3))
    cond_arret_cov <- FALSE
    while(!cond_arret_cov){
      innov_cov <- standard_dev_list$cov * rnorm(3)
      out$cov <- par_list$cov + 
        matrix(c(innov_cov[1], innov_cov[2], 
                 innov_cov[2], innov_cov[3]), nrow = 2)
      cond_arret_cov <- def_pos(out$cov)
    }
    cond_arret_gamma <- FALSE
    while(!cond_arret_gamma){
      innov_gamma <- standard_dev_list$gam * rnorm(3)
      out$gam <- par_list$gam + 
        matrix(c(innov_gamma[1], innov_gamma[2], 
                 innov_gamma[2], innov_gamma[3]), nrow = 2)
      cond_arret_gamma <- def_pos(out$gam)
    }
    return(out)
  }
  
  ## ------------------------------------------------------------------------
  
  EM_function <- function(obs, obs_times, initial_param, initial_sd,
                          n_cands, n_iter,
                          seed,
                          name = NULL){
    param_0 <- initial_param
    out <- list(param_0)
    set.seed(seed)
    final_E_steps <- matrix(NA, nrow = n_iter, ncol = n_cands + 1)
    for(i in 1:n_iter){
      print(paste("Iteration", i))
      par_list <- c(purrr::rerun(n_cands, 
                                 generate_candidate_list(param_0,
                                                         initial_sd)),
                    list(param_0))
      E_step_evals <- get_E_step(obs_ = t(obs), 
                                 obsTimes_ = obs_times, 
                                 myParams = param_0, 
                                 testedParams = par_list, 
                                 n_part = 100, n_dens_samp = 30)
      initial_sd <- purrr::map(initial_sd,
                               function(x) 0.9 * x)
      param_0 <- par_list[[which.max(E_step_evals)]]
      out <- c(out, list(param_0))
      final_E_steps[i, ] <- E_step_evals 
      save(out, final_E_steps, initial_sd, 
           file = paste0("EM_repl_", seed, "_", name, "_diff_obs.RData"))
    }
    return(out)
  }
  param0 <- true_param
  EM_function(obs = observations, 
              obs_times = observation_times,
              initial_param = param0, initial_sd = sd_list, n_cands = 10,
              n_iter = 30, seed = my_seed,
              name = my_name)
}

