## ----chosen_parameters, eval = T-----------------------------------------
rm(list = ls()) # Cleaning environment
raw_data <-  "Year	 Hares	 Lynx
1900	30	4
1901	47.2	6.1
1902	70.2	9.8
1903	77.4	35.2
1904	36.3	59.4
1905	20.6	41.7
1906	18.1	19
1907	21.4	13
1908	22	8.3
1909	25.4	9.1
1910	27.1	7.4
1911	40.3	8
1912	57	12.3
1913	76.6	19.5
1914	52.3	45.7
1915	19.5	51.1
1916	11.2	29.7
1917	7.6	15.8
1918	14.6	9.7
1919	16.2	10.1
1920	24.7	8.6"

hares_lynx_data <- read.table(text = raw_data, header = TRUE)

## Data sourced from
## http://www-rohan.sdsu.edu/~jmahaffy/courses/f00/math122/labs/labj/q3v1.htm
## Originally published in E. P. Odum (1953), Fundamentals of Ecology,
## Philadelphia, W. B. Saunders.
## Initial values inspired from https://gist.github.com/mages/1f0f0d5bbe50af81cc19

# Dynamics parameters

a1 <- c(0.695617134, 0.001147425, 0.030333635) # Prey parameters
a2 <- c(0.630905980, 0.021756719, 0.001352425) # Predator parameters
Gamma <- matrix(c(0.22, 0, 0, 0.2), nrow = 2) # Diffusion parameters
mu0 <- c(30, 4) # Mean of the inital distribution
Sigma0 <- diag(1, 2) # Variance of the initial distribution

# Observation process parameters

Sigma_obs <- matrix(c(0.01, 0.00, 0.00, 0.01), ncol = 2) # Observation noise
q_values <- c(1, 1) # Known q values
initial_param <- list(a1 = a1, a2 = a2, mu0 = mu0,
                      sigma0 = Sigma0, RWC = diag(0.005, 2), 
                      qs = q_values, cov = Sigma_obs,
                      wD = 1, wO = 1,
                      gam = Gamma)

# Model creation

library(LotkaVoltR) # Dedicated library
# Creation of a model instance in R
initial_param <- list(a1 = a1, a2 = a2, mu0 = mu0,
                      sigma0 = Sigma0, RWC = diag(0.005, 2), 
                      qs = q_values, cov = Sigma_obs,
                      wD = 1, wO = 1,
                      gam = Gamma)

# Function generation candidat --------------------------------------------


generate_candidate_list <- function(par_list,
                                    standard_dev_list,
                                    up_lim_list){
  out <- par_list
  # Generation a1, dans [0, 1] * [0, 0.05] * [0, 0.05]
  sigmoid <- function(x){
    1 / (1 + exp(-x))
  }
  logit <- function(x){
    log(x / (1 - x))
  }
  # Version dégeu (version tidy en bas)
  out$a1 <- sigmoid(rnorm(n = 3,
                          mean = logit(par_list$a1 / up_lim_list$a1),
                          sd = standard_dev_list$a1)) * up_lim_list$a1
  # out$a1 <- (par_list$a1 / up_lim_list$a1) %>% # Retour dans [0, 1]
  #   logit() %>% # Retour dans le monde réel
  #   rnorm(n = 3, sd = standard_dev_list$a1) %>% # Déplacement 
  #   sigmoid() %>% # Retour dans [0,1]
  #   {. * up_lim_list$a1} # Retour dans l'espace contraint
  out$a2 <- sigmoid(rnorm(n = 3,
                          mean = logit(par_list$a2 / up_lim_list$a2),
                          sd = standard_dev_list$a2)) * up_lim_list$a2
  # out$a2 <- (par_list$a2 / up_lim_list$a2) %>% # Retour dans [0, 1]
  #   logit() %>% # Retour dans le monde réel
  #   rnorm(n = 3, sd = standard_dev_list$a2) %>% # Déplacement 
  #   sigmoid() %>% # Retour dans [0,1]
  #   {. * up_lim_list$a2} # Retour dans [0, 1] * [0, 0.1] * [0, 0.1]
  diag(out$gam) <- sigmoid(rnorm(n = 2,
                                 mean = logit(diag(par_list$gam) / up_lim_list$gam),
                                 sd = standard_dev_list$gam)) * up_lim_list$gam
  # diag(out$gam) <- (diag(par_list$gam) / up_lim_list$gam) %>% # Retour dans [0, 1]
  #   logit() %>% # Retour dans le monde réel
  #   rnorm(n = 2, sd = standard_dev_list$gam) %>% # Déplacement 
  #   sigmoid() %>% # Retour dans [0,1]
  #   {. * up_lim_list$gam}
  diag(out$cov) <-  sigmoid(rnorm(n = 2,
                                  mean = logit(diag(par_list$cov) / up_lim_list$cov),
                                  sd = standard_dev_list$cov)) * up_lim_list$cov
  # diag(out$cov) <- (diag(par_list$cov) / up_lim_list$cov) %>% # Retour dans [0, 1]
  #   logit() %>% # Retour dans le monde réel
  #   rnorm(n = 2, sd = standard_dev_list$cov) %>% # Déplacement 
  #   sigmoid() %>% # Retour dans [0,1]
  #   {. * up_lim_list$cov}
  return(out)
}


EM_function <- function(obs, obs_times, initial_param, initial_sd, up_lims,
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
                                                       initial_sd,
                                                       up_lims)),
                  list(param_0))
    E_step_evals <- get_E_step(obs_ = t(obs), 
                               obsTimes_ = obs_times, 
                               myParams = param_0, 
                               testedParams = par_list, 
                               n_part = 200, n_dens_samp = 100)
    initial_sd <- purrr::map(initial_sd,
                             function(x) 0.9 * x)
    param_0 <- par_list[[which.max(E_step_evals)]]
    out <- c(out, list(param_0))
    final_E_steps[i, ] <- E_step_evals 
    save(out, final_E_steps, initial_sd, 
         file = paste0("EM_", name, "_seed_", seed, ".RData"))
  }
  return(out)
}

observations <- as.matrix(hares_lynx_data[, c("Hares", "Lynx")])
observation_times <- hares_lynx_data$Year


sd_list <- list(a1 = c(1, 1, 1),
                a2 = c(1, 1, 1),
                cov = 0.5,
                gam = 0.5)
upper_lims_list <- list(a1 = c(1, 0.01, 0.05),
                        a2 = c(1, 0.05, 0.01),
                        cov = 0.5,
                        gam = 0.5)
library(parallel)
my_seed <- sample(1:10000, size = 1)
EM_function(obs = observations, 
            obs_times = observation_times,
            initial_param = initial_param, 
            initial_sd = sd_list,
            up_lims = upper_lims_list,
            n_cands = 10,
            n_iter = 30, seed = my_seed,
            name = "lynx")

