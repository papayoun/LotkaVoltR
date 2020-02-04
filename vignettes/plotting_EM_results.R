rm(list = ls())
library(tidyverse)
results_path <- "vignettes/"

pattern <- "EM_lynx_seed"
my_results <- dir(path = results_path, pattern = pattern)

all_results <- map_dfr(1:length(my_results), function(i){
  load(paste0(results_path, my_results[i]))
  a1s <- map_dfr(out, function(x){
    z <- matrix(x$a1, ncol = 3) %>% as.data.frame()
    colnames(z) = c("a[11]", "a[12]", "a[13]")
    return(z)})
  a2s <- map_dfr(out, function(x){
    z <- matrix(x$a2, ncol = 3) %>% as.data.frame()
    colnames(z) = c("a[21]", "a[22]", "a[23]")
    return(z)
  })
  gammas <- map_dfr(out, function(x){
    z <- matrix(diag(x$gam), ncol = 2) %>% as.data.frame()
    colnames(z) = c("gamma[11]", "gamma[22]")
    return(z)
  })
  sigmas <- map_dfr(out, function(x){
    z <- matrix(diag(x$cov), ncol = 2) %>% as.data.frame()
    colnames(z) = c("sigma[11]", "sigma[22]")
    return(z)
  })
  bind_cols(a1s, a2s, gammas, sigmas, iteration = 1:nrow(a1s)) %>% 
    mutate(replicate = i)
}) %>% 
  mutate(replicate = factor(replicate)) %>% 
  gather(-replicate, -iteration,
         key = "Parameter", value = "Estimate")
ggplot(all_results) +
  aes(x = iteration, y = Estimate, col = replicate) +
  facet_wrap(~Parameter, scales = "free_y", labeller = label_parsed) +
  geom_point() +
  geom_line()
