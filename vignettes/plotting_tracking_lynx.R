rm(list = ls())
library(tidyverse)
source("vignettes/lynx_data_parameters.R")
results_path <- "vignettes/"

pattern <- "smoothing_from_th0_seed_"
my_results <- dir(path = results_path, pattern = pattern)
my_seeds <- my_results %>% 
  str_extract_all("[0-9]+") %>% # On extrait tous les numeros
  map_chr(2) # On prend le 2e
all_results <- map(1:length(my_results), function(i){
  load(paste0(results_path, my_results[i]))
  n_particles <- nrow(weights)
  n_obs <- ncol(weights)
  PM_smooth <- map_dfr(1:n_obs,
                       function(k){
                         smooth_parts <- smoothing_traj[,, k]
                         tibble(Year = 1899 + k,
                                particle_index  = 1:n_particles,
                                Hares = smooth_parts[, 1],
                                Lynx = smooth_parts[, 2],
                                weight = weights[, n_obs],
                                seed = my_seeds[i],
                                type = "PM_smooth")
                       })
  Paris_smooth <- tibble(Year = hares_lynx_data$Year[1:20],
                         seed = my_seeds[i],
                         type = "Paris_smooth",
                         Hares = smoothing_exp[, 1],
                         Lynx = smoothing_exp[, 2])
  filtering <- map_dfr(1:n_obs,
                       function(k){
                         filt_parts <- particles[,, k]
                         tibble(Year = 1899 + k,
                                particle_index  = 1:n_particles,
                                Hares = filt_parts[, 1],
                                Lynx = filt_parts[, 2],
                                weight = weights[, k],
                                seed = my_seeds[i],
                                type = "filtering")
                       })
  
  list(filter = filtering, Poor_man = PM_smooth, Paris = Paris_smooth)
})

filter <- map_dfr(all_results, "filter")
Poor_man <- map_dfr(all_results, "Poor_man")
filter_estimate <- filter %>% 
  group_by(Year, seed, type) %>% 
  summarise(Hares = sum(Hares * weight),
            Lynx = sum(Lynx * weight))
Poor_man_estimate <- Poor_man%>% 
  group_by(Year, seed, type) %>% 
  summarise(Hares = sum(Hares * weight),
            Lynx = sum(Lynx * weight))
Paris_estimate <- map_dfr(all_results, "Paris")
estimates <- bind_rows(filter_estimate,
                       Poor_man_estimate,
                       Paris_estimate)
theme_set(theme_gray() +
            theme(
              panel.border = element_rect(colour = "black", 
                                          fill = rgb(0, 0, 0, 0)),
              panel.grid = element_line(linetype = 2),
              legend.background = element_blank(), 
              text = element_text(family = "Symbola", size = 24),
              strip.background = element_rect(color = "black", 
                                              fill = "lightgoldenrod1")))
ggplot(filter(estimates, type == "Paris_smooth")) +
  aes(x = Hares, y = Lynx) +
  geom_point(aes(group = interaction(seed, type)),
             alpha = 0.5, color = "darkgreen") +
  geom_path(aes(group = interaction(seed, type)),
                alpha = 0.5, color = "darkgreen", size = 0.1) +
  theme(legend.position = "none") +
  geom_point(data = hares_lynx_data[1:20,], shape = 4,
             color = "red", size = 3) + 
  # geom_path(data = hares_lynx_data) +
  geom_text_repel(data = hares_lynx_data[1:20,], aes(label = Year),
                  size = 8)
