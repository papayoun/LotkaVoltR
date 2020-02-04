## ----hares_lynx, echo = F------------------------------------------------
hares_lynx_data <- tribble(~Year,	~Hares,	~Lynx,
                           1900,	30,	4,
                           1901,	47.2,	6.1,
                           1902,	70.2,	9.8,
                           1903,	77.4,	35.2,
                           1904,	36.3,	59.4,
                           1905,	20.6,	41.7,
                           1906,	18.1,	19,
                           1907,	21.4,	13,
                           1908,	22,	8.3,
                           1909,	25.4,	9.1,
                           1910,	27.1,	7.4,
                           1911,	40.3,	8,
                           1912,	57,	12.3,
                           1913,	76.6,	19.5,
                           1914,	52.3,	45.7,
                           1915,	19.5,	51.1,
                           1916,	11.2,	29.7,
                           1917,	7.6,	15.8,
                           1918,	14.6,	9.7,
                           1919,	16.2,	10.1,
                           1920,	24.7,	8.6)


## ----showing_data, echo = F----------------------------------------------
knitr::kable(hares_lynx_data)


## ----pressure, echo=FALSE------------------------------------------------
a1 <- c(0.695617134, 0.001147425, 0.030333635) # Prey parameters
a2 <- c(0.630905980, 0.021756719, 0.001352425) # Predator parameters
Gamma <- matrix(c(0.22, 0, 0, 0.2), nrow = 2) # Diffusion parameters
mu0 <- c(30, 4) # Mean of the inital distribution
Sigma0 <- diag(1, 2) # Variance of the initial distribution

# Observation process parameters

Sigma_obs <- matrix(c(0.01, 0.00, 0.00, 0.01), ncol = 2) # Observation noise
q_values <- c(1, 1) # Known q values