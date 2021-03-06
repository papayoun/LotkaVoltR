% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LV_model.R
\name{LV_create}
\alias{LV_create}
\title{Create}
\usage{
LV_create(a1 = c(12, 0.05, 1), a2 = c(2, 0.2, 0.1), mu0 = c(50, 20),
  gam = matrix(c(0.2, -0.1, -0.1, 0.2), nrow = 2), sigma0 = diag(1, 2))
}
\arguments{
\item{a1}{parameters for the prey dynamics}

\item{a2}{parameters for the predator dynamics}

\item{mu0}{mean}

\item{gam}{Covariance parameters for the dynamics}

\item{sigma0}{cov}
}
\value{
a LVPOD
}
\description{
function that creates a LVPOD instance
}
