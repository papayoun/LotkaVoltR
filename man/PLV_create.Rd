% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ProposalLVModel.R
\name{PLV_create}
\alias{PLV_create}
\title{Create}
\usage{
PLV_create(L = list(a1 = c(0.5, 0.1, 0.2), a2 = c(0.1, 0.4, 0.1), gam =
  matrix(c(1, -0.3, -0.3, 1), nrow = 2), mu0 = c(50, 20), sigma0 = diag(1,
  2), cov = diag(1, 2), qs = c(0.2, 0.3), RWC = diag(10, 2), wO = 1, wD =
  1))
}
\arguments{
\item{a1}{parameters for the prey dynamics}

\item{a2}{parameters for the predator dynamics}

\item{gam}{Covariance parameters for the dynamics}

\item{cov}{Covariance parameters for the observation model}

\item{qs}{Capturability parameter}

\item{RWC}{random walk covariance}
}
\value{
a Proposal LVPOD
}
\description{
function that creates a LVPOD instance
}
