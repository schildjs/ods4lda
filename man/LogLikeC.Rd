% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions2.R
\name{LogLikeC}
\alias{LogLikeC}
\title{Calculate the conditional likelihood for the univariate and bivariate sampling cases across all subjects (Keep.liC=FALSE) or the subject specific contributions to the conditional likelihood along with the log-transformed ascertainment correction for multiple imputation (Keep.liC=TRUE).}
\usage{
LogLikeC(y, x, z, w.function, id, beta, sigma0, sigma1, rho, sigmae,
  cutpoints, SampProb, SampProbi, Keep.liC = FALSE)
}
\arguments{
\item{y}{response vector}

\item{x}{sum(n_i) by p design matrix for fixed effects}

\item{z}{sum(n_i) by 2 design matric for random effects (intercept and slope)}

\item{w.function}{options include "mean" "intercept" "slope" and "bivar"}

\item{id}{sum(n_i) vector of subject ids}

\item{beta}{mean model parameter p-vector}

\item{sigma0}{std dev of the random intercept distribution}

\item{sigma1}{std dev of the random slope distribution}

\item{rho}{correlation between the random intercept and slope}

\item{sigmae}{std dev of the measurement error distribution}

\item{cutpoints}{cutpoints defining the sampling regions. (a vector of length 4 c(xlow, xhigh, ylow, yhigh))}

\item{SampProb}{Sampling probabilities from within each region (vector of length 2 c(central region, outlying region)).}

\item{SampProbi}{Subject specific sampling probabilities.  A vector of length sum(n_i).  Not used unless using weighted Likelihood}

\item{Keep.liC}{If FALSE, the function returns the conditional log likelihood across all subjects.  If TRUE, subject specific contributions and exponentiated subject specific ascertainment corrections are returned in a list.}
}
\value{
If Keep.liC=FALSE, conditional log likelihood.  If Keep.liC=TRUE, a two-element list that contains subject specific likelihood contributions and exponentiated ascertainment corrections.
}
\description{
Calculate the conditional likelihood for the univariate and bivariate sampling cases across all subjects (Keep.liC=FALSE) or the subject specific contributions to the conditional likelihood along with the log-transformed ascertainment correction for multiple imputation (Keep.liC=TRUE).
}
