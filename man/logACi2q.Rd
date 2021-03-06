% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions2.R, R/Functions5.R
\name{logACi2q}
\alias{logACi2q}
\title{Log of the Ascertainment correction piece for bivariate sampling}
\usage{
logACi2q(yi, xi, zi, wi, beta, sigma.vc, rho.vc, sigma.e, cutpoints,
  SampProb)

logACi2q(yi, xi, zi, wi, beta, sigma.vc, rho.vc, sigma.e, cutpoints,
  SampProb)
}
\arguments{
\item{yi}{n_i-response vector}

\item{xi}{n_i by p design matrix for fixed effects}

\item{zi}{n_i by 2 design matric for random effects (intercept and slope)}

\item{wi}{the pre-multiplier of yi to generate the sampling variable q_i}

\item{beta}{mean model parameter p-vector}

\item{sigma.vc}{vector of variance components on standard deviation scale}

\item{rho.vc}{vector of correlations among the random effects.  The length should be q choose 2}

\item{sigma.e}{std dev of the measurement error distribution}

\item{cutpoints}{cutpoints defining the sampling regions. (a vector of length 4 c(xlow, xhigh, ylow, yhigh))}

\item{SampProb}{Sampling probabilities from within each region (vector of length 2 c(central region, outlying region)).}

\item{sigma0}{std dev of the random intercept distribution}

\item{sigma1}{std dev of the random slope distribution}

\item{rho}{correlation between the random intercept and slope}

\item{sigmae}{std dev of the measurement error distribution}

\item{yi}{n_i-response vector}

\item{xi}{n_i by p design matrix for fixed effects}

\item{zi}{n_i by q design matric for random effects (intercept and slope)}

\item{wi}{the pre-multiplier of yi to generate the sampling variable q_i}

\item{beta}{mean model parameter p-vector}

\item{cutpoints}{cutpoints defining the sampling regions. (a vector of length 4 c(xlow, xhigh, ylow, yhigh))}

\item{SampProb}{Sampling probabilities from within each region (vector of length 2 c(central region, outlying region)).}
}
\value{
log transformed ascertainment correction

log transformed ascertainment correction
}
\description{
Calculate the log transformed ascertainment correction under a bivariate Q_i.  Also return vi

Calculate the log transformed ascertainment correction under a bivariate Q_i.  Also return vi
}
