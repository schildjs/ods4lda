% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions2.R
\name{li.lme.score}
\alias{li.lme.score}
\title{Subject specific contribution to the lme model score (also returns marginal Vi=Cov(Y|X))}
\usage{
li.lme.score(subjectData, beta, sigma0, sigma1, rho, sigmae)
}
\arguments{
\item{subjectData}{a list that contains yi, xi, zi, SampProbi.i.  Note that SampProbi.i is used for IPW only.}

\item{beta}{mean model parameter p-vector}

\item{sigma0}{std dev of the random intercept distribution}

\item{sigma1}{std dev of the random slope distribution}

\item{rho}{correlation between the random intercept and slope}

\item{sigmae}{std dev of the measurement error distribution}
}
\value{
Subject specific contribution to the log-likelihood score (also returns marginal Vi=Cov(Y|X))
}
\description{
Subject specific contribution to the lme model score (also returns marginal Vi=Cov(Y|X))
}
