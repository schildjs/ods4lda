% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions5.R
\name{logACi1q.score2}
\alias{logACi1q.score2}
\title{Gradient of the log of the ascertainment correction piece for sampling based on univariate Q_i}
\usage{
logACi1q.score2(subjectData, beta, sigma.vc, rho.vc, sigma.e)
}
\arguments{
\item{subjectData}{a list containing: yi, xi, zi, Weights.i, w.function.i, SampProb.i, cutpoints.i}

\item{beta}{mean model parameter p-vector}

\item{sigma.vc}{vector of variance components on standard deviation scale}

\item{rho.vc}{vector of correlations among the random effects.  The length should be q choose 2}

\item{sigma.e}{std dev of the measurement error distribution}
}
\value{
gradient of the log transformed ascertainment correction under univariate $Q_i$
}
\description{
Calculate the gradient of the log transformed ascertainment correction for sampling based on univariate Q_i
}
