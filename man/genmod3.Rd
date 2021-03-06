% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genmods.R
\name{genmod3}
\alias{genmod3}
\title{Third generative model}
\usage{
genmod3(n, mu, beta, variance, q)
}
\arguments{
\item{n}{Sample size.}

\item{mu}{A parameter giving the mean (times plus or minus one).}

\item{beta}{A vector of parameters with length equal to the number of covariates.}

\item{variance}{The variance of the error term.}

\item{q}{A parameter to control the proportion that is diseased.}
}
\value{
A list with the following components: \code{A}, a vector of class labels; and
  \code{X}, a matrix of covariates.
}
\description{
\code{genmod3} generates data from a simple linear generative model
  with some nonnormal covariates for simulation experiments.
}
