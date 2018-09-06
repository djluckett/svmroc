# genmods.R

# to delete later: I'm changing Y to A so that it matches notation from the paper (8/23/18)

#' First generative model
#'
#' \code{genmod1} generates data from a simple linear generative model
#'   for simulation experiments.
#'
#' @param n The sample size.
#' @param beta A vector of parameters with length equal to the number of covariates.
#' @param variance The variance of the error term.
#'
#' @return A list with the following components: \code{A}, a vector of class labels; and
#'   \code{X}, a matrix of covariates.
#'
#' @export
genmod1 = function(n, beta, variance = 1) {

  dat = NULL

  for (i in 1:n) {

    X = rnorm(length(beta), 0, 1)
    eps = rnorm(1, 0, sqrt(variance))
    A = sign(sum(X * beta) + eps)

    dat = rbind(dat, c(A, X))

  }  # end loop through observations

  return(list(A = dat[ , 1], X = dat[ , 2:ncol(dat)]))

}  # end function genmod1

#' Second generative model
#'
#' \code{genmod2} generates data from a simple nonlinear generative model for
#'   simulation experiments.
#'
#' @param n The sample size.
#' @param beta A vector of parameters with length equal to the number of covariates.
#' @param variance The variance of the error term.
#'
#' @return A list with the following components: \code{A}, a vector of class labels; and
#'   \code{X}, a matrix of covariates.
#'
#' @export
genmod2 = function(n, beta, variance = 1) {

  dat = NULL

  for (i in 1:n) {

    X = rnorm(length(beta), 0, 1)
    eps = rnorm(1, 0, sqrt(variance))
    A = sign(sum(X * beta) + X[1]^2 - X[2]^2 + X[1]^3 + exp(X[1]) - 1.5 + eps)

    dat = rbind(dat, c(A, X))

  }  # end loop through observations

  return(list(A = dat[ , 1], X = dat[ , 2:ncol(dat)]))

}  # end function genmod2

