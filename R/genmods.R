# genmods.R

#' First generative model
#'
#' \code{genmod1} generates data from a simple linear generative model
#'   for simulation experiments.
#'
#' @param n Sample size.
#' @param mu A parameter giving the mean (times plus or minus one).
#' @param beta A vector of parameters with length equal to the number of covariates.
#' @param variance The variance of the error term.
#' @param q A parameter to control the proportion that is diseased.
#'
#' @return A list with the following components: \code{A}, a vector of class labels; and
#'   \code{X}, a matrix of covariates.
#'
#' @export
genmod1 = function(n, mu, beta, variance, q) {

  # create an empty matrix to fill
  dat = NULL

  # loop through sample size
  for (i in 1:n) {

    # generate a random uniform number and set mean
    z = 2 * (runif(1) <= q) - 1

    # generate covariates
    X = rnorm(length(beta), mu * z, sqrt(variance))

    # generate probability of disease
    p = boot::inv.logit(sum(X * beta))

    # generate outcome
    A = 2 * rbinom(1, 1, p) - 1

    # add data
    dat = rbind(dat, c(A, X))

  }  # end loop through observations

  # return list with X and A
  return(list(A = dat[ , 1], X = dat[ , 2:ncol(dat)]))

}  # end function genmod1

#' Second generative model
#'
#' \code{genmod2} generates data from a simple nonlinear generative model for
#'   simulation experiments.
#'
#' @param n Sample size.
#' @param mu A parameter giving the mean (times plus or minus one).
#' @param beta A vector of parameters with length equal to the number of covariates.
#' @param variance The variance of the error term.
#' @param q A parameter to control the proportion that is diseased.
#'
#' @return A list with the following components: \code{A}, a vector of class labels; and
#'   \code{X}, a matrix of covariates.
#'
#' @export
genmod2 = function(n, mu, beta, variance, q) {

  # create an empty matrix to fill
  dat = NULL

  # loop through sample size
  for (i in 1:n) {

    # generate a random uniform number and set mean
    z = 2 * (runif(1) <= q) - 1

    # generate covariates
    X = rnorm(length(beta), mu * z, sqrt(variance))

    # generate probability of disease
    p = boot::inv.logit(sum(X * beta) + X[1]^2 + X[2]^2 + 4 * X[1] * X[2])

    # generate outcome
    A = 2 * rbinom(1, 1, p) - 1

    # add data
    dat = rbind(dat, c(A, X))

  }  # end loop through observations

  # return list with X and A
  return(list(A = dat[ , 1], X = dat[ , 2:ncol(dat)]))

}  # end function genmod2

#' Third generative model
#'
#' \code{genmod3} generates data from a simple linear generative model
#'   with some nonnormal covariates for simulation experiments.
#'
#' @param n Sample size.
#' @param mu A parameter giving the mean (times plus or minus one).
#' @param beta A vector of parameters with length equal to the number of covariates.
#' @param variance The variance of the error term.
#' @param q A parameter to control the proportion that is diseased.
#'
#' @return A list with the following components: \code{A}, a vector of class labels; and
#'   \code{X}, a matrix of covariates.
#'
#' @export
genmod3 = function(n, mu, beta, variance, q) {

  # create an empty matrix to fill
  dat = NULL

  # loop through sample size
  for (i in 1:n) {

    # generate a random uniform number and set mean
    z = 2 * (runif(1) <= q) - 1

    # generate normal covariates
    X_norm = rnorm(length(beta) - 2, mu * z, sqrt(variance))

    # generate nonnormal covariates
    X_nonnorm = rep(NA, 2)
    X_nonnorm[1] = rbinom(1, 1, boot::inv.logit(X_norm[1]))
    X_nonnorm[2] = rpois(1, exp(X_norm[1] / 4))

    # assemble covariates
    X = c(X_nonnorm, X_norm)

    # generate probability of disease
    p = boot::inv.logit(sum(X * beta))

    # generate outcome
    A = 2 * rbinom(1, 1, p) - 1

    # add data
    dat = rbind(dat, c(A, X))

  }  # end loop through observations

  # return list with X and A
  return(list(A = dat[ , 1], X = dat[ , 2:ncol(dat)]))

}  # end function genmod3

