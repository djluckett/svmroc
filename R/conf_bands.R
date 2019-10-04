# conf_bands.R

#' SVM ROC confidence bands
#'
#' \code{conf_band} constructs bootstrap confidence bands for the SVM ROC curve from
#'   an object of class \code{svmroc}.
#'
#' @param object An object of class \code{svmroc}.
#' @param num_boot Number of bootstrap replications. Defaults to 1000.
#' @param gamma Complement of confidence level, e.g., \code{gamma = 0.1}
#'   will produce 90\% confidence bands. Defaults to 0.1.
#' @param x Values used for interpolation. Defaults to \code{seq(0.01, 0.99, 0.01)}.
#'
#' @return An object of class \code{conf_band}, a list with the following components:
#'   \code{lower}, a vector of values for the lower confidence band; \code{upper}, a
#'   vector of values for the upper confidence band; \code{y}, values of the ROC curve;
#'   and \code{x}, values used for interpolation.
#'
#' @export
conf_band = function(object, num_boot = 1000, gamma = 0.1, x = seq(0.01, 0.99, 0.01)) {

  # extract estimated sensitivities and specificities
  sens = object$sens
  spec = object$spec
  actual = object$new_A

  # create list of predicted values
  pred = list()
  for (i in 1:length(object$weights)) {

    new_pred = predict(object, object$new_X, object$weights[i])
    pred[[i]] = new_pred

  }  # end loop through weights

  # ensure matrices and get n
  sens = as.matrix(sens)
  spec = as.matrix(spec)
  n = length(actual)

  # interpolate sens and spec
  y_hat = approx(x = 1 - spec, y = sens, xout = x, yleft = 0, yright = 1)$y

  # compute values needed for confidence band
  sens_tilde = matrix(NA, nrow = length(pred), ncol = num_boot)
  spec_tilde = matrix(NA, nrow = length(pred), ncol = num_boot)
  y_tilde = matrix(NA, nrow = length(x), ncol = num_boot)

  # loop through bootstrap replications
  for (b in 1:num_boot) {

    # create bootstrap weights
    weights = rexp(n, 1)
    weights = weights / mean(weights)

    # calculate bootstrap weighted sensitivity and specificity
    for (k in 1:length(pred)) {

      sens_tilde[k, b] = mean(weights * as.numeric(actual == levels(actual)[2]) * as.numeric(pred[[k]] == levels(actual)[2])) / mean(weights * as.numeric(actual == levels(actual)[2]))
      spec_tilde[k, b] = mean(weights * as.numeric(actual == levels(actual)[1]) * as.numeric(pred[[k]] == levels(actual)[1])) / mean(weights * as.numeric(actual == levels(actual)[1]))

    }  # end loop through alphas

    # linearly interpolate bootstrap weighted ROC curve
    y_tilde[ , b] = approx(x = 1 - spec_tilde[ , b], y = sens_tilde[ ,  b], xout = x,
                           yleft = 0, yright = 1)$y

  }  # end loop through bootstrap samples

  # sort over bootstrap samples for each x
  y_tilde_ordered = t(apply(y_tilde, 1, sort))

  # take medians over bootstrap samples for each x
  y_check = apply(y_tilde_ordered, 1, median)

  # initialize ell and u
  ell = rep(NA, length(x))
  u = rep(NA, length(x))

  # loop through steps toward median
  for (s in (num_boot / 2):1) {

    old_ell = ell
    old_u = u

    ell = y_tilde_ordered[ , num_boot / 2 - s + 1]
    u = y_tilde_ordered[ , num_boot / 2 + s]

    # determine what proportion of bootstrap samples are captured by s steps away from the median
    cover = rep(1, num_boot)
    for (b in 1:num_boot) {

      for (k in 1:length(x)) {

        if (y_tilde[k, b] < ell[k] | y_tilde[k, b] > u[k]) {
          cover[b] = 0
          break
        }

      }  # end loop through x values

    }  # end loop through bootstrap samples

    # check coverage proportion
    cover_prob = mean(cover)
    if (cover_prob < 1 - gamma) {
      break
    }

  }  # end loop through steps away from median

  # set limits equal to the last ones to obtain coverage 1 - gamma across bootstrap replications
  up = old_u
  low = old_ell

  upper = up
  lower = low

  # truncate upper and lower values at 0 and 1
  upper = pmin(upper, 1)
  upper = pmax(upper, 0)
  lower = pmin(lower, 1)
  lower = pmax(lower, 0)
  if (sum(lower >= 0.95) > 0) {
    temp_inds = which(lower >= 0.95)
    lower[temp_inds] = approx(x = c(x[min(temp_inds)], 1), y = c(lower[min(temp_inds)], 1), xout = x[temp_inds], yleft = 0, yright = 1)$y
  }

  # create object to return
  to_return = list(lower = lower, upper = upper, y = y_hat, x = x)
  class(to_return) = "conf_band"

  return(to_return)

}  # end function conf_band

#' Plot SVM ROC curve confidence bands
#'
#' \code{plot.conf_band} produces a plot of the ROC curve with confidence bands
#'   from an object of class conf_band.
#'
#' @param object An object of class \code{conf_band}.
#' @param xlab Label for X axis. Defaults to "One minus specificity".
#' @param ylab Label for Y axis. Defaults to "Sensitivity".
#' @param include_opt Logical. If \code{TRUE}, the optimal point on the ROC curve
#'   (the closest to (0, 1) in Euclidean distance) is marked on the plot.
#'   Defaults to \code{TRUE}.
#'
#' @return A plot as a object of class \code{ggplot}.
#'
#' @export
plot.conf_band = function(object, xlab = "One minus specificity", ylab = "Sensitivity",
                          include_opt = T) {

  # create data frame for plotting
  dat = cbind(c(0, object$y, 1), c(0, object$lower, 1), c(0, object$upper, 1), c(0, object$x, 1))
  dat = as.data.frame(dat)
  names(dat) = c("y", "low", "up", "x")

  # create ggplot object
  g = ggplot(dat, aes(x = x, y = y)) + geom_line(aes(x = x, y = y, linetype = "solid")) +
    geom_ribbon(aes(x = x, ymin = low, ymax = up), alpha = 0.3) +
    geom_segment(x = 0, y = 0, xend = 1, yend = 1, linetype = 2) + theme_classic() +
    theme(legend.position = "none") + labs(x = xlab, y = ylab)

  # add optimal point
  if (include_opt) {
    opt = opt_weight(object)
    g = g + geom_point(x = 1 - opt$opt_spec, y = opt$opt_sens)
  }  # end if optimal cutpoint should be added to plot

  return(g)

}  # end function plot.conf_band


#' Calculate area between the curves
#'
#' \code{abc} is used to calculate the area between two confidence band curves.
#'
#' @param object An object to calculate area between the curves.
#'
#' @return The numeric area between the curves.
#'
#' @export
abc = function(object) {

  UseMethod("abc", object)

}  # end function abc

# function to compute area between the curve for a conf_band object
# calls the generic directly, so there is no need for documentation
#' @export
abc.conf_band = function(object) {

  # calculate area under the upper curve
  y = object$upper
  x = object$x
  x = c(0, x, 1)
  y = c(0, y, 1)
  idx = 2:length(x)
  auc_upper = abs(as.double((x[idx] - x[idx - 1]) %*% (y[idx] + y[idx - 1])) / 2)

  # calculate area under the lower curve
  y = object$lower
  x = object$x
  x = c(0, x, 1)
  y = c(0, y, 1)
  idx = 2:length(x)
  auc_lower = abs(as.double((x[idx] - x[idx - 1]) %*% (y[idx] + y[idx - 1])) / 2)

  # return area between curves
  return(auc_upper - auc_lower)

}  # end function abc.conf_band


#' Determine coverage
#'
#' \code{get_coverage} is used to determine if confidence band covers true ROC curve.
#'
#' @param conf_band An object of class \code{conf_band}.
#' @param svmroc An object of class \code{svmroc}.
#' @param new_X New X matrix to determine true ROC curve.
#' @param new_A New class assignments to determine true ROC curve.
#'
#' @return A list with the following components: \code{coverage}, numeric, 1 if
#'   the confidence band covers the true ROC curve and 0 if it does not; \code{true_roc},
#'   a vector of points on the true ROC curve; \code{x}, the x values for the true ROC
#'   curve (taken from the object \code{conf_band}).
#'
#' @export
get_coverage = function(conf_band, svmroc, new_X, new_A) {

  # check that levels(new_A) match levels used in the original fit
  if (!(identical(levels(as.factor(new_A)), levels(as.factor(svmroc$new_A))))) {
    stop("'as.factor(new_A)' must have the same levels as the class labels in the original fit.")
  }

  # check that new_X has the same number of columns as used in the original fit
  if (ncol(new_X) != ncol(svmroc$new_X)) {
    stop("'new_X' must have the same number of columns as the covariate matrix used to obtain 'svmroc'")
  }

  # estimate true ROC curve
  roc = fit_roc(object = svmroc, new_X = new_X, new_A = new_A)

  # linear interpolation to get true ROC curve
  true_roc = approx(x = 1 - roc$spec, y = roc$sens, xout = conf_band$x, yleft = 0, yright = 1)$y

  # check coverage only at indices larger than some delta
  inds = which(conf_band$x > 0.05 & conf_band$x < 0.95)

  # get indicator of coverage
  coverage = 1 - as.numeric(sum(true_roc[inds] > conf_band$upper[inds]) + sum(true_roc[inds] < conf_band$lower[inds]) > 0)

  return(coverage)

}  # end function get_coverage


