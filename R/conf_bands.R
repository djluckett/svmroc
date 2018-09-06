# conf_bands.R

# to delete later: I'm changing Y to A so that it matches notation from the paper (8/23/18)

# function for SVM ROC confidence bands

#' SVM ROC confidence bands
#'
#' \code{conf_band} constructs bootstrap confidence bands for the SVM ROC curve from
#'   an object of class \code{svmroc}
#'
#' @param object An object of class \code{svmroc}.
#' @param num_boot Number of bootstrap replications. Defaults to 1000.
#' @param gamma Complement of confidence level, e.g., \code{gamma = 0.1}
#'   will produce 90\% confidence bands. Defaults to 0.1.
#' @param x Values used for interpolation. Defaults to \code{seq(0.01, 0.99, 0.01)}.
#'
#' @return An object of class \code{conf_band}, a list with the following components:
#'   \code{lower}, a vector of values for the lower confidence band; \code{upper}, a
#'   vectoro of values for the upper confidence band; \code{y}, values of the ROC curve;
#'   and \code{x}, values used for interpolation.
#'
#' @export
conf_band = function(object, num_boot = 1000, gamma = 0.1, x = seq(0.01, 0.99, 0.01)) {

  # extract estimated sensitivities and specificities
  sens = object$sens
  spec = object$spec
  actual = object$new_A

  # create matrix of predicted values
  pred = list()
  for (i in 1:length(object$weights)) {

    new_pred = predict(object, object$new_X, object$weights[i])
    pred[[i]] = new_pred

  }  # end loop through weights

  # ensure matrices and get n
  sens = as.matrix(sens)
  spec = as.matrix(spec)
  #pred = as.matrix(pred)
  n = length(actual)

  # interpolate sens and spec
  y_hat = approx(x = 1 - spec, y = sens, xout = x, yleft = 0, yright = 1)$y

  # get indices where the ROC curve is not at 0 or 1
  inds = which(y_hat < 1 & y_hat > 0)

  sens_tilde = matrix(NA, nrow = length(pred), ncol = num_boot)
  spec_tilde = matrix(NA, nrow = length(pred), ncol = num_boot)

  y_tilde = matrix(NA, nrow = length(x), ncol = num_boot)

  for (b in 1:num_boot) {

    weights = rexp(n, 1)
    weights = weights / mean(weights)

    for (k in 1:length(pred)) {

      sens_tilde[k, b] = mean(weights * as.numeric(actual == levels(actual)[2]) * as.numeric(pred[[k]] == levels(actual)[2])) / mean(weights * as.numeric(actual == levels(actual)[2]))
      spec_tilde[k, b] = mean(weights * as.numeric(actual == levels(actual)[1]) * as.numeric(pred[[k]] == levels(actual)[1])) / mean(weights * as.numeric(actual == levels(actual)[1]))

    }  # end loop through alphas

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

    old_ell = ell; old_u = u

    ell = y_tilde_ordered[ , num_boot / 2 - s + 1]
    u = y_tilde_ordered[ , num_boot / 2 + s]

    # determine what proportion of bootstrap samples are captured by s steps away from the median
    cover = rep(1, num_boot)
    for (b in 1:num_boot) {

      for (k in inds) {  # amended to only look at the part of the ROC curve that is not equal to 0 or 1 (used to be 1:length(x))

        if (y_tilde[k, b] < ell[k] | y_tilde[k, b] > u[k]) {
          cover[b] = 0
          break
        }

      }  # end loop through x values

    }  # end loop through bootstrap samples

    cover_prob = mean(cover)
    if (cover_prob < 1 - gamma) {
      num_steps = s
      up = old_u
      low = old_ell
      break
    }

  }  # end loop through steps away from median

  upper = 2 * y_hat - low
  lower = 2 * y_hat - up

  upper = pmin(upper, 1)
  upper = pmax(upper, 0)
  lower = pmin(lower, 1)
  lower = pmax(lower, 0)

  to_return = list(lower = lower, upper = upper, y = y_hat, x = x)
  class(to_return) = "conf_band"

  return(to_return)

}  # end function conf_band

# function to plot confidence bands

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
#' @export
plot.conf_band = function(object, xlab = "One minus specificity", ylab = "Sensitivity",
                          include_opt = T) {

  dat = cbind(c(0, object$y, 1), c(0, object$lower, 1), c(0, object$upper, 1), c(0, object$x, 1))
  dat = as.data.frame(dat)
  names(dat) = c("y", "low", "up", "x")

  g = ggplot(dat, aes(x = x, y = y)) + geom_line(aes(x = x, y = y, linetype = "solid")) +
    geom_ribbon(aes(x = x, ymin = low, ymax = up), alpha = 0.3) +
    geom_segment(x = 0, y = 0, xend = 1, yend = 1, linetype = 2) + theme_classic() +
    theme(legend.position = "none") + labs(x = xlab, y = ylab)

  if (include_opt) {
    opt = opt_weight(object)
    g = g + geom_point(x = 1 - opt$opt_spec, y = opt$opt_sens)
  }  # end if optimal cutpoint should be added to plot

  return(g)

}  # end function plot.conf_band

