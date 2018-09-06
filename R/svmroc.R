# svmroc.R

# to delete later: I'm changing Y to A so that it matches notation from the paper (8/23/18)

#' SVM ROC curves
#'
#' \code{svmroc} fits a weighted support vector machine to estimate a receiver operating
#'   characteristic curve
#'
#' @param X A matrix of covariates used for fitting the SVM.
#' @param A A factor with two levels or an object that can be coerced to a factor with
#'   two levels; used as the response when fitting the SVM.
#' @param new_X A matrix of covariates used as the testing set for estimating the ROC curve.
#' @param new_A A factor with two levels or an object that can be coerced to a factor with
#'   two levels; used as the response for the testing set for estimating the ROC curve.
#'   Must have the same levels as \code{Y}.
#' @param kernel The kernel used for fitting the SVM, either "linear" or "Gaussian".
#'   Defaults to "linear".
#' @param lambdas A numeric vector of penalty parameters to select from. Defaults to
#'   \code{c(0.1, 0.25, 0.5, 1)}.
#' @param sigmas A numeric vector of bandwidth parameters to select from. Defaults to
#'   \code{c(0.1, 0.25, 0.5, 1)}.
#' @param num_folds Number of folds for cross-validation to select tuning parameters.
#'   Defaults to 10.
#' @param weights A vector of weights to use for the weighted SVM. Defaults to
#'   \code{seq(0.01, 0.99, 0.01)}.
#' @param seed Random number seed to set. If \code{NULL}, no seed is set. Defaults to \code{NULL}.
#'
#' @return An object of class \code{svmroc}, a list with the following components:
#'   \code{sens}, estimated sensitivities across weights; \code{spec}, estimated
#'   specificities across weights; \code{models}, a list of models fit for each weight;
#'   \code{weights}, the vector of weights used in the fit;
#'   \code{new_X}, the matrix of covariates in the testing set; and \code{new_Y}, the
#'   responses in the testing set.
#'
#' @export
svmroc = function(X, A, new_X, new_A, kernel = "linear", lambdas = c(0.1, 0.25, 0.5, 1),
                  sigmas = c(0.1, 0.25, 0.5, 1), num_folds = 10, weights = seq(0.01, 0.99, 0.01),
                  seed = NULL) {

  # a note on using this function: you can pass as A and new_A a factor or anything that can
  # coerced to a factor. If both are factors, they must have the same levels. Otherwise, they
  # must be able to be coerced into factors such that both will have the same levels after
  # applying as.factor(). Throughout, levels(A)[2] is considered "disease" and levels(A)[1]
  # is considered nondisease (when defining sensitivity and specificity). Thus, if A and new_A
  # have values 1 and -1, the factor levels will work ask expected.

  # set seed if one is provided
  if (!is.null(seed)) set.seed(seed)

  # prepare data
  A = as.factor(A)
  new_A = as.factor(new_A)

  # cross-validation to select tuning parameters
  best_lambda = NA
  best_sigma = NA
  best_error = Inf
  class_weights = c(0.5, 0.5)
  for (i in 1:length(lambdas)) {

    if (kernel == "Gaussian") {

      for (j in 1:length(sigmas)) {

        kern = rbfdot(sigma = sigmas[j])
        fit = ksvm(x = X, y = A, class.weights = class_weights, cross = num_folds,
                   kernel = kern, C = lambdas[i])

        if (fit@cross < best_error) {
          best_lambda = lambdas[i]
          best_sigma = sigmas[j]
          best_error = fit@cross
        }  # end if CV error is best so far

      }  # end loop through sigmas

    }  # end Gaussian kernel

    if (kernel == "linear") {

      kern = vanilladot()
      fit = ksvm(x = X, y = A, class.weights = class_weights, cross = num_folds,
                 kernel = kern, C = lambdas[i])

      if (fit@cross < best_error) {
        best_lambda = lambdas[i]
        best_error = fit@cross
      }  # end if CV error is best so far

    }  # end linear kernel

  }  # end loop through penalty parameters

  # vectors to save all sensitivities and specificities
  se = rep(NA, length(weights))
  sp = rep(NA, length(weights))
  models = list()

  # loop through weights
  for (j in 1:length(weights)) {

    if (kernel == "linear") kern = vanilladot()
    if (kernel == "Gaussian") kern = rbfdot(sigma = best_sigma)

    class_weights = c(1 - weights[j], weights[j])
    names(class_weights) = levels(A)
    fit = ksvm(x = X, y = A, class.weights = class_weights, kernel = kern, C = best_lambda)
    pred_vals = predict(object = fit, newdata = new_X, type = "response")

    spec = sum(new_A == levels(new_A)[1] & pred_vals == levels(new_A)[1]) / sum(new_A == levels(new_A)[1])

    sens = sum(new_A == levels(new_A)[2] & pred_vals == levels(new_A)[2]) / sum(new_A == levels(new_A)[2])

    se[j] = sens
    sp[j] = spec

    models[[j]] = fit

  }  # end loop through alphas

  to_return = list(sens = se, spec = sp, models = models, weights = weights,
                   new_X = new_X, new_A = new_A)
  class(to_return) = "svmroc"

  return(to_return)

}  # end function svmroc

# function to plot ROC curve for the SVM

#' Plot SVM ROC curve
#'
#' \code{plot.svmroc} returns a plot of the ROC curve associated
#' with an object of class \code{svmroc}.
#'
#' @param object An object of class \code{svmroc}.
#' @param xlab Label for X axis. Defaults to "One minus specificity".
#' @param ylab Label for Y axis. Defaults to "Sensitivity".
#' @param include_opt Logical. If \code{TRUE}, the optimal point on the ROC curve
#'   (the closest to (0, 1) in Euclidean distance) is marked on the plot.
#'   Defaults to \code{TRUE}.
#'
#' @export
plot.svmroc = function(object, xlab = "One minus specificity", ylab = "Sensitivity",
                       include_opt = T) {

  to_plot = data.frame(fpf = c(0, 1 - object$spec, 1), tpf = c(0, object$sens, 1))

  g = ggplot(to_plot, aes(x = fpf, y = tpf)) + geom_line() + xlim(0, 1) + ylim(0, 1) +
             theme_classic() + labs(x = xlab, y = ylab) +
             geom_segment(x = 0, y = 0, xend = 1, yend = 1, linetype = 2) +
             theme(legend.position = "none")

  if (include_opt) {
    opt = opt_weight(object)
    g = g + geom_point(x = 1 - opt$opt_spec, y = opt$opt_sens)
  }  # end if optimal cutpoint should be added to plot

  return(g)

}  # end function plot.svmroc

# generic function for area under the curve

#' Compute area under the curve
#'
#' \code{auc} is used to compute the area under the ROC curve for
#' an object of class \code{svmroc}.
#'
#' @param object An object of class \code{svmroc}
#'
#' @return The area under the curve.
#'
#' @export
auc = function(object) {

  UseMethod("auc", object)

}  # end function auc

# function to compute area under the curve
# calls the generic directly, so there is no need for documentation
#' @export
auc.svmroc = function(object) {

  tpf = object$sens
  fpf = 1 - object$spec
  if (fpf[length(fpf)] < fpf[1]) fpf = rev(fpf)
  if (tpf[length(tpf)] < tpf[1]) tpf = rev(tpf)
  fpf = c(0, fpf, 1); tpf = c(0, tpf, 1)
  idx = 2:length(fpf)

  return(abs(as.double( (fpf[idx] - fpf[idx - 1]) %*% (tpf[idx] + tpf[idx - 1])) / 2))

}  # end function auc.svmroc

#' Predict method for \code{svmroc}
#'
#' \code{predict.svmroc} produces predicted class assignments based on the fitted SVM.
#'
#' @param object An object of class \code{svmroc}
#' @param newdata New data on which to predict.
#' @param weight The weight to used for prediction. If \code{NULL}, the optimal weight
#'   is calculated and used. Defaults to \code{NULL}.
#'
#' @return The vector of predicted class assignments.
#'
#' @export
predict.svmroc = function(object, newdata, weight = NULL) {

  if (!is.null(weight)) {
    if(!(weight %in% object$weights)) {
      stop("'weight' must be in the vector of weights used when creating 'object'")
    }
  } else {
    opt = opt_weight(object)
    weight = opt$weight
  }

  index = which(weight == object$weights)
  model = object$models[[index]]
  to_return = predict(model, newdata, type = "response")

  return(to_return)

}  # end function predict.svmroc

# generic function for optimal cutpoint

#' Calculate optimal weight
#'
#' \code{opt_weight} computes the optimal point on the ROC curve, calculated as the point
#'   closest to (0, 1) in Euclidean distance.
#'
#' @param object An object of class \code{svmroc}.
#'
#' @return A list with the following components: \code{weight}, the optimal weight;
#'   \code{opt_sens}, the sensitivity at the optimal weight; \code{opt_spec}, the
#'   specificity at the optimal weight.
#'
#' @export
opt_weight = function(object) {

  UseMethod("opt_weight", object)

}  # end function opt_weight

#' Calculate optimal weight
#'
#' \code{opt_weight.svmroc} computes the optimal point on the ROC curve, calculated as the point
#'   closest to (0, 1) in Euclidean distance from an object of class \code{svmroc}.
#'
#' @param object An object of class \code{svmroc}.
#'
#' @return A list with the following components: \code{weight}, the optimal weight;
#'   \code{opt_sens}, the sensitivity at the optimal weight; \code{opt_spec}, the
#'   specificity at the optimal weight.
#'
#' @export
opt_weight.svmroc = function(object) {

  value = Inf
  for (i in 1:length(object$spec)) {
    new = sqrt((1 - object$sens[i])^2 + (0 - (1 - object$spec[i]))^2)
    if (new < value) {
      value = new; weight = object$weights[i]; index = i
    }
  }

  return(list(weight = weight, opt_sens = object$sens[index], opt_spec = object$spec[index]))

}  # end function opt_weight.svmroc

# function to select optimal ROC point after bootstrap


#' Calculate optimal weight
#'
#' \code{opt_weight.conf_band} computes the optimal point on the ROC curve,
#'   calculated as the point closest to (0, 1) in Euclidean distance for an
#'   object of class \code{conf_band}.
#'
#' @param object An object of class \code{conf_band}.
#'
#' @return A list with the following components: \code{opt_sens}, the sensitivity
#'   at the optimal weight; \code{opt_spec}, the specificity at the optimal weight.
#'
opt_weight.conf_band = function(object) {

  value = Inf
  for (i in 1:length(object$x)) {
    new = sqrt((1 - object$y[i])^2 + (0 - object$x[i])^2)
    if (new < value) {
      value = new; index = i
    }
  }

  return(list(opt_sens = object$y[index], opt_spec = 1 - object$x[index]))

}  # end function opt_weight.conf_band
