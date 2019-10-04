# svmroc.R

#' SVM ROC curves
#'
#' \code{svmroc} fits a weighted support vector machine to estimate a receiver operating
#'   characteristic curve.
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
#'   \code{c(0.5, 1, 2, 4, 8)}.
#' @param sigmas A numeric vector of bandwidth parameters to select from. Defaults to
#'   \code{c(0.05, 0.1, 0.5, 1, 2)}.
#' @param num_folds Number of folds for cross-validation to select tuning parameters.
#'   Defaults to 5.
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
svmroc = function(X, A, new_X, new_A, kernel = "linear", lambdas = c(0.5, 2, 4, 8, 12, 16),
                  sigmas = c(0.05, 0.1, 0.5, 1, 2), num_folds = 5, weights = seq(0.01, 0.99, 0.01),
                  seed = NULL) {

  # a note on using this function: you can pass as A and new_A a factor or anything that can be
  # coerced to a factor. If both are factors, they must have the same levels. Otherwise, they
  # must be able to be coerced into factors such that both will have the same levels after
  # applying as.factor(). Throughout, levels(A)[2] is considered "disease" and levels(A)[1]
  # is considered nondisease (when defining sensitivity and specificity). Thus, if A and new_A
  # have values 1 and -1, the factor levels will work ask expected.

  # set seed if one is provided
  if (!is.null(seed)) set.seed(seed)

  # coerce response to factor
  A = as.factor(A)
  new_A = as.factor(new_A)

  # check that A and new_A have the same levels
  if(!(identical(levels(A), levels(new_A)))) {
    stop("'as.factor(A)' and 'as.factor(new_A)' must have the same levels.")
  }

  # check that X and new_X have the same number of columns
  if (ncol(X) != ncol(new_X)) {
    stop("'X' and 'new_X' must have the same number of columns.")
  }

  # check that a valid kernel has been supplied
  if (!(kernel %in% c("linear", "Gaussian"))) {
    stop("'kernel' must be either 'linear' or 'Gaussian'")
  }

  # cross-validation to select tuning parameters
  if (length(lambdas) == 1) {
    best_lambda = lambdas
  } else {
    best_lambda = NA
    best_sigma = NA
    best_error = Inf
    class_weights = c(0.5, 0.5)

    # loop through provided lambdas
    for (i in 1:length(lambdas)) {

      # Gaussian kernel
      if (kernel == "Gaussian") {

        # loop through provided sigmas
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

      # linear kernel
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

  }  # end cross-validation

  # vectors to save all sensitivities and specificities
  se = rep(NA, length(weights))
  sp = rep(NA, length(weights))
  models = list()

  # loop through weights
  for (j in 1:length(weights)) {

    # set kernel function
    if (kernel == "linear") kern = vanilladot()
    if (kernel == "Gaussian") kern = rbfdot(sigma = best_sigma)

    # set class weights, fit SVM, and compute predicted class labels on testing set
    class_weights = c(1 - weights[j], weights[j])
    names(class_weights) = levels(A)
    fit = ksvm(x = X, y = A, class.weights = class_weights, kernel = kern, C = best_lambda)
    pred_vals = predict(object = fit, newdata = new_X, type = "response")

    # calculate sensitivity and specificity
    spec = sum(new_A == levels(new_A)[1] & pred_vals == levels(new_A)[1]) / sum(new_A == levels(new_A)[1])
    sens = sum(new_A == levels(new_A)[2] & pred_vals == levels(new_A)[2]) / sum(new_A == levels(new_A)[2])

    # save sensitivity and specificity
    se[j] = sens
    sp[j] = spec

    # save model fit
    models[[j]] = fit

  }  # end loop through alphas

  # create object to return
  to_return = list(sens = se, spec = sp, models = models, weights = weights,
                   new_X = new_X, new_A = new_A, best_lambda = best_lambda,
                   best_sigma = best_sigma)
  class(to_return) = "svmroc"

  return(to_return)

}  # end function svmroc

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
#' @return A plot as a object of class \code{ggplot}.
#'
#' @export
plot.svmroc = function(object, xlab = "One minus specificity", ylab = "Sensitivity",
                       include_opt = T) {

  # create data frame for plotting
  to_plot = data.frame(fpf = c(0, 1 - object$spec, 1), tpf = c(0, object$sens, 1))

  # create ggplot object
  g = ggplot(to_plot, aes(x = fpf, y = tpf)) + geom_line() + xlim(0, 1) + ylim(0, 1) +
             theme_classic() + labs(x = xlab, y = ylab) +
             geom_segment(x = 0, y = 0, xend = 1, yend = 1, linetype = 2) +
             theme(legend.position = "none")

  # add optimal point
  if (include_opt) {
    opt = opt_weight(object)
    g = g + geom_point(x = 1 - opt$opt_spec, y = opt$opt_sens)
  }  # end if optimal cutpoint should be added to plot

  return(g)

}  # end function plot.svmroc

#' Compute area under the curve
#'
#' \code{auc} is used to compute the area under the ROC curve for
#' an object of class \code{svmroc}.
#'
#' @param object An object of class \code{svmroc} or \code{roc}.
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

  # calculate AUC based on sensitivity and specificity vectors returned by svmroc
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
#' @param newdata New data on which to predict, a matrix with the same columns as
#'   \code{X} used for the fit producing \code{object}.
#' @param weight The weight to used for prediction. If \code{NULL}, the optimal weight
#'   is calculated and used. Defaults to \code{NULL}.
#'
#' @return The vector of predicted class assignments.
#'
#' @export
predict.svmroc = function(object, newdata, weight = NULL) {

  # check that weight is in the vector of weights used to fit and compute weight if needed
  if (!is.null(weight)) {
    if(!(weight %in% object$weights)) {
      stop("'weight' must be in the vector of weights used when creating 'object'")
    }
  } else {
    opt = opt_weight(object)
    weight = opt$weight
  }

  # check that newdata has the same number of columns as used in the original fit
  if (ncol(newdata) != ncol(object$new_X)) {
    stop("'newdata' must have the same number of columns as the covariate matrix used to obtain 'object'")
  }

  # selected the corresponding model and get predicted values
  index = which(weight == object$weights)
  model = object$models[[index]]
  to_return = predict(model, newdata, type = "response")

  return(to_return)

}  # end function predict.svmroc

#' Calculate optimal weight
#'
#' \code{opt_weight} computes the optimal point on the ROC curve, calculated as the point
#'   closest to (0, 1) in Euclidean distance.
#'
#' @param object An object of class \code{svmroc} or \code{conf_band}.
#'
#' @return A list with the following components: \code{weight}, the optimal weight;
#'   \code{opt_sens}, the sensitivity at the optimal weight; \code{opt_spec}, the
#'   specificity at the optimal weight.
#'
#' @export
opt_weight = function(object) {

  UseMethod("opt_weight", object)

}  # end function opt_weight

# function to compute the optimal weight
# calls the generic directly, so there is no need for documentation
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

# function to compute the optimal weight
# calls the generic directly, so there is no need for documentation
#' @export
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

#' Estimate ROC curve for a fitted classifier using a new testing set
#'
#' \code{fit_roc} estimates the ROC curve for a fitted classifier (an object
#' of class \code{svmroc}) using a new testing set.
#'
#' @param object An object of class \code{svmroc}.
#' @param new_X A new X matrix.
#' @param new_A A new set of class labels.
#'
#' @return An object of class \code{roc}, a list with components: \code{sens},
#'   a vector of sensitivity values, and \code{spec}, a vector of specificity values.
#'
#' @export
fit_roc = function(object, new_X, new_A) {

  # coerce new_A to factor
  new_A = as.factor(new_A)

  # check that levels(new_A) match levels used in the original fit
  if (!(identical(levels(as.factor(new_A)), levels(as.factor(object$new_A))))) {
    stop("'new_A' must have the same levels as the class labels in the original fit.")
  }

  # check that new_X has the same number of columns as used in the original fit
  if (ncol(new_X) != ncol(object$new_X)) {
    stop("'new_X' must have the same number of columns as the covariate matrix used to obtain 'object'")
  }

  # create vectors to store sensitivity and specificity of ROC curve
  se = rep(NA, length(object$weights))
  sp = rep(NA, length(object$weights))

  # loop through weights
  for (w in 1:length(object$weights)) {

    # get current weight
    weight = object$weights[w]

    # predict treatment assignments
    pred_vals = predict(object = object, newdata = new_X, weight = weight)

    # calculate sensitivity and specificity
    spec = sum(new_A == levels(new_A)[1] & pred_vals == levels(new_A)[1]) / sum(new_A == levels(new_A)[1])
    sens = sum(new_A == levels(new_A)[2] & pred_vals == levels(new_A)[2]) / sum(new_A == levels(new_A)[2])

    # save sensitivity and specificity
    se[w] = sens
    sp[w] = spec

  }  # end loop through weights

  # creates vectors of sensitivity and specificity to return
  to_return = list(sens = se, spec = sp)
  class(to_return) = "roc"

  return(to_return)

}  # end function fit_roc

# function to compute area under the curve for an object of class roc
# calls the generic directly, so there is no need for documentation
#' @export
auc.roc = function(object) {

  # calculate AUC based on sensitivity and specificity vectors returned by fit_roc
  tpf = object$sens
  fpf = 1 - object$spec
  if (fpf[length(fpf)] < fpf[1]) fpf = rev(fpf)
  if (tpf[length(tpf)] < tpf[1]) tpf = rev(tpf)
  fpf = c(0, fpf, 1); tpf = c(0, tpf, 1)
  idx = 2:length(fpf)

  return(abs(as.double( (fpf[idx] - fpf[idx - 1]) %*% (tpf[idx] + tpf[idx - 1])) / 2))

}  # end function auc.roc

#' Estimate ROC curve for fitted classifier with linear interpolation
#'
#' \code{true_roc} allows for plotting an approximation to the true ROC curve along with
#'   confidence bands.
#'
#' @param object An object of class svmroc
#' @param new_X New X matrix to determine true ROC curve.
#' @param new_A New class assignments to determine true ROC curve.
#' @param points Points at which to interpolate the ROC curve.
#'
#' @return A list with components: \code{x} and \code{y} giving the coordinates of
#'   the ROC curve.
#'
#' @export
true_roc = function(object, new_X, new_A, points) {

  # check that levels(new_A) match levels used in the original fit
  if (!(identical(levels(as.factor(new_A)), levels(as.factor(object$new_A))))) {
    stop("'new_A' must have the same levels as the class labels in the original fit.")
  }

  # check that new_X has the same number of columns as used in the original fit
  if (ncol(new_X) != ncol(object$new_X)) {
    stop("'new_X' must have the same number of columns as the covariate matrix used to obtain 'object'")
  }

  # estimate true ROC curve
  roc = fit_roc(object = object, new_X = new_X, new_A = new_A)

  # linear interpolation to get true ROC curve
  true_roc = approx(x = 1 - roc$spec, y = roc$sens, xout = points, yleft = 0, yright = 1)$y

  # create object to return
  to_return = list(x = points, y = true_roc)

  return(to_return)

}  # end true_roc
