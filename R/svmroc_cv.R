# svmroc_cv.R

# A function to estimate the SVM ROC curve using cross-validation

#' Estimate SVM ROC curve using cross-validation
#'
#' \code{svmroc_cv} estimates the SVM ROC curve using cross-validation.
#'
#' @param X A matrix of covariates used for fitting the SVM.
#' @param A A factor with two levels or an object that can be coerced to a factor with
#'   two levels; used as the response when fitting the SVM.
#' @param kernel The kernel used for fitting the SVM, either "linear" or "Gaussian".
#'   Defaults to "linear".
#' @param fold_id A vector indicating which fold each observation belows to when estimating
#'   the ROC curve. If \code{NULL}, folds are generated randomly. Defaults to \code{NULL}
#' @param num_folds_roc Number of folds used to estimate ROC curve. Ignored if \code{fold_id}
#'   is not \code{NULL}. Defaults to 10.
#' @param lambdas A numeric vector of penalty parameters to select from. Defaults to
#'   \code{c(0.1, 0.25, 0.5, 1)}.
#' @param sigmas A numeric vector of bandwidth parameters to select from. Defaults to
#'   \code{c(0.1, 0.25, 0.5, 1)}.
#' @param num_folds_params Number of folds for cross-validation to select tuning parameters.
#'   Defaults to 10.
#' @param weights A vector of weights to use for the weighted SVM. Defaults to
#'   \code{seq(0.01, 0.99, 0.01)}.
#' @param seed Random number seed to set. If \code{NULL}, no seed is set. Defaults to \code{NULL}.
#'
#' @return An object of class \code{svmroc}, a list with only the following components:
#'   \code{sens}, a vector of sensitivities across weights, averaged across folds;
#'   \code{spec}, a vector of specificities across weights, averaged across folds.
#'
#' @export
svmroc_cv = function(X, A, kernel = "linear", fold_id = NULL, num_folds_roc = 10,
                     lambdas = c(0.1, 0.25, 0.5, 1), sigmas = c(0.1, 0.25, 0.5, 1),
                     num_folds_params = 10, weights = seq(0.01, 0.99, 0.01), seed = NULL) {

  # if seed is supplied, set it here
  if (!is.null(seed)) set.seed(seed)

  # if fold ID is not provided, create fold ID
  if (is.null(fold_id)) {
    if (!requireNamespace("caret", quietly = TRUE)) {
      stop("Package \"caret\" is needed to create folds. Please install it.")
    }
    fold_id = caret::createFolds(y = A, k = num_folds_roc, list = F)
  }  # end if fold_id is NULL

  # create objects to save sensitivities and specificities across folds
  all_sens = matrix(NA, nrow = length(weights), ncol = length(unique(fold_id)))
  all_spec = matrix(NA, nrow = length(weights), ncol = length(unique(fold_id)))

  # coerce data to matrices
  X = as.matrix(X)
  A = as.matrix(A)

  # loop through folds
  for (fold in unique(fold_id)) {

    # create training data
    X_train = X[-which(fold_id == fold), ]
    A_train = A[-which(fold_id == fold), ]

    # create testing data
    X_test = X[which(fold_id == fold), ]
    A_test = A[which(fold_id == fold), ]

    # SVM ROC on current fold
    fold_out = svmroc(X = X_train, A = A_train, new_X = X_test, new_A = A_test, kernel = kernel,
                      lambdas = lambdas, sigmas = sigmas, num_folds = num_folds_params,
                      weights = weights, seed = NULL)

    # save sensitivities and specificities for current fold
    all_sens[ , fold] = fold_out$sens
    all_spec[ , fold] = fold_out$spec

  }  # end loop through folds

  # average sensitivities and specificities across folds
  avg_sens = apply(all_sens, 1, mean)
  avg_spec = apply(all_spec, 1, mean)

  # create object to return
  to_return = list(sens = avg_sens, spec = avg_spec)
  class(to_return) = "svmroc"

  return(to_return)

}  # end function svmroc_cv
