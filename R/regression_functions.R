fit_lasso_poly <- function(X, y, degree, lambda = NULL) {
  #' Fit a lasso regression using quadratic polynomial basis, with interactions.
  #'
  #' Compute regression function and derivative estimates based
  #' on polynomial basis lasso with penalty parameter chosen by
  #'  cross validation (CV).
  #'
  #' @param X matrix of covariates.
  #' @param y vector of responses.
  #' @param degree maximum degree of polynomial terms.
  #' @param lambda optional scalar tuning parameter, if "NULL" chosen via
  #'     cross-validation.
  #'
  #' @return List containing: A function "fit" which takes matrix input of the
  #'     same width as X, and returns a vector of y-predictions. A
  #'     scalar "lambda" the tuning parameter.

  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package \"glmnet\" needed for this function to work.
         Please install it.",
         call. = FALSE)
  }

  if (degree==1){
    basis <- X
  }
  else{
    basis <- as.matrix(stats::poly(X,degree = degree,raw=TRUE))
  }

  # Fit model
  if (is.null(lambda)){
    fitted_model <- glmnet::cv.glmnet(basis,y)
    lambda <- fitted_model$lambda.1se
  } else{
    fitted_model <- glmnet::glmnet(basis, y, lambda = lambda)
  }

  if (degree==1){
    fit <- function(X){
      basis <- X
      return(as.vector(stats::predict(fitted_model, basis)))
    }
  } else{
    fit <- function(X){
      basis <- as.matrix(stats::poly(X,degree = degree,raw=TRUE))
      return(as.vector(stats::predict(fitted_model, basis)))
    }
  }

  return(list("fit"=fit,
              "lambda"=lambda))
}



fit_xgboost <- function(X, y, params, derivative = FALSE) {
  #' Fit pre-tuned XGBoost regression for use in simulations.
  #'
  #' @param X matrix of covariates.
  #' @param y vector of responses.
  #' @param params XGBoost hyperparameters.
  #' @param derivative logical determining if numerical difference derivative
  #'     estimate (wrt the first predictor) should also be returned.
  #'
  #' @return list containing a function "fit" which takes matrix input of the
  #' same width as X, and returns a vector of predictions. Optionally the list
  #' also contains a function "deriv_fit" for numerical difference derivative
  #' estimates.

  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("Package \"xgboost\" needed for this function to work.
         Please install it.",
         call. = FALSE)
  }

  nrounds <- params$nrounds
  params$nrounds <- NULL

  xgb <- xgboost::xgboost(data = X,
                          label = y,
                          params = params,
                          nrounds = nrounds,
                          verbose = 0,
                          nthread = 2,
                          objective = "reg:squarederror")

  fit <- function(X){
    stats::predict(xgb, newdata = as.matrix(X))
  }

  if (derivative){

    default <- stats::sd(X[,1])/4

    deriv_fit <- function(X, D=NULL){
      if (is.null(D)){
        D <- default
      }
      n <- dim(X)[1]
      newcovariates <- rbind(X, X)

      newcovariates[,1] <- rep(X[,1], 2) + rep(c(D/2, -D/2), each=n)

      newpred <- stats::predict(xgb, newdata = as.matrix(newcovariates))

      return((newpred[1:n] - newpred[(n+1):(2*n)])/D)
    }

    return(list("fit"=fit,"deriv_fit"=deriv_fit))

  } else {
    return(list("fit"=fit))
  }
}
