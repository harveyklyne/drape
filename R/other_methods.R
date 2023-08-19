# Partially linear model

partially_linear = function(X, y, g_params, m_params){
  #' Fit a doubly-robust partially linear regression using the DoubleML package
  #' and pre-tuned XGBoost regressions, for use in simulations.
  #'
  #' @param X matrix of covariates.
  #' @param y vector of responses.
  #' @param g_params XGBoost hyperparameters for partially linear regression of y on X.
  #' @param m_params XGBoost hyperparameters for predictor regression of the first column of X on the others.
  #'
  #' @return List containing the linear parameter estimate and the
  #'  corresponding standard error estimate.

  if (!requireNamespace(c("DoubleML", "paradox"), quietly = TRUE)) {
    stop("Packages \"DoubleML\" and \"paradox\" needed for this function to work. Please install them.",
         call. = FALSE)
  }

  obj_dml_data = DoubleML::double_ml_data_from_matrix(X=X[,-1], y=y, d=X[,1])


  learner_g <- mlr3::lrn("regr.xgboost")
  names(g_params)[names(g_params)=="max.depth"] <- "max_depth"
  learner_g$param_set$values <- g_params

  learner_m <- mlr3::lrn("regr.xgboost")
  names(m_params)[names(m_params)=="max.depth"] <- "max_depth"
  learner_m$param_set$values <- m_params


  doubleml = DoubleML::DoubleMLPLR$new(obj_dml_data,
                             ml_m = learner_m,
                             ml_g = learner_g,
                             n_folds = 5,
                             apply_cross_fitting=TRUE)

  doubleml$fit()
  est <- as.vector(doubleml$all_coef)
  se <- as.vector(doubleml$all_se)

  return(list("est" = est, "se" = se))

}


# Rothenhausler + Yu

rothenhausler_yu <- function(X, y, f_lambda, m_lambda){
  #' Estimate the average partial effect of using the debiased lasso method of Rothenhausler and Yu,
  #'  using pre-tuned lasso penalties, for use in simulations.
  #'
  #' @param X matrix of covariates.
  #' @param y vector of responses.
  #' @param f_lambda lasso penalty for regression of y on X.
  #' @param m_lambda lasso penalty for predictor regression of the first column of X on the others.
  #'
  #' @return List containing the linear parameter estimate and the
  #'  corresponding standard error estimate.

  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package \"glmnet\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  basis_list <- rothenhausler_basis(X)
  mod_basis <- basis_list$mod_basis
  dbasis <- basis_list$dbasis

  f_fit <- glmnet::glmnet(x=mod_basis, y=y, lambda=f_lambda)
  m_fit <- glmnet::glmnet(x=as.matrix(mod_basis[,-1]), y=mod_basis[,1], lambda=m_lambda)

  zed <- mod_basis[,1] - stats::predict(m_fit, mod_basis[,-1])
  eps <- y - stats::predict(f_fit, mod_basis)

  # debiased lasso estimate
  est <- f_fit$beta[1] - sum(zed * eps) / sum(zed * mod_basis[,1])
  # standard error estimate
  var_vec <- as.vector(zed * eps / mean(zed * mod_basis[,1]) + dbasis %*% f_fit$beta)
  se <- stats::sd(var_vec) / sqrt(length(y))

  return(list("est" = est, "se" = se))

}


rothenhausler_basis <- function(X){
  #' Generate the modified quadratic basis of Rothenhausler and Yu.
  #'
  #' @param X matrix of covariates.
  #'
  #' @return List containing the modified basis matrices for regression and derivative estimation.

  p <- ncol(X)
  n <- nrow(X)

  basis <- stats::poly(X,degree = 2,raw=TRUE)
  dbasis <- matrix(0, nrow=n, ncol=ncol(basis))

  dbasis[,c(1:2, 1+(2:p)*(1 + 2:p)/2)] <- cbind(rep(1,nrow(X)), 2*X[,1], X[,-1])
  # bhat(x)[1] = b(x)[1] = x[1]
  # bhat(x)[-1] = b(x)[-1] - x[1] * E[b'(x)[-1]]
  mod_basis <- matrix(0, nrow=n, ncol=ncol(basis))
  mod_basis[,1] <- X[,1]
  mod_basis[,-1] <- basis[,-1] - X[,1] * rep(colMeans(dbasis)[-1], each=n)

  return(list("mod_basis" = mod_basis,
              "dbasis" = dbasis))

}

