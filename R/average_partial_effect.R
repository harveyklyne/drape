drape <- function(y, x, z,
                  response_regression,
                  predictor_regression,
                  resmooth_bw = NULL,
                  spline_df = NULL,
                  nfolds = 5L,
                  foldid = NULL,
                  verbose = FALSE) {
  #' Estimate the doubly-robust average partial effect estimate of X on Y, in the presence of Z.
  #'
  #' @param y vector of responses.
  #' @param x vector of the predictor of interest.
  #' @param z matrix of additional predictors.
  #' @param response_regression function which takes input data of the form
  #'     (X,y), where X=cbind(x,z), and returns a prediction function f:X -> y and optionally a
  #'     similar derivative estimation function (in this case no resmoothing is done).
  #' @param predictor_regression function which takes input data of the form
  #'     (z,x), and returns a prediction function m:z -> x.
  #' @param resmooth_bw optional numeric to be used as resmoothing bandwidth,
  #'     otherwise chosen via cross-validation. Only used if
  #'     response_regression doesn't predict derivatives.
  #' @param spline_df optional double, a smoothing parameter for the
  #'     unconditional spline score estimator, corresponding to the effective degrees
  #'     of freedom for a smoothing spline. If NULL, chosen via
  #'     cross-validation.
  #' @param nfolds integer, number of sample-splits. If set to
  #'     one, then all data is used for both training and evaluation.
  #' @param foldid optional vector with components in 1:nfolds indicating the folds in which
  #'      each observation fell. Overwrites nfolds.
  #' @param verbose boolean controlling level of information outputted.
  #'
  #' @return list containing the average partial effect estimate and the
  #'  corresponding standard error estimate. If verbose=TRUE, additionally contains
  #'  variables used in computations.
  #'
  #' @export
  #'
  #' @examples
  #' set.seed(0)
  #' data <- simulate_data(200, "normal", "plm")
  #' response_regression <- function(X,y){
  #'     df <- data.frame(y,X)
  #'     colnames(df) <- c("y", paste0("X", 1:10))
  #'     lm1 <- stats::lm(y~X1+sin(X2), data=df)
  #'     fit <- function(newX){
  #'         newdf <- data.frame(newX)
  #'         colnames(newdf) <- paste0("X", 1:10)
  #'         return(as.vector(stats::predict(lm1, newdata=newdf)))}
  #'     return(list("fit"=fit))
  #' }
  #' predictor_regression <- function(z,x){
  #'     df <- data.frame(x,z)
  #'     colnames(df) <- c("x", paste0("Z", 1:9))
  #'     lm1 <- stats::lm(x~Z1+Z2, data=df)
  #'     fit <- function(newz){
  #'         newdf <- data.frame(newz)
  #'         colnames(newdf) <- paste0("Z", 1:9)
  #'         return(as.vector(stats::predict(lm1, newdata=newdf)))}
  #'     return(list("fit"=fit))
  #' }
  #' drape(data$y, data$x, data$z, response_regression, predictor_regression, nfolds=2)

  if(!is.vector(y)){
    warning("y must be a vector. Updating y <- as.vector(y).")
    y <- as.vector(y)
  }

  # Converting (x,z) to matrix form for compatibility.
  X <- cbind(x,z)
  d <- 1

  n <- length(y)

  if(dim(X)[1]!=n){stop("Must have length(y) == dim(X)[1].")}

  if (is.null(foldid)){
    foldid <- sample(rep(seq(nfolds), length.out=n))
  } else{ nfolds <- max(foldid) }

  if(nfolds < 1){
    warning("nfolds must be a positive integer. Setting nfolds=1.")
    nfolds <- 1L
  }
  if(round(nfolds)!=nfolds){
    warning("nfolds must be an integer. Rounding.")
    nfolds <- round(nfolds)
  }
  foldid <- sample(rep(seq(nfolds), length.out=n))

  if (verbose){out_all <- list("folds" = foldid)}


  # Train regression procedures

  regression_list <- list()
  predictor_list <- list()

  for (fold in seq(nfolds)){

    te <- (foldid == fold)
    if (nfolds > 1){tr <- !te} else{tr <- te}
    train <- list("X" = matrix(X[tr,], ncol=ncol(X)), "y" = y[tr])

    regression_list[[fold]] <- response_regression(train$X, train$y)
    predictor_list[[fold]] <- predictor_regression(train$X[,-d], train$X[,d])

  }

  # See if regression method produces usable derivative estimates on all folds
  do_resmooth <- !all(sapply(regression_list, function(x){"deriv_fit" %in% names(x)}))

  # Choose resmoothing bandwidth
  if (do_resmooth & is.null(resmooth_bw)){
    cv_resmooth_output <- cv_resmooth(X=X,
                                      y=y,
                                      d=1,
                                      regression=regression_list,
                                      tol=2,
                                      prefit=TRUE,
                                      foldid=foldid)
    resmooth_bw <- cv_resmooth_output$bw_opt[1]
    if (verbose) {out_all["cv_resmooth_output"] = list(cv_resmooth_output)}
  }

  # Choose spline df (using in-sample residuals from whole data set)
  if (is.null(spline_df)){
    m_fit <- predictor_regression(X[,-d], X[,d])$fit

    # Variance estimation
    resid <- X[,d] - m_fit(X[,-d])
    tree <- partykit::ctree(resid^2 ~ .,
                            data=data.frame("resid"=resid, "z"=X[,-d]))
    sigma <- sqrt(unname(stats::predict(tree,
                                        newdata = data.frame("z"=X[,-d]))))
    eps <- resid / sigma
    cv_spline_score_output <- cv_spline_score(x=eps)
    spline_df <- cv_spline_score_output$df_1se # Or df_min

    if (verbose) {out_all["cv_spline_score_output"] = list(cv_spline_score_output)}
  }

  # Evaluate f, df, and score
  f_pred <- rep(NA, n)
  df_pred <- rep(NA, n)
  score_pred <- rep(NA, n)
  for (fold in seq(nfolds)){

    if (verbose) {out_fold = list()}

    te <- (foldid == fold)
    if (nfolds > 1){tr <- !te} else{tr <- te}
    train <- list("X" = matrix(X[tr,], ncol=ncol(X)), "y" = y[tr])
    test <- list("X" = matrix(X[te,], ncol=ncol(X)), "y" = y[te])


    # # # f and df evaluation
    f_fit <- regression_list[[fold]]$fit

    if ("deriv_fit" %in% names(regression_list[[fold]])) {
      f_pred[te] <- f_fit(test$X)
      df_fit <- regression_list[[fold]]$deriv_fit
      df_pred[te] <- df_fit(test$X)
    }  else{
      smooth_list <- resmooth(fit = f_fit,
                              X = test$X,
                              d = d,
                              bw = resmooth_bw)
      if (verbose){out_fold["smooth_list"] <- list(smooth_list)}
      f_pred[te] <- smooth_list$pred[[1]]
      df_pred[te] <- smooth_list$deriv[[1]]
    }

    # # # Score evaluation
    m_fit <- predictor_list[[fold]]$fit

    # Variance estimation
    resid_train <- train$X[,d] - m_fit(train$X[,-d])
    tree <- partykit::ctree(resid_train^2 ~ .,
                            data=data.frame("resid_train"=resid_train,
                                            "z"=train$X[,-d]))
    if (verbose){out_fold["tree"] <- list(tree)}
    sigma_train <- sqrt(unname(stats::predict(tree,
                                              newdata = data.frame("z"=train$X[,-d]))))

    resid_test <- test$X[,1] - m_fit(test$X[,-1])
    sigma_test <- sqrt(unname(stats::predict(tree,
                                             newdata = data.frame("z"=test$X[,-d]))))

    # If any variance estimates are close to zero, assume homogeneous (more stable)
    if (any(sigma_test < 0.01)){
      warning(paste0("Minimal heterogeneous variance equals ",
                     round(min(sigma_test),4), ", reducing to homogeneous
                   case for stability."))
      sigma_train <- rep(sqrt(mean(resid_train^2)), length(train$y))
      sigma_test <- rep(sqrt(mean(resid_test^2)), length(test$y))
    }

    # Score estimation
    eps_train <- resid_train / sigma_train
    rho <- spline_score(eps_train, df=spline_df)$rho
    if (verbose){out_fold["rho"] <- list(rho)}

    eps_test <- resid_test / sigma_test
    rho_test <- rho(eps_test)

    # If any score estimates are large, assume Gaussian (more stable)
    if (any(abs(rho_test) > 10)){
      warning(paste0("Maximal absolute spline score prediction equals ",
                     round(max(abs(rho_test)),2), ", reducing to Gaussian
                   case for stability."))
      rho_test <- - eps_test
    }

    score_pred[te] <- rho_test/sigma_test # rho_P(x,z)=1/sigma(z) * rho_\varepsilon( (x-m(z))/sigma(z) )

    if (verbose){
      str <- paste0("fold", fold)
      out_all[[str]] <- out_fold
    }

  } # End of fold

  evaluations <- df_pred - score_pred * (y - f_pred)
  est <- mean(evaluations)
  se <- stats::sd(evaluations) / sqrt(n)

  out <- list("est" = est,
              "se" = se)

  if (verbose) {
    out[["computations"]] <- out_all
    out[["f_pred"]] <- f_pred
    out[["df_deriv"]] <- df_pred
    out[["score_pred"]] <- score_pred
    out[["evaluations"]] <- evaluations
  }

  return(out)
}
