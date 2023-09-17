resmooth <- function(fit,
                     X,
                     d=1,
                     bw = exp(seq(-1,1))/(2*sqrt(3))*stats::sd(X[,d]),
                     n_points=101,
                     sd_trim=5){
  #' Resmooth the predictions of a fitted model
  #'
  #' Smooth the predictions of a fitted model by convolving them
  #' with a Gaussian kernel along the d'th covariate.
  #'
  #' @param fit a prediction function taking matrix input and
  #'     returning a vector.
  #' @param X matrix of covariates.
  #' @param d integer index of covariate to be smoothed along.
  #' @param bw vector of bandwidths for the Gaussian kernel.
  #' @param n_points integer number of gridpoints to be used for convolution.
  #' @param sd_trim float number of standard deviations at which to trim the
  #'     Gaussian distribution.
  #'
  #' @return List with the following elements. A list "pred" of the same length
  #'     as "bw". Each element is a vector of predictions which are smooth with
  #'     respect to the dth column of X, with smoothness corresponding to the
  #'     respective element of "bw". A similar list "deriv" of corresponding
  #'     vectors of first derivatives. Vectors "gridpoints" and "prob_weights"
  #'     of equally spaced gridpoints and corresponding normal density weights.
  #'     Vector "bw" of bandwidths used.
  #' @export
  #'
  #' @examples
  #' # Single bandwidth
  #' X <- matrix(seq(-2,2,by=0.05))
  #' fit <- function(Y){1*(rowMeans(Y)<0)}
  #' sm <- resmooth(fit=fit, X=X, d=1, bw=0.2)
  #' sm$pred[[1]]
  #'
  #' # Multiple bandwidths simultaneously
  #' X <- matrix(stats::rnorm(200), ncol=2)
  #' y <- X[,1] + sin(X[,2]) + 0.5 * stats::rnorm(nrow(X))
  #' df <- data.frame(y,X)
  #' lm1 <- stats::lm(y~X1+sin(X2), data=df)
  #' fit <- function(Y){as.vector(stats::predict(lm1, newdata=data.frame(Y)))}
  #' resmooth(fit=fit, X=X, d=2)

  gridpoints <- seq(-sd_trim, sd_trim, length.out=n_points)
  prob_weights <- stats::dnorm(gridpoints)
  prob_weights <- prob_weights / sum(prob_weights)

  n <- dim(X)[1]
  nbw <- length(bw)

  newX <- new_X(X=X, d=d, MC_variates=gridpoints, bw=bw)
  pred <- fit(newX)
  dpred <- rep(gridpoints,times=n*nbw)*pred/rep(bw,each=n_points*n)

  fhat <- MC_sums(pred*prob_weights, n=n, nMC=n_points, nbw=nbw)
  dfhat <- MC_sums(dpred*prob_weights, n=n, nMC=n_points, nbw=nbw)

  return(list('pred'=fhat,
              'deriv'=dfhat,
              'gridpoints'=gridpoints,
              'prob_weights'=prob_weights,
              'bw'=bw
  ))

}


MC_sums <- function(a, n, nMC, nbw){
  #' Compute sums of a Monte Carlo vector for use in resmoothing.
  #'
  #' @param a vector of length (n x nMC x nbw).
  #' @param n integer.
  #' @param nMC integer.
  #' @param nbw integer.
  #'
  #' @return list with nbw elements. The j'th element of which is
  #' a vector of length n, the i'th element being the sum of the
  #' (((j-1)n + (i-1)) x nMC + 1) to (((j-1)n + i) x nMC)
  #' elements of a inclusive.

  if (length(a)!=n*nMC*nbw){stop("MC_sums requires length(a)=n*nMC*nbw.")}
  vec <- .colSums(a, nMC, n*nbw)
  # vec = (fhat_1(X1),...,fhat_1(Xn),...,fhat_nbw(X1),...,fhat_nbw(Xn))
  # length(n*nbw)

  return(unname(split(vec, rep(1:nbw, each=n))))
  # split into nbw-list of n-vectors (fhat_j(X1),...,fhat_j(Xn))
}

new_X <- function(X, d, MC_variates, bw){
  #' Generate a matrix of covariates for use in resmoothing, in which
  #' the d'th column of X is translated successively by the Kronecker
  #' product of bw and MC_variates.
  #'
  #' @param X matrix of covariates.
  #' @param d integer index of covariate to be smoothed along.
  #' @param MC_variates vector of standard Gaussian rvs.
  #' @param bw vector of bandwidths for the Gaussian kernel.
  #'
  #' @return matrix with ncol(X) columns and (nrow(X)length(MC_variates)
  #' length(bw)) rows.

  if(!is.matrix(X)){
    warning("new_X requires X to be a matrix. Updating X <- matrix(X).")
    X <- matrix(X)
  }

  n <- dim(X)[1]
  nMC <- length(MC_variates)
  nbw <- length(bw)

  d_col <- rep(X[,d],each=nMC) + as.vector(rep(bw,each=n) %x% MC_variates)
  rest <- X[rep(rep(1:n, each=nMC), times=nbw),-d]

  out <- unname(as.matrix(cbind(rest, d_col))) # need to reorder columns
  p <- dim(X)[2]
  ind <- append(seq_len(p-1), p, d-1) # (1:d-1, p, d:p-1), works with d=0 and p
  return(as.matrix(out[, ind]))
}


cv_resmooth <- function(X, y, d=1,
                        regression,
                        tol=2,
                        prefit=FALSE,
                        foldid=NULL,
                        bw=exp(seq(-5, 2, 0.2))/(2*sqrt(3)) * stats::sd(X[,d]),
                        nfolds=5L,
                        n_points=101,
                        sd_trim=5){
  #' K-fold cross-validation for resmoothing bandwidth.
  #'
  #' Picks the largest resmoothing bandwidth achieving a cross-validation score within some specified tolerance of the original regression.
  #'
  #' @param X matrix of covariates.
  #' @param y vector of responses.
  #' @param d integer index of covariate to be smoothed along.
  #' @param regression If prefit = FALSE this is a function which takes input data of the form (X,y),
  #'     and returns a prediction function. This prediction function
  #'     itself accepts matrix input same width as X,
  #'     and returns a vector of y-predictions,
  #'     and optionally a vector of derivative predictions.
  #'     If prefit = TRUE then this is a list of length nfolds with each entry containing a component "fit" consisting
  #'     of a prediction function taking matrix input and returning a vector.
  #' @param tol vector of tolerances controlling the degree of permissible cross-validation error increase.
  #'     Larger values lead to a larger amount of smoothing being selected.
  #' @param prefit boolean signifying if the regressions are already fit to the training data for each fold.
  #' @param foldid optional vector with components in 1:nfolds indicating the folds in which
  #'      each observation fell. Overwrites nfolds.
  #' @param bw vector of bandwidths for the Gaussian resmoothing kernel.
  #' @param nfolds integer number of cross-validation folds.
  #' @param n_points integer number of gridpoints to be used for convolution.
  #' @param sd_trim float number of standard deviations at which to trim the
  #'     Gaussian distribution.
  #'
  #' @return list. Vector "bw" of bandwidths used. Vectors "cv" of
  #'     cross-validation scores and numeric "cv_unsm" for the cross-validation
  #'     without any smoothing. Vector "bw_opt_inds" for the indices of the selected bandwidths
  #'     under various tolerances. Vector "bw_opt" for the corresponding bandwidths.
  #' @export
  #'
  #' @examples
  #' X <- matrix(stats::rnorm(200), ncol=2)
  #' y <- X[,1] + sin(X[,2]) + 0.5 * stats::rnorm(nrow(X))
  #' reg <- function(X,y){
  #'     df <- data.frame(y,X)
  #'     colnames(df) <- c("y", "X1", "X2")
  #'     lm1 <- stats::lm(y~X1+sin(X2), data=df)
  #'     fit <- function(newX){
  #'         newdf = data.frame(newX)
  #'         colnames(newdf) <- c("X1", "X2")
  #'         return(as.vector(stats::predict(lm1, newdata=newdf)))}
  #'     return(list("fit"=fit))
  #' }
  #' cv_resmooth(X=X, y=y, d=2, regression=reg, tol = c(0.5, 1, 2))

  macheps <- 10^(-8)
  n <- dim(X)[1]
  nbw <- length(bw)

  if (is.null(foldid)){
    foldid <- sample(rep(seq(nfolds), length.out=n))
  } else{ nfolds <- max(foldid) }

  # Matrix to hold MSE scores for the desired bandwidths.
  mse_bws <- matrix(NA,nrow=n,ncol=nbw)
  # Vector to hold MSE scores for the unsmoothed regression.
  mse_unsm <- rep(NA,n)

  for (fold in seq(nfolds)){

    which <- (foldid == fold)

    if (prefit){fit <- regression[[fold]]$fit} else{fit <- regression(as.matrix(X[!which,]), y[!which])$fit}

    # Fold CV score for unsmoothed regression
    pred <- fit(as.matrix(X[which,]))
    mse_unsm[which] <- (y[which]-pred)^2

    # Resmoothed fold CV scores by bandwidth.
    pred_list <- resmooth(fit=fit,
                          X=as.matrix(X[which,]),
                          d=d,
                          bw=bw,
                          n_points=n_points,
                          sd_trim=sd_trim)$pred
    foo <- function(x){(y[which]-x)^2}
    mse_bws[which,] <- sapply(pred_list, foo)
  }

  cv_unsm <- mean(mse_unsm)
  cv <- colMeans(mse_bws)

  cv_all <- c(cv, cv_unsm)

  cv_min_ind <- which.min(cv_all)
  cv_min <- cv_all[cv_min_ind]

  if (cv_min_ind < length(cv_all)) {
    # minimiser is among bandwidths
    # we end up considering only bandwidths larger than the minimiser
    sd_diff <- apply(mse_bws - mse_bws[, cv_min_ind], 2, stats::sd) / sqrt(n)
    # Have sd_diff[cv_min_ind] = 0
  } else {
    # minimiser is cv_unsm
    sd_diff <- apply(mse_bws - mse_unsm, 2, stats::sd) / sqrt(n)
  }
  thresh <- outer(sd_diff, tol, "*") + cv_min + macheps

  bw_opt_inds <- integer(ncol(thresh))
  for (j in 1L:ncol(thresh)) {
    bw_opt_inds[j] <- max(which(cv <= thresh[, j]))
    # If minimiser is among bandwidths we always have cv_min <= thresh, so bw_opt_inds >= cv_min_ind.
  }
  names(bw_opt_inds) <- as.character(tol)

  # bw_opt_inds_raw <- bw_opt_inds
  if (any(is.infinite(bw_opt_inds))){
    warning("No bandwidth has CV score within tolerance.")
    bw_opt_inds[is.infinite(bw_opt_inds)] <- 1L # take minimium bandwidth
  }

  return(list('bw'=bw,
              'cv'=cv,
              'cv_unsm'=cv_unsm,
              'bw_opt_inds'=bw_opt_inds,
              'bw_opt'=bw[bw_opt_inds]))
}
