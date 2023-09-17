simulate_data <- function(n, ex_setting, f_setting){
  #' Generate simulation data.
  #'
  #' If ex_setting = "401k" then 401k data set is used for (X,Z). Otherwise:
  #'
  #' Z ~ N_{9}(0,Sigma),
  #' where Sigma_{jj} = 1, Sigma_{jk} = corr for all j not equal to k.
  #' X = m(Z) + s(Z)*ex
  #' where m and sigma are step functions of z_1 and z_3 respectively.
  #' Y = f(X,Z) + N(0,1)
  #'
  #' @param n integer number of samples. For "401k" ex_setting this is ignored
  #'     and the whole data set is used.
  #' @param ex_setting string "normal", "mixture2", "mixture3",
  #'     "logistic", "t4", "401k".
  #' @param f_setting string "plm", "additive", "interaction".
  #'
  #' @return list containing y, x, z. Additionally contains the population
  #'     nuisance parameters evaluated on the data, and the sample version of
  #'     the average partial effect.
  #' @export
  #'
  #' @examples
  #' simulate_data(100, "normal", "plm")

  if (!(ex_setting %in% c("normal", "mixture2", "mixture3", "logistic", "t4", "401k"))){
    stop("ex_setting not recognised")
  }
  if (!(f_setting %in% c("plm", "additive", "interaction"))){
    stop("f_setting not recognised")
  }

  if (ex_setting == "401k"){

    data <- readRDS("data-raw/401k/401k_data.rds") # Data file generated using data-raw/401k/data_preparation.R

    sc <- stats::sd(data$x) * 2/ sqrt(5)
    xmn <- mean(data$x)
    zsc <- stats::sd(data$z[,1])
    zmn <- mean(data$z[,1])

    n <- length(data$x)
    x <- data$x
    z <- data$z

    xscale <-  (x-xmn)/sc
    z2scale <- (z[,1]-zmn)/zsc

    m <- NULL
    rho <- NULL

  }

  if (ex_setting %in% c("normal", "mixture2", "mixture3", "logistic", "t4")){

    z <- z_correlated_normal(n, 9, corr=0.5)
    m <- 1*(z[,1]>0)
    sigma <- 1/sqrt(2) + (sqrt(3)-1)/sqrt(2)*(z[,3]<0)

    # eps_x
    if (ex_setting=="normal"){
      eps_x <- stats::rnorm(n)
      rho <- -eps_x
    }
    if (ex_setting=="mixture2"){
      eps_x <- rmixture(n, sd=1/sqrt(2))
      rho <- mixture_score(eps_x, sd=1/sqrt(2))
    }
    if (ex_setting=="mixture3"){
      eps_x <- rmixture(n, sd=1/sqrt(3))
      rho <- mixture_score(eps_x, sd=1/sqrt(3))
    }
    if (ex_setting=="logistic"){
      eps_x <- stats::rlogis(n, scale=sqrt(3)/pi)
      rho <- -pi/sqrt(3) * tanh(pi*eps_x/(2*sqrt(3)))
    }
    if (ex_setting=="t4"){
      eps_x <- 1/sqrt(2) * stats::rt(n=n, df=4)
      rho <- - 5/2 * eps_x / (1 + eps_x^2 / 2)
    }

    # x
    x <- m + sigma * eps_x

    sc <- 1
    xscale <-  x
    z2scale <- z[,1]

  }


  if (f_setting %in% c("plm","m")){
    f <- xscale + sigmoidal(z2scale, s=1) + sinosoidal(z2scale, a=1)
    f <- sc * f
    df <- rep(1,n)
  }
  if (f_setting=="additive"){
    f <- sigmoidal(xscale, s=1) + sinosoidal(xscale, a=1) + sinosoidal(z2scale, a=3)
    f <- sc * f
    df <- dsigmoidal(xscale, s=1) + dsinosoidal(xscale, a=1)
  }
  if (f_setting=="interaction"){
    f <- sigmoidal(xscale, s=3) + sinosoidal(xscale, a=3) + sinosoidal(z2scale, a=3) + xscale * z2scale
    f <- sc * f
    df <- dsigmoidal(xscale, s=3) + dsinosoidal(xscale, a=3) + z2scale
  }
  y <- f + sc * stats::rnorm(n)

  ground_truth_eval <- df - rho * (y - f)
  if (any(is.null(ground_truth_eval)) | (length(ground_truth_eval) == 0)){
    sample_theta <- mean(df)
    sample_std_dev <- NULL
  }
  else{
    sample_theta <- mean(ground_truth_eval)
    sample_std_dev <- stats::sd(ground_truth_eval)
  }

  return(list("y"=y, "x"=x, "z"=z, "sample_theta"=sample_theta, "sample_std_dev"=sample_std_dev,
              "f"=f, "df"=df, "rho"=rho, "m"=m))

}


z_correlated_normal <- function(n, p, corr){
  #' Generate n copies of Z ~ N_{p}(0,Sigma), where Sigma_{jj} = 1, Sigma_{jk} = corr for all j not equal to k.
  #'
  #' @param n integer number of samples.
  #' @param p integer number of dimensions.
  #' @param corr float correlation in (-1,1).
  #'
  #' @return n by p matrix.

  Sigma <- matrix(corr, nrow=p, ncol=p)
  diag(Sigma) <- 1
  sqrt_Sigma <- chol(Sigma)
  return(matrix(stats::rnorm(n*(p)), nrow=n) %*% sqrt_Sigma)
}

sigmoidal <- function(x,s){1/(1+exp(-s*x))}
dsigmoidal <- function(x,s=1){s*exp(-s*x)/(1+exp(-s*x))^2}

sinosoidal <- function(x,a){exp(-x^2/2)*sin(a*x)}
dsinosoidal <- function(x, a){(-x*sin(a*x)+a*cos(a*x))*exp(-x^2/2)}

rmixture <- function(n, sd){
  #' Symmetric mixture two Gaussian random variables.
  #'
  #' The resulting distribution is mean zero, variance one.
  #' X ~ N(-sqrt(1-sd^2),sd^2) wp 0.5,
  #'     N(sqrt(1-sd^2),sd^2) wp 0.5.
  #'
  #' @param n number of observations.
  #' @param sd standard deviation of each Gaussian.
  #'
  #' @return vector of length n

  u <- stats::runif(n)
  z <- stats::rnorm(n)
  # Gaussian mixture
  mu <- sqrt(1-sd^2)
  x <- (u<0.5)*(-mu + sd*z) + (u>=0.5)*(mu + sd*z)
  return(x)
}

mixture_score <- function(x, sd){
  #' Population score function for the symmetric mixture two Gaussian random variables.
  #'
  #' @param x vector of observations.
  #' @param sd standard deviation of each Gaussian.
  #'
  #' @return vector of length n
  mu <- sqrt(1-sd^2)
  p <- 0.5 * stats::dnorm(x,mean = -mu, sd=sd) + 0.5 * stats::dnorm(x, mean= mu, sd=sd)
  dp <- -1/2 * ( (x+mu)/(sd^2) * stats::dnorm(x, mean=-mu, sd=sd)  +
                   (x-mu)/(sd^2) * stats::dnorm(x, mean=mu, sd=sd))
  return(dp/p)
}


compare <- function(n, ex_setting, f_setting, nfold=5) {
  #' Generate simulation data and evaluate estimators, with sample splitting.
  #'
  #' @param n integer number of samples. For "401k" ex_setting this is ignored
  #'     and the whole data set is used.
  #' @param ex_setting string "normal", "mixture2", "mixture3",
  #'     "logistic", "t4", "401k".
  #' @param f_setting string "plm", "additive", "interaction".
  #' @param nfold integer number of cross-validation folds.
  #'
  #' @return list containing estimates, standard error estimates, and sample theta (for debugging).


  data <- simulate_data(n = n,
                        ex_setting = ex_setting,
                        f_setting = f_setting)

  y <- data$y
  X <- cbind(data$x, data$z)

  # Do sample splitting.
  foldid <- sample(rep(seq(nfold), length.out=n))

  out <- list()
  evaluate <- list()

  # Setup regressions
  tune <- rjson::fromJSON(file = paste0("data-raw/tuning/tune_", f_setting, "_", ex_setting, "_results.json"))

  ### Regression and derivative estimation
  resp_reg <- function(X,y){fit_xgboost(X = X,
                                        y = y,
                                        params = tune$xgb_f,
                                        derivative = TRUE)}

  regression_list <- vector("list", nfold)
  # Perform fits
  for (fold in seq(nfold)) {

    fold_index <- (foldid == fold)
    train <- list("X" = matrix(X[!fold_index,], ncol=ncol(X)),
                  # ensures correct dimension when sum(fold_index=1) and ncol(X)=1
                  "y" = y[!fold_index])
    regression_list[[fold]] <- resp_reg(train$X, train$y)
  }

  # Find bandwidths
  sm_bw_out <- cv_resmooth(X, y, d=1, regression=regression_list,
                           prefit = TRUE,
                           foldid = foldid)

  for (fold in seq(nfold)) {

    fold_index <- (foldid == fold)
    train <- list("X" = matrix(X[!fold_index,], ncol=ncol(X)),
                  "y" = y[!fold_index])
    test <- list("X" = matrix(X[fold_index,], ncol=ncol(X)),
                 # ensures correct dimension when sum(fold_index=1) and ncol(X)=1
                 "y" = y[fold_index])

    out[[fold]] <- compare_evaluate(train = train,
                                    test = test,
                                    ex_setting = ex_setting,
                                    f_setting = f_setting,
                                    regression = regression_list[[fold]],
                                    sm_bw_out = sm_bw_out)

    expand <- expand.grid(out[[fold]]$f, out[[fold]]$score)
    evaluate[[fold]] <- mapply(function(df, f, score){df - score * (test$y-f)},
                               df = rep(out[[fold]]$df, length(out[[fold]]$score)),
                               f = expand$Var1,
                               score = expand$Var2)


  } # End of fold

  estimate <- matrix(NA, nrow=n, ncol=length(out[[1]]$f)*length(out[[1]]$score))
  for (i in seq(nfold)){
    estimate[foldid==i,] <- evaluate[[i]]
  }

  est <- colMeans(estimate)
  se <- apply(estimate, 2, stats::sd) / sqrt(n)

  bw_names <- names(out[[1]]$df)
  score_names <- c("basis", "min", "1se")

  names(est) <- names(se) <- paste(paste0("deriv_", rep(bw_names,times=length(score_names))),
                                   paste0("score_",rep(score_names, each=length(bw_names))), sep="_")

  results <- list("est" = est,
                  "se" = se,
                  "sample_theta" = data$sample_theta,
                  "sample_std_err" = data$sample_std_dev / sqrt(n))

  return(results)

}


# needs to accept bw argument

compare_evaluate <- function(train, test, ex_setting, f_setting, regression, sm_bw_out){
  #' Evaluate estimators by training nuisance functions on training set and evaluating them on test set.
  #'
  #' @param train list containing vector of responses y and matrix of predictors X = (x,z).
  #' @param test list containing vector of responses y and matrix of predictors X = (x,z).
  #' @param ex_setting string "normal", "mixture2", "mixture3",
  #'     "logistic", "t4", "401k".
  #' @param f_setting string "plm", "additive", "interaction".
  #' @param regression Optional fitted regression.
  #' @param sm_bw_out Output of cv_resmooth.

  #' @return list containing f, df, and score estimates evaluated on the test set.

  tune <- rjson::fromJSON(file = paste0("data-raw/tuning/tune_", f_setting, "_", ex_setting, "_results.json"))

  if (missing(regression)) {
    ### Regression and derivative estimation
    resp_reg <- function(X,y){fit_xgboost(X = X,
                                          y = y,
                                          params = tune$xgb_f,
                                          derivative = TRUE)}

    regression <- resp_reg(train$X, train$y)
  }

  fit <- regression$fit

  f <- list("unsm"=fit(test$X))

  difference_fit <- regression$deriv_fit

  df <- list("diff"=difference_fit(test$X, D=tune$difference))

  sm_bw <- sm_bw_out$bw_opt
  names(sm_bw) <- names(sm_bw_out$bw_opt_inds)

  smooth_list <- resmooth(fit=fit,
                          X=test$X,
                          d=1,
                          bw=sm_bw)


  f <- c(f, smooth_list$pred)
  names(f) <- c("unsm", names(sm_bw))
  df <- c(df, smooth_list$deriv)
  names(df) <- c("diff", names(sm_bw))

  ### Score estimation

  # Basis score
  basis_lambda <- tune$basis$lambda_min
  basis_degree <- 2
  score <- list("basis_poly" = basis_poly(X=train$X, d=1,
                                          degree=basis_degree,
                                          lambda=basis_lambda)$rho(test$X))

  pred_reg <- function(X,y){fit_xgboost(X = X,
                                        y = y,
                                        params = tune$xgb_m)}
  m_fit <- pred_reg(X=train$X[,-1], y=train$X[,1])$fit

  # Variance estimation
  resid_train <- train$X[,1] - m_fit(train$X[,-1])
  tree <- partykit::ctree(resid_train^2 ~ .,
                          data=data.frame("resid_train"=resid_train,
                                          "z"=train$X[,-1]))
  sigma_train <- sqrt(unname(stats::predict(tree,
                                            newdata = data.frame("z"=train$X[,-1]))))
  eps_train <- resid_train / sigma_train


  # Spline score
  spl_df <- c(tune$spl$df_min, tune$spl$df_1se)
  rho <- spline_score(eps_train, df=spl_df, tol=1e-3)$rho

  # Evaluate
  resid_test <- test$X[,1] - m_fit(test$X[,-1])
  sigma_test <- sqrt(unname(stats::predict(tree,
                                           newdata = data.frame("z"=test$X[,-1]))))

  # If any variance estimates are close to zero, assume homogeneous (more stable)
  if (any(sigma_test < 0.01)){
    warning(paste0("Minimal heterogeneous variance equals ",
                   round(min(sigma_test),4), ", reducing to homogeneous
                   case for stability."))
    sigma_test <- rep(sqrt(mean(resid_test^2)), length(test$y))
  }

  eps_test <- resid_test / sigma_test
  rho_test <- rho(eps_test)

  # If any score estimates are large, assume Gaussian (more stable)
  for (i in seq(length(rho_test))){
    if (any(abs(rho_test[[i]]) > 10)){
      warning(paste0("Maximal absolute spline score prediction equals ",
                     round(max(abs(rho_test[[i]])),2), ", reducing to Gaussian
                   case for stability."))
      rho_test[[i]] <- - eps_test
    }
  }

  score_test <- lapply(rho_test, function(x){x/sigma_test}) # rho_P(x,z)=1/sigma(z) * rho_\varepsilon( (x-m(z))/sigma(z) )


  score <- c(score, score_test)

  names(score) <- c("basis", "min", "1se")

  out <- list("f" = f,
              "df" = df,
              "score" = score)

  return(out)

}


simulate_train_test <- function(n_tr, n_te, ex_setting, f_setting){
  data <- simulate_data(n = n_tr + n_te,
                        ex_setting = ex_setting,
                        f_setting = f_setting)

  smpl <- sample(length(data$x), size=n_tr, replace=FALSE)
  train <- list("x" = data$x[smpl],
                "z" = data$z[smpl,],
                "y" = data$y[smpl],
                "f" = data$f[smpl])
  test <- list("x" = data$x[-smpl],
               "z" = data$z[-smpl,],
               "y" = data$y[-smpl],
               "f" = data$f[-smpl])

  return(list("train" = train,
              "test" = test))

}


compare_partially_linear <- function(n, ex_setting, f_setting) {
  #' Generate simulation data and evaluate partially linear estimator.
  #'
  #' @param n integer number of samples. For "401k" ex_setting this is ignored
  #'     and the whole data set is used.
  #' @param ex_setting string "normal", "mixture2", "mixture3",
  #'     "logistic", "t4", "401k".
  #' @param f_setting string "plm", "additive", "interaction".
  #'
  #' @return list containing estimate, standard error estimate, and sample theta (for debugging).

  data <- simulate_data(n = n,
                        ex_setting = ex_setting,
                        f_setting = f_setting)

  y <- data$y
  X <- cbind(data$x, data$z)


  tune <- rjson::fromJSON(file = paste0("data-raw/tuning/tune_", f_setting, "_", ex_setting, "_results.json"))


  pl <- partially_linear(X=X, y=y,
                         g_params = tune$xgb_f_pl,
                         m_params = tune$xgb_m)



  results <- list("est"=pl$est,
                  "se"=pl$se,
                  "sample_theta"=data$sample_theta)

  return(results)

}

compare_rothenhausler <- function(n, ex_setting, f_setting) {
  #' Generate simulation data and evaluate Rothenhausler estimator.
  #'
  #' @param n integer number of samples. For "401k" ex_setting this is ignored
  #'     and the whole data set is used.
  #' @param ex_setting string "normal", "mixture2", "mixture3",
  #'     "logistic", "t4", "401k".
  #' @param f_setting string "plm", "additive", "interaction".
  #'
  #' @return list containing estimate, standard error estimate, and sample theta (for debugging).

  data <- simulate_data(n = n,
                        ex_setting = ex_setting,
                        f_setting = f_setting)

  y <- data$y
  X <- cbind(data$x, data$z)


  tune <- rjson::fromJSON(file = paste0("data-raw/tuning/tune_", f_setting, "_", ex_setting, "_results.json"))


  roth <- rothenhausler_yu(X=X, y=y,
                           f_lambda = tune$rothenhausler_f$lambda_min,
                           m_lambda = tune$rothenhausler_m$lambda_min)



  results <- list("est"=roth$est,
                  "se"=roth$se,
                  "sample_theta"=data$sample_theta)

  return(results)

}

compare_lm <- function(n, ex_setting, f_setting) {
  #' Generate simulation data and evaluate OLS estimator.
  #'
  #' @param n integer number of samples. For "401k" ex_setting this is ignored
  #'     and the whole data set is used.
  #' @param ex_setting string "normal", "mixture2", "mixture3",
  #'     "logistic", "t4", "401k".
  #' @param f_setting string "plm", "additive", "interaction".
  #'
  #' @return list containing estimate, standard error estimate, and sample theta (for debugging).

  data <- simulate_data(n = n,
                        ex_setting = ex_setting,
                        f_setting = f_setting)

  y <- data$y
  X <- cbind(data$x, data$z)

  lm1 <- stats::lm(y ~ X)



  results <- list("est" = unname(lm1$coefficients[2]),
                  "se"= stats::coef(summary(lm1))[2, "Std. Error"],
                  "sample_theta" = data$sample_theta)

  return(results)

}
