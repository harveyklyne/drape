xgb_tune_grid <- function(n_tr, n_te, f_setting, ex_setting, partially_linear = FALSE){

  if (f_setting=="m"){
    data <- simulate_train_test(n_tr = n_tr,
                                n_te = n_te,
                                ex_setting = ex_setting,
                                f_setting = "plm")
    train_data <- data$train
    test_data <- data$test

    train <- xgboost::xgb.DMatrix(data = train_data$z,
                                  label = train_data$x)
    test <- xgboost::xgb.DMatrix(data = test_data$z,
                                 label = test_data$x)
  } else{
    data <- simulate_train_test(n_tr = n_tr,
                                n_te = n_te,
                                ex_setting = ex_setting,
                                f_setting = f_setting)
    train_data <- data$train
    test_data <- data$test

    if (partially_linear){ # partial out x term using population theta

      if (f_setting == "plm"){ theta <- 1 }
      else {theta <- readRDS(paste0("data-raw/variables/theta/theta_", ex_setting, "_", f_setting, ".rds"))}

      train <- xgboost::xgb.DMatrix(data = train_data$z,
                                    label = train_data$y - theta*train_data$x)
      test <- xgboost::xgb.DMatrix(data = test_data$z,
                                   label = test_data$f - theta*test_data$x)

    } else{
      train <- xgboost::xgb.DMatrix(data = cbind(train_data$x,train_data$z),
                                    label = train_data$y)
      test <- xgboost::xgb.DMatrix(data = cbind(test_data$x,test_data$z),
                                   label = test_data$f)
    }

  }


  watchlist <- list(train=train, test=test)

  depths <- 1:8
  gammas <- seq(0,4, by=0.2)
  maxrounds <- 5000
  test_rmse <- array(NA, c(length(depths), length(gammas), maxrounds))

  for (i in seq(length(depths))){
    for (j in seq(length(gammas))){
      f_fit <- xgboost::xgb.train(data=train,
                                  max.depth = depths[i],
                                  eta = 0.01,
                                  alpha=0,
                                  gamma = gammas[j],
                                  nthread = 2,
                                  nrounds = maxrounds,
                                  watchlist=watchlist,
                                  objective = "reg:squarederror",
                                  verbose = 0)
      test_rmse[i,j,] <- f_fit$evaluation_log$test_rmse
    }
  }

  return(test_rmse)
}


spline_tune <- function(n_tr, n_te, ex_setting, xgb_params, df){

  data <- simulate_train_test(n_tr = n_tr,
                              n_te = n_te,
                              ex_setting = ex_setting,
                              f_setting = "plm")
  train_data <- data$train
  test_data <- data$test

  train <- spline_tune_generate(data = train_data,
                                xgb_params = xgb_params)
  test <- spline_tune_generate(data = test_data,
                               xgb_params = xgb_params)

  spline <- spline_score(x = train, df = df)

  score <- spline$rho(test)
  dscore <- spline$drho(test)

  # If any score estimates are large, assume Gaussian (more stable)
  for (i in seq(length(score))){
    if (any(abs(score[[i]]) > 10)){
      warning(paste0("Maximal absolute spline score prediction equals ",
                     round(max(abs(score[[i]])),2), ", reducing to Gaussian
                   case for stability."))
      score[[i]] <- - test
      dscore[[i]] <- rep(-1 , length(test))
    }
  }



  eval <- mapply(function(x,y){x^2+2*y},
                 x = score,
                 y = dscore)

  cv <- colMeans(eval)
  se <- apply(eval,2,stats::sd)/sqrt(n_te)

  dfmin_index <- which.min(cv)
  cv_target <- cv[dfmin_index]+se[dfmin_index]
  df1se_index <- which(df==min(df[cv < cv_target]))

  df_min <- df[dfmin_index]
  df_1se <- df[df1se_index]

  out <- c(df_min, df_1se, cv)
  names(out) <- c("df_min", "df_1se", paste0("cv", df))

  return(out)

}

spline_tune_generate  <- function(data, xgb_params){

  fit <- fit_xgboost(X = data$z,
                     y = data$x,
                     params = xgb_params)$fit

  pred <- fit(data$z)
  resid <- data$x - pred

  tree <- partykit::ctree(resid^2 ~ ., data=data.frame("resid"=resid, "z"=data$z))
  sigma <- sqrt(unname(stats::predict(tree, newdata = data.frame("z"=data$z))))

  # If any variance estimates are close to zero, assume homogeneous (more stable)
  if (any(sigma < 0.01)){
    warning(paste0("Minimal heterogeneous variance equals ",
                   round(min(sigma),4), ", reducing to homogeneous
                   case for stability."))
    sigma <- rep(sqrt(mean(resid^2)), length(resid))
  }

  eps <- resid / sigma

  return(eps)
}


basis_tune <- function(n_tr, n_te, ex_setting, lambda){

  data <- simulate_train_test(n_tr = n_tr,
                              n_te = n_te,
                              ex_setting = ex_setting,
                              f_setting = "plm")
  train_data <- data$train
  test_data <- data$test

  test_output <- basis_tune_generate(train_data = train_data,
                                     test_data = test_data,
                                     degree = 2,
                                     lambda = lambda)

  eval <- test_output$test_rho^2 + 2*test_output$test_drho

  cv <- colMeans(eval)
  se <- apply(eval,2,stats::sd)/sqrt(n_te)

  lambdamin_index <- which.min(cv)
  cv_target <- cv[lambdamin_index]+se[lambdamin_index]
  lambda1se_index <- which(lambda==max(lambda[cv < cv_target]))

  lambda_min <- lambda[lambdamin_index]
  lambda_1se <- lambda[lambda1se_index]

  out <- c(lambda_min, lambda_1se, cv)
  names(out) <- c("lambda_min", "lambda_1se", paste0("cv", 1:length(lambda)))

  return(out)

}

basis_tune_generate <- function(train_data, test_data, degree, lambda){

  if (degree == 2){
    train_basis <- as.matrix(stats::poly(cbind(train_data$x, train_data$z),
                                         degree = degree,
                                         raw=TRUE))[,-1]

    test_basis <- as.matrix(stats::poly(cbind(test_data$x, test_data$z),
                                        degree = degree,
                                        raw=TRUE))[,-1]
    test_dbasis <- matrix(0, nrow(test_basis), ncol(test_basis))
    test_dbasis[,c(1, (2:(1+ncol(test_data$z)))*(1 + 2:(1+ncol(test_data$z)))/2)] <- cbind(2*test_data$x, test_data$z)
  }

  fitted_model <- glmnet::glmnet(train_basis, train_data$x, lambda = lambda)
  train_pred <- stats::predict(fitted_model, train_basis)
  train_var <- mean(train_data$x^2) - colMeans(train_data$x * train_pred)

  test_pred <- stats::predict(fitted_model, test_basis)
  test_dpred <- test_dbasis %*% as.matrix(fitted_model$beta)

  test_rho <- - (test_data$x - test_pred) / rep(train_var, each=length(test_data$x))
  test_drho <- - (1 - test_dpred) / rep(train_var, each=length(test_data$x))

  return(list("test_rho" = unname(test_rho),
              "test_drho" = unname(test_drho)))

}


getStddev <- function(f_setting, ex_setting, n){

  data <- simulate_data(n = n,
                        ex_setting = ex_setting,
                        f_setting = "plm")

  stddev <- stats::sd(data$x)

  return(stddev)
}


rothenhausler_tune <- function(n_tr, n_te, ex_setting, f_setting, lambda){

  if (f_setting=="m"){
    data <- simulate_train_test(n_tr = n_tr,
                                n_te = n_te,
                                ex_setting = ex_setting,
                                f_setting = "plm")

    train <- list("X" = data$train$z,
                  "y" = data$train$x)
    test <- list("X" = data$test$z,
                  "y" = data$test$x)

  } else{
    data <- simulate_train_test(n_tr = n_tr,
                                n_te = n_te,
                                ex_setting = ex_setting,
                                f_setting = f_setting)

    train <- list("X" = rothenhausler_basis(cbind(data$train$x, data$train$z))$mod_basis,
                  "y" = data$train$y)
    test <- list("X" = rothenhausler_basis(cbind(data$test$x, data$test$z))$mod_basis,
                 "y" = data$test$y)
  }

  fitted_model <- glmnet::glmnet(x = train$X, y = train$y, lambda = lambda)
  test_pred <- stats::predict(fitted_model, test$X)

  eval <- (test$y - test_pred)^2

  cv <- colMeans(eval)
  se <- apply(eval,2,stats::sd)/sqrt(n_te)

  lambdamin_index <- which.min(cv)
  cv_target <- cv[lambdamin_index]+se[lambdamin_index]
  lambda1se_index <- which(lambda==max(lambda[cv < cv_target]))

  lambda_min <- lambda[lambdamin_index]
  lambda_1se <- lambda[lambda1se_index]

  out <- c(lambda_min, lambda_1se, cv)
  names(out) <- c("lambda_min", "lambda_1se", paste0("cv", 1:length(lambda)))

  return(out)

}
