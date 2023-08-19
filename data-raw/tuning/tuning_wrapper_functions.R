tune_xgb <- function(f_setting, ex_setting, reps, n_tr, n_te, partially_linear = FALSE){

  if (partially_linear){print(paste0("Tuning XGB partially linear ", f_setting, " ", ex_setting))}
  else {print(paste0("Tuning XGB ", f_setting, " ", ex_setting))}

  sim_df <- expand.grid(reps=1:reps)

  with_progress( {
    prog_bar <- progressor(along=1:(nrow(sim_df)))
    sim_res <- future_apply(sim_df, MARGIN=1, future.seed=TRUE, future.scheduling = Inf, simplify=FALSE, FUN = function(x) {
      prog_bar()

      out <- suppressWarnings(xgb_tune_grid(n_tr = n_tr,
                                            n_te = n_te,
                                            f_setting = f_setting,
                                            ex_setting = ex_setting,
                                            partially_linear = partially_linear))
      return(out)
    })
  })

  mean_mse <- Reduce("+", sim_res) / length(sim_res)
  mean_mse_out <- melt(mean_mse)
  mean_mse_out[,2] <- seq(0,4, by=0.2)[mean_mse_out[,2]]
  names(mean_mse_out) = c("depth", "gamma", "nround", "mse")

  write.csv(mean_mse_out, paste0("xgb_tune_", f_setting, "_", ex_setting, "_results.csv"), row.names=FALSE)

  tmp <- arrayInd(which.min(mean_mse), dim(mean_mse))
  depth <- seq(8)[tmp[1]]
  gamma <- seq(0,4, by=0.2)[tmp[2]]
  nrounds <- seq(5000)[tmp[3]]

  result <- list("eta" = 0.01,
                 "max.depth" = depth,
                 "gamma" = gamma,
                 "nrounds" = nrounds)

  return(result)

}


tune_basis <- function(ex_setting, reps, n_tr, n_te){

  print(paste0("Tuning basis ", ex_setting))

  # Get lambda path

  lambda <- getLambdas(n = n_tr + n_te,
                       ex_setting = ex_setting)


  # Evaluate cv scores

  sim_df <- expand.grid(reps=1:reps)

  with_progress( {
    prog_bar <- progressor(along=1:(nrow(sim_df)))
    sim_res <- future_apply(sim_df, MARGIN=1, future.seed=TRUE, future.scheduling = Inf, simplify=FALSE, FUN = function(x) {
      prog_bar()

      out <- suppressWarnings(basis_tune(n_tr = n_tr,
                                         n_te = n_te,
                                         ex_setting = ex_setting,
                                         lambda = lambda))
      return(out)
    })
  })


  sim_res_df <- cbind(sim_df, data.frame(Reduce(rbind, sim_res), row.names = NULL))
  write.csv(sim_res_df, paste0("basis_tune_", ex_setting, ".csv"), row.names=FALSE)

  cv <- colMeans(sim_res_df[,4:ncol(sim_res_df)])
  se <- apply(sim_res_df[,4:ncol(sim_res_df)],2,stats::sd)/sqrt(reps)

  lambdamin_index <- which.min(cv)
  cv_target <- cv[lambdamin_index]+se[lambdamin_index]
  lambda1se_index <- which(lambda==max(lambda[cv < cv_target]))

  lambda_min <- lambda[lambdamin_index]
  lambda_1se <- lambda[lambda1se_index]

  result <- list("lambda_min" = lambda_min,
                 "lambda_1se" = lambda_1se)

  return(result)

}

getLambdas <- function(n, ex_setting){
  # from https://stackoverflow.com/questions/23686067/default-lambda-sequence-in-glmnet-for-cross-validation

  data <- simulate_data(n = n,
                        f_setting = "plm",
                        ex_setting = ex_setting)

  X <- cbind(data$x, data$z)
  basis <- as.matrix(stats::poly(X, degree = 2, raw=TRUE))
  basis <- basis[,-1]

  n <- length(data$x)

  ## Standardize variables: (need to use n instead of (n-1) as denominator)
  mysd <- function(z) sqrt(sum((z-mean(z))^2)/n)
  sx <- as.matrix(scale(basis, scale = apply(basis, 2, mysd)))

  ## Calculate lambda path (first get lambda_max):
  lambda_max <- max(abs(colSums(sx*data$x)))/n
  epsilon <- .0001
  K <- 100
  lambdapath <- round(exp(seq(log(lambda_max),
                              log(lambda_max*epsilon),
                              length.out = K)), digits = 10)
  return(lambdapath)

}


tune_spline <- function(ex_setting,
                        xgb_params,
                        reps, n_tr, n_te){

  print(paste0("Tuning spline ", ex_setting))

  # Get df path

  df <- seq(2, 20, 0.5)


  # Evaluate cv scores

  sim_df <- expand.grid(reps=1:reps)

  with_progress( {
    prog_bar <- progressor(along=1:(nrow(sim_df)))
    sim_res <- future_apply(sim_df, MARGIN=1, future.seed=TRUE, future.scheduling = Inf, simplify=FALSE, FUN = function(x) {
      prog_bar()

      out <- suppressWarnings(spline_tune(n_tr = n_tr,
                                          n_te = n_te,
                                          ex_setting = ex_setting,
                                          xgb_params = xgb_params,
                                          df = df))
      return(out)
    })
  })

  sim_res_df <- cbind(sim_df, data.frame(Reduce(rbind, sim_res), row.names = NULL))
  write.csv(sim_res_df, paste0("spline_tune_", ex_setting, ".csv"), row.names=FALSE)

  cv <- colMeans(sim_res_df[,4:ncol(sim_res_df)])
  se <- apply(sim_res_df[,4:ncol(sim_res_df)],2,stats::sd)/sqrt(reps)

  dfmin_index <- which.min(cv)
  cv_target <- cv[dfmin_index]+se[dfmin_index]
  df1se_index <- which(df==min(df[cv < cv_target]))

  df_min <- df[dfmin_index]
  df_1se <- df[df1se_index]

  result <- list("df_min" = df_min,
                 "df_1se" = df_1se)

  return(result)

}

tune_stddev <- function(ex_setting, reps, n){

  print(paste0("Tuning stddev ", ex_setting))

  sim_df <- expand.grid(reps=1:reps)

  with_progress( {
    prog_bar <- progressor(along=1:(nrow(sim_df)))
    sim_res <- future_apply(sim_df, MARGIN=1, future.seed=TRUE, future.scheduling = Inf, simplify=FALSE, FUN = function(x) {
      prog_bar()

      out <- suppressWarnings(getStddev(ex_setting = ex_setting,
                                        n = n))
      return(out)
    })
  })

  devs <- Reduce(c, sim_res)

  stddev <- sqrt(mean(devs^2))

  return(stddev)

}



tune_rothenhausler <- function(ex_setting, f_setting, reps, n_tr, n_te){

  print(paste0("Tuning Rothenhausler ", f_setting, " ", ex_setting))

  # Get lambda path

  lambda <- getRothenhauslerLambdas(n = n_tr + n_te,
                       ex_setting = ex_setting,
                       f_setting = f_setting)


  # Evaluate cv scores

  sim_df <- expand.grid(reps=1:reps)

  with_progress( {
    prog_bar <- progressor(along=1:(nrow(sim_df)))
    sim_res <- future_apply(sim_df, MARGIN=1, future.seed=TRUE, future.scheduling = Inf, simplify=FALSE, FUN = function(x) {
      prog_bar()

      out <- suppressWarnings(rothenhausler_tune(n_tr = n_tr,
                                         n_te = n_te,
                                         ex_setting = ex_setting,
                                         f_setting = f_setting,
                                         lambda = lambda))
      return(out)
    })
  })


  sim_res_df <- cbind(sim_df, data.frame(Reduce(rbind, sim_res), row.names = NULL))
  write.csv(sim_res_df, paste0("rothenhausler_tune_", f_setting, "_", ex_setting, ".csv"), row.names=FALSE)

  cv <- colMeans(sim_res_df[,4:ncol(sim_res_df)])
  se <- apply(sim_res_df[,4:ncol(sim_res_df)],2,stats::sd)/sqrt(reps)

  lambdamin_index <- which.min(cv)
  cv_target <- cv[lambdamin_index]+se[lambdamin_index]
  lambda1se_index <- which(lambda==max(lambda[cv < cv_target]))

  lambda_min <- lambda[lambdamin_index]
  lambda_1se <- lambda[lambda1se_index]

  result <- list("lambda_min" = lambda_min,
                 "lambda_1se" = lambda_1se)

  return(result)

}


getRothenhauslerLambdas <- function(n, ex_setting, f_setting){
  # from https://stackoverflow.com/questions/23686067/default-lambda-sequence-in-glmnet-for-cross-validation


  data <- simulate_data(n = n,
                        f_setting = f_setting,
                        ex_setting = ex_setting)
  X <- rothenhausler_basis(cbind(data$x, data$z))

  if (f_setting == "m"){
    y <- data$x
    X <- X[,-1]
  } else{ y <- data$y }

  n <- length(data$x)

  ## Standardize variables: (need to use n instead of (n-1) as denominator)
  mysd <- function(z) sqrt(sum((z-mean(z))^2)/n)
  sx <- as.matrix(scale(X, scale = apply(X, 2, mysd)))

  ## Calculate lambda path (first get lambda_max):
  lambda_max <- max(abs(colSums(sx*y)))/n
  epsilon <- .0001
  K <- 100
  lambdapath <- round(exp(seq(log(lambda_max),
                              log(lambda_max*epsilon),
                              length.out = K)), digits = 10)
  return(lambdapath)

}
