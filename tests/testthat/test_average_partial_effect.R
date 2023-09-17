test_that("drape works with simple nuisance estimators", {
  set.seed(0)
  data <- simulate_data(200, "normal", "plm")
  response_regression <- function(X,y){
    df <- data.frame(y,X)
    colnames(df) <- c("y", paste0("X", 1:10))
    lm1 <- stats::lm(y~X1+sin(X2), data=df)
    fit <- function(newX){
      newdf <- data.frame(newX)
      colnames(newdf) <- paste0("X", 1:10)
      return(as.vector(stats::predict(lm1, newdata=newdf)))}
    return(list("fit"=fit))
  }
  predictor_regression <- function(z,x){
    df <- data.frame(x,z)
    colnames(df) <- c("x", paste0("Z", 1:9))
    lm1 <- stats::lm(x~Z1+Z2, data=df)
    fit <- function(newz){
      newdf <- data.frame(newz)
      colnames(newdf) <- paste0("Z", 1:9)
      return(as.vector(stats::predict(lm1, newdata=newdf)))}
    return(list("fit"=fit))
  }
  out <- drape(data$y, data$x, data$z, response_regression, predictor_regression, nfolds=2)
  expect_vector(out$est, 1)
  expect_vector(out$se, 1)
})

test_that("drape works with xgboost nuisance estimators", {
  set.seed(0)
  data <- simulate_data(1000, "normal", "plm")
  params <- list("eta" = 0.1, "max_depth" = 2, "nrounds" = 100)
  response_regression <- predictor_regression <- function(X,y){
      fit_xgboost(X = X, y = y, params = params)}
  out <- drape(data$y, data$x, data$z, response_regression, predictor_regression)
  expect_vector(out$est, 1)
  expect_vector(out$se, 1)
})
