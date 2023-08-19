test_that("fit_lasso_poly acts correctly", {
  X <- matrix(seq(15), ncol=3)
  y <- seq(5)
  out <- fit_lasso_poly(X=X, y=y, degree=2, lambda=1)
  expect_vector(out$fit(X), 5)
})

test_that("fit_xgboost acts correctly", {
  X <- matrix(stats::rnorm(30), ncol=3)
  y <- X[,1] + 0.1 * stats::rnorm(nrow(X))
  params <- list("eta"=1,
                 "max.depth"=2,
                 "gamma"=1,
                 "nrounds"=5)
  out <- fit_xgboost(X=X, y=y, params=params, derivative=TRUE)
  expect_vector(out$fit(X), 10)
  expect_vector(out$deriv_fit(X, D=1), 10)
})
