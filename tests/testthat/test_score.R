test_that("spline_score acts correctly for single df", {
  x <- seq(10)
  spl <- spline_score(x, df=6)
  expect_vector(spl$rho(x), size=10)
  expect_vector(spl$drho(x), size=10)
})

test_that("spline_score acts correctly for multiple df", {
  x <- seq(10)
  spl <- spline_score(x, df=c(2,3,4))
  out <- spl$rho(seq(5))
  expect_equal(sapply(out, length), rep(5,3))
})

test_that("ng_pseudo_y gives correct form of output", {
  expect_vector(ng_pseudo_response(1:5), size=5)
})

test_that("basis_poly acts correctly for single df", {
  X <- matrix(stats::rnorm(20), ncol=4)
  bs <- basis_poly(X=X, d=1, degree=2, lambda=1)
  expect_vector(bs$rho(X), size=5)
})

