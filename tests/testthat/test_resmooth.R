test_that("resmooth acts correctly", {
  X <- matrix(seq(-5,5))
  fit <- function(Y){1*(rowMeans(Y)<0)}
  sm <- resmooth(fit=fit, X=X, d=1, bw=0.2)
  expect_vector(sm$pred[[1]], 11)
})

test_that("new_X acts correctly", {
  expect_equal(new_X(X=matrix(1),
                     d=1,
                     MC_variates=-1:1,
                     bw=1),
               matrix(0:2))
  expect_equal(new_X(X=matrix(1:6, ncol=3),
                     d=1,
                     MC_variates=c(0.1,0.2),
                     bw=1),
               matrix(c(1.1,1.2,2.1,2.2, rep(3:6, each=2)), ncol=3))
  expect_equal(new_X(X=matrix(1:4,ncol=2),
                     d=1,
                     MC_variates=c(0.1,0.2),
                     bw=c(1,0.1)),
               matrix(c(1.1,1.2,2.1,2.2,1.01,1.02,2.01,2.02,3,3,4,4,3,3,4,4),
                      ncol=2))
  expect_equal(new_X(X=matrix(1:4,ncol=2),
                     d=2,
                     MC_variates=c(0.1,0.2),
                     bw=c(1,0.1)),
               matrix(c(1,1,2,2,1,1,2,2,3.1,3.2,4.1,4.2,3.01,3.02,4.01,4.02),
                      ncol=2))
  expect_equal(new_X(X=matrix(1:3,ncol=3),
                     d=2,
                     MC_variates=c(0.1),
                     bw=c(1,0.1)),
               matrix(c(1,1,2.1,2.01,3,3),
                      ncol=3))
  expect_equal(dim(new_X(X=matrix(rep(1,6), ncol=3),
                         d=1,
                         MC_variates=1:3,
                         bw=1:5)), c(30,3))
})

test_that("MC_sums acts correctly", {
  expect_equal(MC_sums(1:10,2,5,1), list(c(15,40)))
  expect_equal(MC_sums(1:30,3,5,2), list(c(15,40,65),c(90,115,140)))
})
