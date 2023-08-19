test_that("simulate_data works in simulated settings", {
  n <- 5
  for (ex in c("normal", "mixture2", "mixture3", "logistic", "t4")){
    for (f in c("plm", "additive", "interaction")){
      out <- simulate_data(n, ex, f)
      expect_vector(out$y, size=n)
      expect_vector(out$x, size=n)
      expect_equal(nrow(out$z), n)
    }
  }
})
