test_that("Prop_SA returns a PropSA object that sums to 1", {
  result <- Prop_SA(Da = 26, scs = "XY", mud = 0.4)
  expect_s3_class(result, "PropSA")
  expect_equal(result$P_SA + result$P_AA + result$P_SS, 1, tolerance = 1e-10)
})

test_that("Prop_SA returns numeric values between 0 and 1", {
  result <- Prop_SA(Da = 26, scs = "XY", mud = 0.4)
  expect_true(is.numeric(result$P_SA))
  expect_true(result$P_SA >= 0 && result$P_SA <= 1)
})

test_that("Prop_SA works with XO system", {
  result <- Prop_SA(Da = 20, scs = "XO", mud = 0.5)
  expect_s3_class(result, "PropSA")
  expect_true(result$P_SA >= 0 && result$P_SA <= 1)
})

test_that("Prop_SA works with ZW system", {
  result <- Prop_SA(Da = 20, scs = "ZW", mud = 0.5)
  expect_s3_class(result, "PropSA")
  expect_true(result$P_SA >= 0 && result$P_SA <= 1)
})
