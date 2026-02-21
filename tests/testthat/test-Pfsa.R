test_that("Pfsa, Pfaa, Pfss sum to 1", {
  pfsa <- Pfsa(Da = 26, scs = "XY", mud = 0.4)
  pfaa <- Pfaa(Da = 26, scs = "XY", mud = 0.4)
  pfss <- Pfss(Da = 26, scs = "XY", mud = 0.4)
  expect_equal(pfsa + pfaa + pfss, 1, tolerance = 1e-10)
})

test_that("Pfsa returns numeric of length 1", {
  result <- Pfsa(Da = 26, scs = "XY", mud = 0.4)
  expect_true(is.numeric(result))
  expect_equal(length(result), 1)
  expect_true(result >= 0 && result <= 1)
})

test_that("Pfsa works with XO system", {
  result <- Pfsa(Da = 20, scs = "XO", mud = 0.5)
  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})

test_that("Pfsa works with ZW system", {
  result <- Pfsa(Da = 20, scs = "ZW", mud = 0.5)
  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})
