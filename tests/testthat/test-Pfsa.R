test_that("pSAF returns a numeric value between 0 and 1", {
  result <- pSAF(n_auto = 13, n_x = 1, n_y = 1, mu = 0.4)
  expect_true(is.numeric(result))
  expect_length(result, 1)
  expect_true(result >= 0 && result <= 1)
})

test_that("pSAF occurrence model returns value between 0 and 1", {
  result <- pSAF(n_auto = 13, n_x = 1, n_y = 1, mu = 0.4, weighted = FALSE)
  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})

test_that("pSAF works with XO system", {
  result <- pSAF(n_auto = 10, n_x = 1, n_y = 0)
  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})

test_that("pSAF works with ZW system", {
  result <- pSAF(n_auto = 10, n_z = 1, n_w = 1)
  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})

test_that("pSAF is vectorized over n_auto", {
  result <- pSAF(n_auto = 1:10, n_x = 1, n_y = 1)
  expect_length(result, 10)
  expect_true(all(result >= 0 & result <= 1))
})

test_that("pSAF fixation-weighted >= occurrence for standard XY", {
  p_fix <- pSAF(n_auto = 10, n_x = 1, n_y = 1, weighted = TRUE)
  p_occ <- pSAF(n_auto = 10, n_x = 1, n_y = 1, weighted = FALSE)
  expect_true(p_fix >= p_occ)
})
