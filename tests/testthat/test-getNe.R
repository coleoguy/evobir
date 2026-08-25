test_that("getNe returns numeric value", {
  result <- getNe(males = 100, females = 200)
  expect_true(is.numeric(result))
  expect_equal(length(result), 1)
})

test_that("getNe with equal sex ratio returns expected value", {
  # With equal numbers, Ne = N
  result <- getNe(males = 100, females = 100)
  expect_equal(result, 200)
})

test_that("getNe works for X-linked locus", {
  result <- getNe(males = 100, females = 200, locus = "X")
  expect_true(is.numeric(result))
  expect_equal(length(result), 1)
})

# Regression test for evobir-18 (invalid locus gave "object 'ne' not
# found" instead of an informative error)
test_that("getNe rejects an unrecognized locus", {
  expect_error(getNe(100, 100, locus = "a"), "arg")
})
