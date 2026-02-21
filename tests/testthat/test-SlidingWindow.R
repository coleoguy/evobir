test_that("SlidingWindow works on vectors", {
  x <- 1:100
  result <- SlidingWindow(FUN = "mean", data = x, window = 10, step = 5, strict = TRUE)
  expect_true(is.numeric(result))
  expect_true(length(result) > 0)
  # first window should be mean of 1:10 = 5.5
  expect_equal(result[1], 5.5)
})

test_that("SlidingWindow works on matrices", {
  x <- matrix(1:100, 10, 10)
  result <- SlidingWindow(FUN = "mean", data = x, window = 5, step = 3, strict = TRUE)
  expect_true(is.matrix(result))
})

test_that("SlidingWindow handles different functions", {
  x <- 1:100
  result_sd <- SlidingWindow(FUN = "sd", data = x, window = 10, step = 10, strict = TRUE)
  expect_true(is.numeric(result_sd))
  expect_true(all(result_sd > 0))
})
