test_that("ResSel returns a list with high and low lines", {
  data <- read.csv(file = system.file("horn.beetle.csv", package = "evobiR"))
  result <- ResSel(data = data, traits = c(2, 3), percent = 15,
                   identifier = 1, model = "linear")
  expect_true(is.list(result))
  expect_true("high line" %in% names(result))
  expect_true("low line" %in% names(result))
})
