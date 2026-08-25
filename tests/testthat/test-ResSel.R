# Regression test for evobir-07 (low tail selected percent-1
# individuals due to asymmetric > vs < comparisons)
test_that("ResSel selects the same number in both tails", {
  set.seed(9)
  df <- data.frame(
    id = paste0("ind_", 1:100),
    trait1 = rnorm(100, 50, 10),
    trait2 = rnorm(100, 30, 5)
  )
  pdf(NULL)
  on.exit(dev.off())
  res <- ResSel(df, traits = c(2, 3), percent = 10)
  expect_equal(length(res$`high line`), 10L)
  expect_equal(length(res$`low line`), 10L)
})

# evobir-07: unsupported model must error informatively
test_that("ResSel rejects unsupported model", {
  df <- data.frame(id = 1:20, a = rnorm(20), b = rnorm(20))
  expect_error(ResSel(df, traits = c(2, 3), model = "quadratic"))
})

test_that("ResSel returns a list with high and low lines", {
  data <- read.csv(file = system.file("horn.beetle.csv", package = "evobiR"))
  result <- ResSel(data = data, traits = c(2, 3), percent = 15,
                   identifier = 1, model = "linear")
  expect_true(is.list(result))
  expect_true("high line" %in% names(result))
  expect_true("low line" %in% names(result))
})
