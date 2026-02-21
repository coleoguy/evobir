test_that("CalcD returns a numeric D-statistic", {
  result <- CalcD(alignment = system.file("1.fasta", package = "evobiR"),
                  sig.test = "N")
  expect_true(is.numeric(result))
  expect_true(result >= -1 && result <= 1)
})

test_that("CalcPopD returns a list with required elements", {
  result <- CalcPopD(alignment = system.file("3.fasta", package = "evobiR"),
                     sig.test = "N")
  expect_true(is.list(result))
  expect_true("d.stat" %in% names(result))
})

test_that("CalcD with bootstrap returns numeric", {
  result <- CalcD(alignment = system.file("1.fasta", package = "evobiR"),
                  sig.test = "B", replicate = 10)
  expect_true(is.numeric(result))
})
