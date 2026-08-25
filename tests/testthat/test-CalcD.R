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

# Regression test for evobir-09 (jackknife dropped block.size + 1 sites
# via an inclusive range and used fractional drop positions)
test_that("CalcD jackknife runs with exact block sizes", {
  capture.output(
    result <- CalcD(alignment = system.file("1.fasta", package = "evobiR"),
                    sig.test = "J", block.size = 50, replicate = 20)
  )
  expect_true(is.numeric(result))
  expect_true(result >= -1 && result <= 1)
})

# Regression test for evobir-13 (vectorized biallelic test must be
# exactly equivalent to the per-column unique-count test)
test_that("vectorized biallelic test matches per-column unique counts", {
  set.seed(31)
  m <- matrix(sample(c("a", "c", "g", "t"), 4 * 500, replace = TRUE),
              nrow = 4)
  expected <- vapply(seq_len(ncol(m)), function(i) {
    length(unique(m[, i])) == 2L
  }, logical(1))
  got <- evobiR:::.biallelic4(m[1, ], m[2, ], m[3, ], m[4, ])
  expect_identical(got, expected)
})
