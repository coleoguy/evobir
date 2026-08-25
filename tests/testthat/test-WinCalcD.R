# Regression tests for evobir-08 (WinCalcD returned all-character
# columns and dropped the last window)

test_that("WinCalcD returns numeric columns and keeps the last window", {
  aln <- system.file("1.fasta", package = "evobiR")
  capture.output(
    res <- WinCalcD(aln, win.size = 100, step.size = 50)
  )
  expect_true(is.data.frame(res))
  expect_true(is.numeric(res$abba))
  expect_true(is.numeric(res$baba))
  expect_true(is.numeric(res$d))
  expect_true(is.character(res$range))

  # window count must include a final window ending at the last site
  # when (total - win.size) is divisible by step.size
  aln.len <- nchar(gsub("\\s", "",
    paste(readLines(aln)[!startsWith(readLines(aln), ">")][1], collapse = "")))
  # ranges parse back to integers; final start must be > total - win.size - step.size + 1
  starts <- as.integer(sub(":.*", "", res$range))
  expected.starts <- seq(1, aln.len - 100 + 1, by = 50)
  expect_equal(starts, expected.starts)
})

test_that("WinCalcD D values match CalcD on a single full-length window", {
  aln <- system.file("1.fasta", package = "evobiR")
  ln <- readLines(aln)
  aln.len <- nchar(ln[which(!startsWith(ln, ">"))[1]])
  capture.output({
    res <- WinCalcD(aln, win.size = aln.len, step.size = aln.len)
    d.full <- CalcD(aln, sig.test = "N")
  })
  expect_equal(nrow(res), 1L)
  expect_equal(res$d[1], d.full, tolerance = 1e-12)
})
