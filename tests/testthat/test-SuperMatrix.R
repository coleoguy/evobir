# Regression tests for evobir-11 (SuperMatrix read every file in the
# directory, including its own previous outputs on a second run)

make_fasta <- function(path, taxa, len) {
  seqs <- vapply(taxa, function(t)
    paste(sample(c("a", "c", "g", "t"), len, replace = TRUE),
          collapse = ""), character(1))
  writeLines(paste0(">", taxa, "\n", seqs), path)
}

test_that("SuperMatrix ignores its own previous outputs on re-run", {
  wd <- tempfile("supermatrix_test")
  dir.create(wd)
  old.wd <- setwd(wd)
  on.exit({setwd(old.wd); unlink(wd, recursive = TRUE)})

  set.seed(21)
  make_fasta("gene1.fas", c("sp1", "sp2", "sp3"), 12)
  make_fasta("gene2.fas", c("sp1", "sp2", "sp4"), 9)

  capture.output(res1 <- SuperMatrix(input = ".fas", save = TRUE))
  ncol1 <- ncol(res1[[2]])
  expect_equal(ncol1, 21L)

  # Second run in the same directory: the concatenated.fasta written by
  # run 1 must NOT be swept into the result
  capture.output(res2 <- SuperMatrix(input = ".fas", save = TRUE))
  expect_equal(ncol(res2[[2]]), ncol1)
  expect_equal(nrow(res2[[2]]), nrow(res1[[2]]))
})

test_that("SuperMatrix treats input as a literal string, not a regex", {
  wd <- tempfile("supermatrix_regex")
  dir.create(wd)
  old.wd <- setwd(wd)
  on.exit({setwd(old.wd); unlink(wd, recursive = TRUE)})

  set.seed(22)
  make_fasta("geneA.fasta", c("sp1", "sp2"), 6)
  # this file would match the regex ".fasta" interpreted with "." as a
  # wildcard ("xfasta") but not the literal string ".fasta"
  make_fasta("genexfasta", c("sp1", "sp2"), 6)

  capture.output(res <- SuperMatrix(input = ".fasta", save = FALSE))
  expect_equal(ncol(res[[2]]), 6L)
})

test_that("SuperMatrix errors when no files match", {
  wd <- tempfile("supermatrix_empty")
  dir.create(wd)
  old.wd <- setwd(wd)
  on.exit({setwd(old.wd); unlink(wd, recursive = TRUE)})
  expect_error(SuperMatrix(input = ".fasta"), "no alignment files")
})
