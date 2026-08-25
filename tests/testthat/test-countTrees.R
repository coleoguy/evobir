# Regression tests for evobir-02 (countTrees breaks on single-tree files)

test_that("countTrees works with a single-tree reference file", {
  set.seed(3)
  trees <- c(ape::rtree(5), ape::rtree(5), ape::rtree(5))
  class(trees) <- "multiPhylo"
  coll <- tempfile(fileext = ".nwk")
  ref <- tempfile(fileext = ".nwk")
  ape::write.tree(trees, file = coll)
  # single reference topology: the first tree
  ape::write.tree(trees[[1]], file = ref)
  res <- countTrees(coll, ref, classes = TRUE, verbose = FALSE)
  expect_equal(length(res[[1]]), 1L)
  expect_equal(length(res[[2]]), 3L)
  expect_equal(res[[2]][1], 1)
  unlink(c(coll, ref))
})

test_that("countTrees works with a single-tree collection file", {
  set.seed(4)
  tr <- ape::rtree(5)
  coll <- tempfile(fileext = ".nwk")
  ref <- tempfile(fileext = ".nwk")
  ape::write.tree(tr, file = coll)
  refs <- c(tr, ape::rtree(5))
  class(refs) <- "multiPhylo"
  ape::write.tree(refs, file = ref)
  res <- countTrees(coll, ref, classes = TRUE, verbose = FALSE)
  # exactly one tree in the collection, classified against topology 1
  expect_equal(length(res[[2]]), 1L)
  expect_equal(res[[1]][1], 1)
  unlink(c(coll, ref))
})
