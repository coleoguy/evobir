# Regression tests for evobir-01 (state labels used as numeric indices)
# and evobir-05 (multiPhylo branch dropped rejmax/rejint/monitor)

# helper: a simmap is valid if per-edge map times sum to the edge
# lengths and every map segment is named by a real state label
expect_valid_simmap <- function(res, state.labels) {
  expect_true(inherits(res, "simmap"))
  map.sums <- vapply(res$maps, sum, numeric(1))
  expect_equal(unname(map.sums), unname(res$edge.length), tolerance = 1e-8)
  expect_true(all(names(unlist(res$maps)) %in% state.labels))
  expect_equal(unname(rowSums(res$mapped.edge)), unname(res$edge.length),
               tolerance = 1e-8)
}

test_that("make.simmap2 works with a 0/1-coded binary trait", {
  skip_if_not_installed("phytools")
  set.seed(42)
  tree <- ape::rtree(12)
  x <- setNames(sample(0:1, 12, replace = TRUE), tree$tip.label)
  capture.output(
    res <- make.simmap2(tree, x, model = "ER", nsim = 1, message = FALSE)
  )
  expect_valid_simmap(res, c("0", "1"))
})

test_that("make.simmap2 works with a 1/2-coded binary trait", {
  skip_if_not_installed("phytools")
  set.seed(42)
  tree <- ape::rtree(12)
  x <- setNames(sample(1:2, 12, replace = TRUE), tree$tip.label)
  capture.output(
    res <- make.simmap2(tree, x, model = "ER", nsim = 1, message = FALSE)
  )
  expect_valid_simmap(res, c("1", "2"))
})

test_that("make.simmap2 works with chromosome-number-like state labels", {
  skip_if_not_installed("phytools")
  set.seed(7)
  tree <- ape::rtree(15)
  x <- setNames(sample(c(5, 12), 15, replace = TRUE), tree$tip.label)
  capture.output(
    res <- make.simmap2(tree, x, model = "ER", nsim = 1, message = FALSE)
  )
  expect_valid_simmap(res, c("5", "12"))
})

test_that("make.simmap2 honors rejmax for multiPhylo input", {
  skip_if_not_installed("phytools")
  set.seed(1)
  trees <- c(ape::rtree(8), ape::rtree(8))
  class(trees) <- "multiPhylo"
  x <- setNames(sample(0:1, 8, replace = TRUE), trees[[1]]$tip.label)
  # rejmax = 0 forces every branch to fail immediately; the old code
  # silently dropped rejmax for multiPhylo input (called plain
  # make.simmap), so no "fail" segments were produced
  capture.output(
    res <- make.simmap2(trees, x, model = "ER", nsim = 1,
                        rejmax = 0, message = FALSE)
  )
  expect_true(inherits(res, "multiSimmap"))
  expect_equal(length(res), 2L)
  has.fail <- vapply(res, function(tr) {
    "fail" %in% names(unlist(tr$maps))
  }, logical(1))
  expect_true(all(has.fail))
})
