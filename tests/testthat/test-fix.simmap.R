# Regression tests for evobir-10 (fix.simmap: root-parent tips
# misaligned start states; stale `branches` leaked across iterations
# via exists())

# helper: mark an edge of a simmap as failed the same way
# make.simmap2 does when rejmax is exceeded
fail_edge <- function(smap, e) {
  smap$maps[[e]] <- setNames(smap$edge.length[e], "fail")
  smap$mapped.edge[e, ] <- 0
  smap
}

test_that("fix.simmap repairs a failed tip branch whose parent is the root", {
  skip_if_not_installed("phytools")
  skip_if_not_installed("igraph")
  set.seed(13)
  # t4 hangs directly off the root: the evobir-10 misalignment case
  tree <- ape::read.tree(text = "(((t1:1,t2:1):1,t3:2):1,t4:3);")
  x <- setNames(c(0, 1, 0, 1), tree$tip.label)
  capture.output(
    res <- make.simmap2(tree, x, model = "ER", nsim = 1, message = FALSE)
  )
  root.tip.edge <- which(res$edge[, 1] == 5 & res$edge[, 2] <= 4)
  other.tip.edge <- which(res$edge[, 2] == 1)  # edge leading to t1
  res <- fail_edge(res, root.tip.edge)
  res <- fail_edge(res, other.tip.edge)
  expect_true("fail" %in% names(unlist(res$maps)))

  tips <- data.frame(sp = names(x), state = as.character(x),
                     stringsAsFactors = FALSE)
  tm <- matrix(c(0, 1, 1, 0), 2, 2,
               dimnames = list(c("0", "1"), c("0", "1")))
  fixed <- fix.simmap(res, tips, tm)
  # every fail segment must be gone and edge sums preserved
  expect_false("fail" %in% names(unlist(fixed$maps)))
  map.sums <- vapply(fixed$maps, sum, numeric(1))
  expect_equal(unname(map.sums), unname(fixed$edge.length),
               tolerance = 1e-8)
  # the repaired branches must END in the observed tip states
  end.state <- function(m) names(m)[length(m)]
  expect_equal(end.state(fixed$maps[[root.tip.edge]]),
               as.character(x[fixed$tip.label[fixed$edge[root.tip.edge, 2]]]))
  expect_equal(end.state(fixed$maps[[other.tip.edge]]),
               as.character(x[fixed$tip.label[fixed$edge[other.tip.edge, 2]]]))
})

test_that("fix.simmap does not leak state across multiple simulations", {
  skip_if_not_installed("phytools")
  skip_if_not_installed("igraph")
  set.seed(14)
  tree <- ape::rtree(8)
  x <- setNames(sample(0:1, 8, replace = TRUE), tree$tip.label)
  capture.output(
    res <- make.simmap2(tree, x, model = "ER", nsim = 3, message = FALSE)
  )
  # fail a different tip edge in each simulation
  tip.edges <- which(res[[1]]$edge[, 2] <= 8)
  for (i in 1:3) res[[i]] <- fail_edge(res[[i]], tip.edges[i])

  tips <- data.frame(sp = names(x), state = as.character(x),
                     stringsAsFactors = FALSE)
  tm <- matrix(c(0, 1, 1, 0), 2, 2,
               dimnames = list(c("0", "1"), c("0", "1")))
  fixed <- fix.simmap(res, tips, tm)
  for (i in seq_along(fixed)) {
    expect_false("fail" %in% names(unlist(fixed[[i]]$maps)))
    map.sums <- vapply(fixed[[i]]$maps, sum, numeric(1))
    expect_equal(unname(map.sums), unname(fixed[[i]]$edge.length),
                 tolerance = 1e-8)
  }
})
