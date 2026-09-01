# Tests for AncCond with missing data: x and y observed on different
# (possibly disjoint) subsets of tips.

# Independent BM conditional expectation for every node of the tree given
# the observed tips, via the full covariance matrix (C[i,j] = height of the
# MRCA of i and j) and a GLS root estimate. This is what fastAnc computes at
# internal nodes and what .anc.subset must reproduce for on-path and
# dataless nodes.
bm.cond.expect <- function(tree, x.obs) {
  n.tips  <- length(tree$tip.label)
  n.total <- n.tips + tree$Nnode
  h  <- ape::node.depth.edgelength(tree)
  M  <- ape::mrca(tree, full = TRUE)
  C  <- matrix(h[M], n.total, n.total)
  obs <- match(names(x.obs), tree$tip.label)
  Coo <- C[obs, obs]
  Cinv <- solve(Coo)
  one <- rep(1, length(obs))
  mu  <- as.numeric((one %*% Cinv %*% x.obs) / (one %*% Cinv %*% one))
  as.numeric(mu + C[, obs] %*% Cinv %*% (x.obs - mu))
}

test_that(".pruning2 marginalises NA tips exactly (equals dropping the tip)", {
  skip_if_not_installed("phytools")
  set.seed(3)
  tree <- ape::rtree(20)
  y <- setNames(sample(0:1, 20, replace = TRUE), tree$tip.label)
  miss <- c(3, 8, 15)
  y.na <- y; y.na[miss] <- NA
  q01 <- 0.7; q10 <- 0.3

  ll.na <- .pruning2(reorder(tree, "postorder"), unname(y.na), q01, q10)$logL
  tree.d <- ape::drop.tip(tree, tree$tip.label[miss])
  y.d <- y[tree.d$tip.label]
  ll.d <- .pruning2(reorder(tree.d, "postorder"), unname(y.d), q01, q10)$logL
  expect_equal(ll.na, ll.d, tolerance = 1e-10)
})

test_that(".anc.subset with complete data equals fastAnc", {
  skip_if_not_installed("phytools")
  set.seed(4)
  tree <- ape::rtree(15)
  x <- setNames(rnorm(15), tree$tip.label)
  anc <- .anc.subset(tree, x)
  fa  <- phytools::fastAnc(tree, x)
  expect_equal(anc[1:15], unname(x))
  expect_equal(anc[16:29], unname(fa[as.character(16:29)]), tolerance = 1e-8)
})

test_that(".anc.subset returns the BM conditional expectation at every node", {
  skip_if_not_installed("phytools")
  set.seed(5)
  for (rep.i in 1:5) {
    tree <- ape::rtree(25)
    x <- setNames(phytools::fastBM(tree), tree$tip.label)
    keep <- sort(sample(25, 11))
    x.obs <- x[keep]
    anc <- .anc.subset(tree, x.obs)
    expect_false(anyNA(anc))
    expect_equal(anc, bm.cond.expect(tree, x.obs), tolerance = 1e-6)
  }
})

test_that(".anc.subset handles data confined to one clade (nodes above pruned root)", {
  skip_if_not_installed("phytools")
  set.seed(6)
  tree <- ape::rtree(20)
  x <- setNames(rnorm(20), tree$tip.label)
  # keep only the tips of one non-root clade
  root <- 21L
  kid <- tree$edge[tree$edge[, 1] == root, 2]
  kid <- kid[kid > 20][1]
  if (is.na(kid)) skip("no internal child of root")
  clade.tips <- ape::extract.clade(tree, kid)$tip.label
  if (length(clade.tips) < 3) skip("clade too small")
  anc <- .anc.subset(tree, x[clade.tips])
  expect_equal(anc, bm.cond.expect(tree, x[clade.tips]), tolerance = 1e-6)
  # the full root inherits the pruned-root estimate
  expect_equal(anc[root], anc[kid])
})

test_that("AncCond runs with x and y observed on disjoint tip sets", {
  skip_if_not_installed("phytools")
  set.seed(7)
  tree <- phytools::pbtree(n = 60, scale = 1)
  x <- phytools::fastBM(tree)
  y <- setNames(rbinom(60, 1, 0.4), tree$tip.label)
  x.sub <- x[1:35]
  y.sub <- y[30:60]
  res <- AncCond(tree, x.sub, y.sub, n.maps = 10, n.sims = 30)
  expect_s3_class(res, "anccond")
  expect_equal(unname(res$n.obs), c(60, 35, 31))
  expect_true(is.finite(res$T.obs.01) || is.finite(res$T.obs.10))
  expect_true(all(res$p.01 >= 0 & res$p.01 <= 1, na.rm = TRUE))
  expect_true(all(res$p.10 >= 0 & res$p.10 <= 1, na.rm = TRUE))
  expect_equal(length(res$T.null.01), 30)
  # NA-coded and absent tips are treated identically
  x.na <- x; x.na[36:60] <- NA
  y.na <- y; y.na[1:29] <- NA
  set.seed(11); r1 <- AncCond(tree, x.sub, y.sub, n.maps = 5, n.sims = 10)
  set.seed(11); r2 <- AncCond(tree, x.na,  y.na,  n.maps = 5, n.sims = 10)
  expect_equal(r1$T.obs.01, r2$T.obs.01)
  expect_equal(r1$T.null.01, r2$T.null.01)
})

test_that("AncCond prune=TRUE excludes unknown-y tips from the reconstruction", {
  skip_if_not_installed("phytools")
  set.seed(8)
  tree <- phytools::pbtree(n = 40, scale = 1)
  x <- phytools::fastBM(tree)
  y <- setNames(rbinom(40, 1, 0.4), tree$tip.label)
  y[1:5] <- NA
  res <- AncCond(tree, x, y, n.maps = 5, n.sims = 10, prune = TRUE)
  expect_s3_class(res, "anccond")
  expect_true(res$prune)
})

test_that("AncCond errors informatively on unusable data", {
  skip_if_not_installed("phytools")
  set.seed(9)
  tree <- ape::rtree(10)
  x <- setNames(rnorm(10), tree$tip.label)
  y <- setNames(rep(0:1, 5), tree$tip.label)
  expect_error(AncCond(tree, x[1:2], y), "At least 3 tips")
  expect_error(AncCond(tree, x, setNames(rep(0L, 10), tree$tip.label)),
               "Both states")
  expect_error(AncCond(tree, setNames(x, paste0("z", 1:10)), y),
               "do not match")
})
