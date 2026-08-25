# Regression tests for evobir-06 (GetTipRates: `if (NA)` crash on
# unmatched tip labels; proceeded on misordered data after a print)

test_that("GetTipRates errors informatively on unmatched tip labels", {
  skip_if_not_installed("castor")
  set.seed(5)
  tree <- ape::rcoal(10)
  Q <- matrix(c(-0.2, 0.1, 0.1,
                 0.1, -0.2, 0.1,
                 0.1, 0.1, -0.2), 3, 3)
  colnames(Q) <- rownames(Q) <- c("1", "2", "3")
  tip.states <- setNames(sample(1:3, 10, replace = TRUE), tree$tip.label)
  # rename one entry so it no longer matches a tip label
  names(tip.states)[1] <- "not_a_tip"
  expect_error(
    GetTipRates(tree, Q, tip.states),
    "missing"
  )
})

test_that("GetTipRates works on correctly named data", {
  skip_if_not_installed("castor")
  set.seed(6)
  tree <- ape::rcoal(10)
  Q <- matrix(c(-0.2, 0.1, 0.1,
                 0.1, -0.2, 0.1,
                 0.1, 0.1, -0.2), 3, 3)
  colnames(Q) <- rownames(Q) <- c("1", "2", "3")
  tip.states <- setNames(sample(1:3, 10, replace = TRUE), tree$tip.label)
  # scramble order: function must reorder by name, not position
  tip.states <- tip.states[sample(seq_along(tip.states))]
  rates <- GetTipRates(tree, Q, tip.states)
  expect_equal(length(rates), 10L)
  expect_equal(names(rates), tree$tip.label)
  expect_true(all(is.finite(rates)))
})
