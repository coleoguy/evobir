# Regression tests for evobir-03 (plot() on scaleTreeRates output never
# dispatched to plot.phyloscaled because the class vector was
# c("phylo", "phyloscaled"))

test_that("scaleTreeRates returns class c('phyloscaled', 'phylo')", {
  skip_if_not_installed("phytools")
  set.seed(11)
  tree <- ape::rcoal(8)
  states <- setNames(sample(1:2, 8, replace = TRUE), tree$tip.label)
  capture.output(
    scaled <- scaleTreeRates(tree, states, model = "ER",
                             nbins = 2, max.ratio = 1.5)
  )
  # "phyloscaled" must come FIRST for S3 dispatch of plot()
  expect_equal(class(scaled)[1], "phyloscaled")
  expect_true("phylo" %in% class(scaled))
  expect_equal(length(scaled$scalar), length(scaled$edge.length))
  # plot() must dispatch to plot.phyloscaled (scalars honored)
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plot(scaled, method = "multiply"))
})
