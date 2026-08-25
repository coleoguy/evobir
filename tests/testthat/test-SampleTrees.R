# Regression tests for evobir-04 (SampleTrees burn-in off-by-one,
# fractional indices, and broken output filenames)

test_that("SampleTrees writes correctly named output files", {
  nex <- system.file("trees.nex", package = "evobiR")
  old.wd <- setwd(tempdir())
  on.exit(setwd(old.wd))
  capture.output(
    SampleTrees(trees = nex, burnin = 0.25, final.number = 10,
                format = "new", prefix = "sampled_trees")
  )
  # old code wrote "sampled_trees .nwk" (with a space)
  expect_true(file.exists("sampled_trees.nwk"))
  expect_false(file.exists("sampled_trees .nwk"))
  out <- ape::read.tree("sampled_trees.nwk")
  n.out <- if (inherits(out, "phylo")) 1L else length(out)
  expect_equal(n.out, 10L)
  unlink("sampled_trees.nwk")
})

test_that("SampleTrees discards the full burn-in fraction", {
  nex <- system.file("trees.nex", package = "evobiR")
  # trees.nex ships 100 trees; burnin = 0.25 leaves exactly 75.
  # Asking for all 75 must work ...
  old.wd <- setwd(tempdir())
  on.exit(setwd(old.wd))
  capture.output(
    SampleTrees(trees = nex, burnin = 0.25, final.number = 75,
                format = "new", prefix = "all_post_burnin")
  )
  out <- ape::read.tree("all_post_burnin.nwk")
  expect_equal(length(out), 75L)
  unlink("all_post_burnin.nwk")
  # ... and asking for 76 must fail (the old code retained the last
  # burn-in tree, so 76 would have "worked")
  expect_error(
    capture.output(
      SampleTrees(trees = nex, burnin = 0.25, final.number = 76,
                  format = "new", prefix = "too_many")
    ),
    "exceeds"
  )
})

test_that("SampleTrees validates burnin", {
  nex <- system.file("trees.nex", package = "evobiR")
  expect_error(
    SampleTrees(trees = nex, burnin = 1, final.number = 5,
                format = "new", prefix = "bad"),
    "burnin"
  )
})
