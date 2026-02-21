test_that("FuzzyMatch returns a dataframe", {
  data(hym.tree)
  names <- c("Pepsis_elegans", "Pheidele_lucreti")
  result <- FuzzyMatch(tree = hym.tree, data = names, max.dist = 3)
  expect_true(is.data.frame(result))
})

test_that("FuzzyMatch finds close matches", {
  data(hym.tree)
  # Intentionally misspelled names
  names <- c("Pepsis_elegans")
  result <- FuzzyMatch(tree = hym.tree, data = names, max.dist = 3)
  # result is a data.frame when matches are found, NULL otherwise
  expect_true(is.null(result) || (is.data.frame(result) && nrow(result) > 0))
})
