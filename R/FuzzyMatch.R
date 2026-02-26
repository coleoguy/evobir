#' Fuzzy Name Matching Between Tree and Data
#'
#' Finds taxon names in your data that are close (but not exact) matches
#' to the tip labels on your phylogenetic tree.
#'
#' \strong{Why do you need this?}  In comparative biology, you often
#' combine a phylogenetic tree (from one source) with trait data (from
#' another source).  If a species name is spelled slightly differently
#' in the two datasets -- for example, "Homo_sapiens" in the tree but
#' "Homo_sapien" in your data -- that species will be silently dropped
#' from your analysis with no warning.  This function catches those
#' near-misses so you can fix them before running your analysis.
#'
#' The function uses Levenshtein distance (the minimum number of
#' single-character edits -- insertions, deletions, or substitutions --
#' needed to change one name into another) to measure similarity.
#'
#' @param tree A phylogenetic tree of class \code{"phylo"} or
#'   \code{"multiPhylo"}.  If a \code{multiPhylo} object is provided, only
#'   the first tree is used (with a warning).
#' @param data A character vector of taxon names from your dataset
#'   (e.g., the species column of your data frame).
#' @param max.dist Integer.  Maximum number of character differences to
#'   still count as a "close match."  For example, \code{max.dist = 2}
#'   will flag name pairs that differ by 1 or 2 characters.
#'
#' @return A data frame with three columns:
#'   \describe{
#'     \item{name.in.data}{The name as it appears in your data.}
#'     \item{name.in.tree}{The closest matching name on the tree.}
#'     \item{differences}{Number of character differences (Levenshtein
#'       distance) between the two names.}
#'   }
#'   Returns \code{NULL} invisibly if no close matches are found.
#'
#' @examples
#' \dontrun{
#' library(ape)
#' tree <- read.tree(text = "((Homo_sapiens, Pan_troglodytes), Gorilla_gorilla);")
#' my_data <- c("Homo_sapien", "Pan_troglodytes", "Gorila_gorilla")
#' FuzzyMatch(tree, my_data, max.dist = 2)
#' # Should flag "Homo_sapien" and "Gorila_gorilla" as close matches
#' }
#'
#' @importFrom utils adist
#' @export
FuzzyMatch <- function(tree, data, max.dist) {
  # If user passed a collection of trees, just use the first one
  if (inherits(tree, "multiPhylo")) {
    warning("Multiple trees were supplied only first is being used")
    tree <- tree[[1]]
  }

  tree.names <- tree$tip.label
  close.taxa <- data.frame()
  counter <- 1
  data.names <- unique(data)

  # Compare each name in the data against all tree tip labels
  for (i in 1:length(data.names)) {
    # adist() computes Levenshtein distance between strings
    name.dist <- min(adist(data.names[i], as.character(tree.names)))
    # Keep it if it's close but NOT an exact match (distance > 0)
    if (name.dist <= max.dist & name.dist > 0) {
      best.match <- which.min(adist(data.names[i], as.character(tree.names)))
      close.taxa[counter, 1] <- data.names[i]
      close.taxa[counter, 2] <- tree.names[best.match]
      close.taxa[counter, 3] <- name.dist
      counter <- counter + 1
    }
  }

  # Report how many close matches were found
  if (counter == 1) {
    cat("Found", counter - 1, "names that were close but imperfect matches\n")
  }
  if (counter > 1) {
    cat("Found", counter - 1, "names that were close but imperfect matches\n")
    colnames(close.taxa) <- c("name.in.data", "name.in.tree", "differences")
    return(close.taxa)
  }
}
