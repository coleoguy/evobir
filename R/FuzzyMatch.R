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
  data.names <- unique(data)

  # Compute full distance matrix once (rows = data names, cols = tree tips).
  # This avoids calling adist() twice per name inside a loop.
  dist_mat <- adist(data.names, as.character(tree.names))

  # For each data name, find the closest tree tip and its distance
  min_dists <- apply(dist_mat, 1, min)
  best_idx  <- apply(dist_mat, 1, which.min)

  # Keep names that are close but NOT exact matches (distance > 0)
  close_idx <- which(min_dists <= max.dist & min_dists > 0)

  if (length(close_idx) == 0) {
    cat("Found 0 names that were close but imperfect matches\n")
  } else {
    close.taxa <- data.frame(
      name.in.data = data.names[close_idx],
      name.in.tree = tree.names[best_idx[close_idx]],
      differences  = min_dists[close_idx],
      stringsAsFactors = FALSE
    )
    cat("Found", nrow(close.taxa), "names that were close but imperfect matches\n")
    return(close.taxa)
  }
}
