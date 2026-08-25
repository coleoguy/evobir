#' Count Tree Topologies in a Collection
#'
#' Classifies each tree in a collection against a set of reference
#' topologies and counts how many trees match each reference.
#'
#' \strong{What is a "topology"?}  A topology is the branching pattern
#' of a tree -- which species are most closely related to which --
#' ignoring branch lengths.  For example, with three species (A, B, C),
#' there are three possible topologies: ((A,B),C), ((A,C),B), and
#' ((B,C),A).  This function counts how many trees in your collection
#' have each branching pattern.
#'
#' \strong{When would you use this?}  After a Bayesian phylogenetic
#' analysis (e.g., MrBayes, BEAST), you have thousands of sampled trees.
#' This function tells you what fraction of those trees support each
#' possible branching arrangement, helping you assess which topology
#' has the most support in the posterior distribution.
#'
#' @param collection Character.  Path to a Newick format file containing
#'   the collection of trees to classify (e.g., posterior trees from
#'   an MCMC analysis).
#' @param ref Character.  Path to a Newick format file containing the
#'   reference topologies to count.  Each tree in this file represents
#'   one possible branching pattern you want to test.
#' @param classes Logical.  If \code{TRUE} (default), returns both the
#'   counts per topology AND the classification of each individual tree.
#'   If \code{FALSE}, returns only the counts.
#' @param verbose Logical.  If \code{TRUE} (default), prints a warning
#'   listing any trees that don't match any reference topology.
#'
#' @return If \code{classes = TRUE}: a list with two elements:
#'   \describe{
#'     \item{[[1]]}{A numeric vector of counts -- how many trees matched
#'       each reference topology.}
#'     \item{[[2]]}{A numeric vector giving the topology class (reference
#'       number) for each tree in the collection.  \code{NA} means the
#'       tree didn't match any of the reference topologies you provided.}
#'   }
#'   If \code{classes = FALSE}: just the count vector.
#'
#' @examples
#' \dontrun{
#' # You have 1000 posterior trees and 3 candidate topologies.
#' # How many trees support each arrangement?
#' result <- countTrees("posterior_trees.nwk", "candidate_topologies.nwk")
#' result[[1]]  # e.g., c(712, 201, 87) = topology 1 has most support
#' result[[2]]  # which topology each individual tree matched
#' }
#'
#' @importFrom ape read.tree all.equal.phylo
#' @export
countTrees <- function(collection = NULL, ref = NULL, classes = TRUE,
                       verbose = TRUE) {
  # ---- Input validation ----
  if (is.null(collection)) {
    stop("supply a path to a collection of trees in a Newick format file")
  }
  if (is.null(ref)) {
    stop("supply a path to a set of topologies to count in a Newick format file")
  }

  # Read tree collections
  trees <- read.tree(collection)
  types <- read.tree(ref)
  # Normalize single-tree files: read.tree returns a "phylo" object for
  # one tree; wrap it so trees[[i]]/types[[i]] index whole trees rather
  # than tree components
  if (inherits(trees, "phylo")) trees <- c(trees)
  if (inherits(types, "phylo")) types <- c(types)
  top.num <- length(types)

  # Storage for classifications
  tree.class <- as.numeric(rep(NA, times = length(trees)))
  classification <- vector(length = top.num)
  bad <- c()  # tracks unmatched trees

  # ---- Compare each tree to each reference topology ----
  for (i in 1:length(trees)) {
    missing <- TRUE
    counter <- 1
    while (missing) {
      # Compare topology only (ignore branch lengths)
      x <- all.equal.phylo(trees[[i]], types[[counter]],
                           use.edge.length = FALSE)
      if (x) {
        # Match found -- record it and move on
        classification[counter] <- classification[counter] + 1
        tree.class[i] <- counter
        missing <- FALSE
      } else if (counter == top.num) {
        # No match among any reference topology
        bad <- c(bad, i)
        missing <- FALSE
      }
      counter <- counter + 1
    }
  }

  # Report unmatched trees
  if (verbose) {
    if (length(bad) > 0) {
      print(paste("Some trees do not match available topologies. You may want to check trees:", bad))
    }
  }

  # Return results
  if (classes == TRUE) {
    classification <- list(classification, tree.class)
    return(classification)
  } else {
    return(classification)
  }
}
