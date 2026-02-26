#' Calculate Evolutionary Rates at Tree Tips
#'
#' Estimates how fast a discrete character has been evolving at each tip
#' of a phylogenetic tree.  For example, if you are studying chromosome
#' number evolution, this tells you which species have recently experienced
#' rapid changes in chromosome number on the branch leading to the present.
#'
#' The "rate" for a tip is calculated as:
#' \code{|tip_state - parent_state| / branch_length}.
#' In plain English: how much the trait changed on the last branch of the
#' tree (between the species and its most recent ancestor), divided by the
#' length of that branch (a proxy for time).  Bigger values mean faster
#' recent evolution; a rate of zero means no change occurred on that branch.
#'
#' @param tree An ultrametric phylogenetic tree of class \code{"phylo"}.
#'   "Ultrametric" means all tips are the same distance from the root
#'   (i.e., the tree represents time, not just branching pattern).
#' @param Q A transition rate matrix (square matrix) with column names
#'   matching the possible character states.  This matrix tells the model
#'   how likely each type of state change is.  It is used internally for
#'   ancestral state reconstruction (estimating what states the ancestors
#'   had).
#' @param tip.states An integer vector of observed states at the tips,
#'   with names matching the tree's tip labels.  Values should be integers
#'   from 1 to N where N is the number of possible states.  For example,
#'   if your trait has 3 states, each tip gets a value of 1, 2, or 3.
#'   If \code{NULL}, states will be derived from \code{p.mat}.
#' @param hyper Logical.  If \code{TRUE}, the model includes a binary
#'   hyperstate (state labels contain "h" which gets stripped during
#'   processing).  This is an advanced option; most users should leave
#'   it as \code{FALSE} (the default).
#' @param p.mat A probability matrix where rows are species (with row
#'   names matching tree tips) and columns are states.  Each row should
#'   have exactly one 1 and the rest 0s (i.e., you know each species'
#'   state with certainty).  This is an alternative way to provide tip
#'   states instead of \code{tip.states}.
#'
#' @return A named numeric vector of tip rates, one value per tip.
#'   Names correspond to tree tip labels.  Higher values indicate
#'   faster evolution on that terminal branch.
#'
#' @details
#' The function works in three steps:
#' \enumerate{
#'   \item \strong{Ancestral state reconstruction}: Uses
#'     \code{\link[castor]{asr_mk_model}} with the supplied transition
#'     matrix to estimate what state each internal node (ancestor)
#'     most likely had.
#'   \item \strong{Compare tip to parent}: For each species (tip), it
#'     looks at the reconstructed state of the immediately ancestral
#'     node (the parent node).
#'   \item \strong{Calculate rate}: The rate is the absolute difference
#'     between the tip state and its parent state, divided by the branch
#'     length connecting them.
#' }
#'
#' @examples
#' \dontrun{
#' library(ape)
#' library(castor)
#' tree <- rcoal(10)
#' # Make a simple 3-state transition matrix
#' Q <- matrix(c(-0.2, 0.1, 0.1,
#'                0.1, -0.2, 0.1,
#'                0.1, 0.1, -0.2), 3, 3)
#' colnames(Q) <- rownames(Q) <- c("1", "2", "3")
#' tip.states <- setNames(sample(1:3, 10, replace = TRUE), tree$tip.label)
#' rates <- GetTipRates(tree, Q, tip.states)
#' # Species with high rates evolved rapidly on their terminal branch
#' }
#'
#' @importFrom ape is.ultrametric
#' @importFrom phytools getParent
#' @importFrom castor asr_mk_model
#' @importFrom stats setNames
#' @export
GetTipRates <- function(tree = NULL,
                        Q = NULL,
                        tip.states = NULL,
                        hyper = FALSE,
                        p.mat = NULL) {

  # ---- Input validation ----
  # Tree must be ultrametric (all tips at same distance from root)
  if (is.ultrametric(tree) == FALSE) {
    stop("\n Tree is not ultrametric. Please use a formal method to ultrametricize your tree.
        If your tree is failing is.ultrametric due to rounding, you can use force.ultrametric,
        but this does not serve as a substitute for formal rate-smoothing methods.")
  }

  # Must have either tip.states or a probability matrix
  if (is.null(tip.states) && is.null(p.mat)) {
    stop("\n Tip states are not present. Please provide tip states as integers or provide
         a probability matrix to create tip states for given data.")
  }

  # Tip states must be whole numbers
  if (!is.null(tip.states)) {
    if (!is.numeric(tip.states) || any(tip.states != round(tip.states))) {
      stop("\n Tip states are not whole numbers. Please provide tip states as integers.")
    }
  }

  # If no explicit tip states given, derive them from the probability matrix
  # (each species gets the state where it has probability = 1)
  if (is.null(tip.states) && !is.null(p.mat)) {
    tip.states <- c()
    for (i in 1:nrow(p.mat)) {
      tip.states[i] <- which(p.mat[i, ] == 1)
    }
    names(tip.states) <- row.names(p.mat)
  }

  # ---- Reorder data to match tree tip order ----
  hit <- c()
  for (j in 1:length(tree$tip.label)) {
    hit[j] <- which(tree$tip.label[j] == names(tip.states))
  }
  tip.states <- tip.states[hit]

  # Verify the reordering worked
  print("Checking the order of tree tip labels and provided data.")
  order <- vapply(1:length(tree$tip.label), function(i) {
    tree$tip.label[i] == names(tip.states)[i]
  }, logical(1))
  if (sum(order) == length(tree$tip.label)) {
    print("Tree tip labels and provided data are in the correct order.")
  } else {
    print("Tree tip labels and provided data are not in the correct order.")
  }

  # ---- Ancestral state reconstruction ----
  # Use the Mk model from castor to estimate states at internal nodes
  recon <- asr_mk_model(tree = tree,
                        tip_states = tip.states,
                        transition_matrix = Q,
                        Nstates = ncol(Q),
                        include_ancestral_likelihoods = TRUE)

  # ---- Convert states back to real labels ----
  # The ASR uses column indices; restore the original state names
  colnames(recon$ancestral_likelihoods) <- colnames(Q)

  # If hyperstates are present, strip the "h" prefix from labels
  if (hyper == TRUE) {
    colnames(recon$ancestral_likelihoods) <-
      as.numeric(gsub("h", "", colnames(recon$ancestral_likelihoods)))
  }

  # Get the actual state names (e.g., chromosome numbers)
  states <- colnames(recon$ancestral_likelihoods)

  # Convert integer tip states back to actual state values
  tip.state <- as.numeric(states[tip.states])

  # ---- Get parent node states ----
  # For each internal node, pick the most likely state (highest posterior)
  est <- c()
  for (i in 1:nrow(recon$ancestral_likelihoods)) {
    est[i] <- which.max(as.vector(recon$ancestral_likelihoods[i, ]))
    names(est) <- seq(from = 1 + length(tree$tip.label),
                      to = length(tree$tip.label) + length(est))
  }

  # Map tip labels to their node numbers
  tips <- tree$tip.label
  nodes <- seq_along(tree$tip.label)
  names(nodes) <- names(tip.states)

  # Get the branch length leading to each tip
  edge.lengths <- setNames(
    tree$edge.length[sapply(nodes, function(x, y) which(y == x), y = tree$edge[, 2])],
    names(nodes)
  )

  # Find the parent node of each tip
  nodepulls <- c()
  for (i in 1:length(tree$tip.label)) {
    nodepulls[i] <- getParent(tree, nodes[i])
  }

  # Look up the state at each parent node
  node.states <- c()
  for (i in 1:length(nodepulls)) {
    node.state.idx <- which(names(est) == nodepulls[i])
    node.states[i] <- as.numeric(states[est[node.state.idx]])
  }

  # ---- Calculate tip rates ----
  # Rate = |tip_state - parent_state| / branch_length
  tip.changes <- abs(tip.state - node.states) / edge.lengths
  names(tip.changes) <- tree$tip.label

  return(tip.changes)
}
