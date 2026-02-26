#' Ancestral Condition Analysis
#'
#' Tests whether transitions in a binary discrete trait tend to happen at
#' unusually high or low values of a continuous trait.  The idea is simple:
#' if evolving a horn (discrete trait = 2) requires being large (continuous
#' trait), then the ancestors that "gave birth" to horned lineages should
#' have been bigger than random ancestors.  This function checks that by
#' comparing real transition nodes against a null distribution built from
#' Monte Carlo resampling.
#'
#' @param tree A phylogenetic tree of class \code{"phylo"}.
#' @param data A data frame with exactly 3 columns:
#'   \enumerate{
#'     \item \strong{Tip labels} -- must match the labels on \code{tree}.
#'     \item \strong{Continuous trait} -- numeric values (e.g., body size).
#'     \item \strong{Discrete trait} -- coded as 1 and 2.  If one state is
#'           ancestral it \emph{must} be coded as 1.
#'   }
#'   No missing data are allowed; prune any taxa that lack data before
#'   calling this function.
#' @param mc Integer.  Number of Monte Carlo replicates used to build the
#'   null distribution (default 1000).
#' @param drop.state \code{NULL} (default) or \code{1} or \code{2}.  If set,
#'   continuous data from taxa in this state are excluded before ancestral
#'   state reconstruction.  Useful when one state is derived and you want to
#'   remove its influence on the ancestral estimates.
#' @param mat A length-4 numeric vector describing the transition rate matrix
#'   for the discrete trait.  Must be one of:
#'   \itemize{
#'     \item \code{c(0, 0, 1, 0)} -- irreversible model (only 1 -> 2 allowed)
#'     \item \code{c(0, 1, 1, 0)} -- equal rates model (1 <-> 2 same rate)
#'     \item \code{c(0, 2, 1, 0)} -- all rates different model (default)
#'   }
#' @param pi Root state prior probabilities.  Options:
#'   \itemize{
#'     \item \code{"equal"} (default) -- equal probability for each state.
#'     \item \code{"estimated"} -- estimated from the data.
#'     \item A numeric vector of length 2 that sums to 1.
#'   }
#' @param n.tails Either 1 (default, one-tailed test) or 2 (two-tailed test).
#'   Use 1 if you have an a priori hypothesis about the direction.
#' @param message Logical.  If \code{TRUE} (default), results are printed
#'   to the console.
#'
#' @return A named list.  For the irreversible model (\code{sum(mat) <= 1})
#'   the list contains:
#'   \describe{
#'     \item{OriginatingVals}{Mean continuous trait value at nodes where
#'       trait 2 originated.}
#'     \item{NTrans}{Number of 1 -> 2 transitions detected.}
#'     \item{NullDist}{Numeric vector of length \code{mc} giving the null
#'       distribution of mean values.}
#'     \item{pval}{P-value from the Monte Carlo test.}
#'   }
#'   For the reversible model (\code{sum(mat) > 1}), the list has
#'   analogous elements for both directions (1->2 and 2->1).
#'
#' @details
#' The function works in three steps:
#' \enumerate{
#'   \item Reconstruct ancestral continuous values via maximum likelihood
#'         under a Brownian motion model (\code{geiger::anc.ML}).
#'   \item Map discrete trait transitions onto the tree using stochastic
#'         character mapping (\code{\link[phytools]{make.simmap}}).
#'   \item Build a null distribution by randomly sampling the same number
#'         of nodes from the pool of eligible ancestors, repeating
#'         \code{mc} times.
#' }
#'
#' @references
#' Blackmon, H., Adams, R.H. (2015). evobiR: Comparative analyses and
#' teaching tools for evolutionary biology.
#'
#' @examples
#' \dontrun{
#' library(ape)
#' library(phytools)
#' # Simulate a tree and data
#' tree <- rcoal(30)
#' data <- data.frame(
#'   species = tree$tip.label,
#'   body_size = rnorm(30, 10, 2),
#'   has_horn = sample(1:2, 30, replace = TRUE)
#' )
#' result <- AncCond(tree, data, mc = 500)
#' }
#'
#' @importFrom phytools anc.ML make.simmap describe.simmap
#' @importFrom stats sd
#' @export
AncCond <- function(tree, data, mc = 1000, drop.state = NULL,
                    mat = c(0, 2, 1, 0), pi = "equal",
                    n.tails = 1, message = TRUE) {

  # ---- Input validation ----
  # Make sure the user provided correct data types before doing any analysis
  if (!inherits(tree, 'phylo')) {
    stop('tree must be class phylo')
  }
  if (!is.data.frame(data) | ncol(data) != 3) {
    stop('data should be a dataframe with 3 columns\n(tip labels, cont data, discrete data)')
  }
  if (!is.numeric(mc) | round(mc) != mc) {
    stop('mc should be a numeric positive integer')
  }
  if (!is.null(drop.state)) {
    if (!drop.state %in% c(1, 2)) {
      stop('drop.state must be NULL, or numeric 1 or 2')
    }
  }
  if (!sum(mat == c(0, 0, 1, 0)) == 4 &
      !sum(mat == c(0, 1, 1, 0)) == 4 &
      !sum(mat == c(0, 2, 1, 0)) == 4) {
    stop('mat must be a vector of the form c(0,0,1,0), c(0,1,1,0), or c(0,2,1,0)')
  }
  if ((!pi %in% c('equal', 'estimated'))[1]) {
    if (!is.numeric(pi)) {
      stop('pi must be equal, estimated or a vector of length 2\nwith probabilities for the state of the discrete character at the root')
    }
    if (length(pi) != 2 | sum(pi) != 1) {
      stop('pi must be equal, estimated or a vector of length 2\nwith probabilities for the state of the discrete character at the root')
    }
  }
  if (n.tails != 1 & n.tails != 2) {
    stop('n.tails should be numeric 1 or 2')
  }

  # ---- Prepare trait vectors ----
  # Create a named vector for the discrete trait (needed by make.simmap)
  dt.vec <- data[, 3]
  names(dt.vec) <- data[, 1]

  # Create a named vector for the continuous trait
  # If drop.state is set, exclude taxa in that state so they don't
  # bias the Brownian motion reconstruction
  if (!is.null(drop.state)) {
    ct.data <- data[(data[, 3] != drop.state), ]
  } else {
    ct.data <- data
  }
  ct.vec <- as.numeric(ct.data[, 2])
  names(ct.vec) <- ct.data[, 1]

  # Check for missing data -- the analysis can't handle NAs
  if (sum(is.na(ct.vec)) > 0 | sum(is.na(dt.vec)) > 0) {
    stop('There exists missing trait data for some species in the phylogeny.\n
         Please remove such taxa from the tree.')
  }

  # ---- Ancestral State Reconstruction ----
  # Continuous trait: ML reconstruction under Brownian motion
  anc.states.cont.trait <- anc.ML(tree, ct.vec, model = "BM")

  # Discrete trait: stochastic character mapping (1 simulation)
  # This places trait changes at specific points along branches
  anc.state.dt <- make.simmap(tree, dt.vec,
                              model = matrix(mat, 2),
                              nsim = 1,
                              pi = pi,
                              message = FALSE)

  # ---- Find transition nodes ----
  # mapped.edge has time in each state per branch; branches with changes
  # have positive values in BOTH columns
  ss_nodes <- anc.state.dt$mapped.edge[, 1] > 0 &
    anc.state.dt$mapped.edge[, 2] > 0

  # Get the node-pair labels for branches that had transitions
  wanted_branches <- ss_nodes[ss_nodes == TRUE]
  wanted_nodes <- names(wanted_branches)

  # ---- Reversible model (both 1->2 and 2->1 transitions possible) ----
  if (sum(mat) > 1) {
    # Separate producing nodes by direction of transition
    producing.nodes12 <- c()
    producing.nodes21 <- c()
    trans.maps <- anc.state.dt$maps[ss_nodes == TRUE]

    # Extract the rootward (parent) node from each branch label
    wanted_nodes <- gsub(",.*", "", wanted_nodes)

    # Classify each transition as 1->2 or 2->1 based on starting state
    for (i in 1:length(wanted_nodes)) {
      if (names(trans.maps[[i]])[1] == '1') {
        producing.nodes12 <- c(producing.nodes12, wanted_nodes[i])
      } else if (names(trans.maps[[i]])[1] == '2') {
        producing.nodes21 <- c(producing.nodes21, wanted_nodes[i])
      }
    }
    # Remove duplicates (one node could appear on multiple branches)
    producing.nodes12 <- unique(producing.nodes12)
    producing.nodes21 <- unique(producing.nodes21)

    # Get the mean continuous trait value at transition nodes
    anc.states <- anc.states.cont.trait
    orig.val12 <- mean(anc.states$ace[names(anc.states$ace) %in%
                                        producing.nodes12])
    orig.val21 <- mean(anc.states$ace[names(anc.states$ace) %in%
                                        producing.nodes21])

    # ---- Build null distributions ----
    # For each direction, resample from the eligible node pool
    null.orig.val12 <- vector(length = mc)
    null.orig.val21 <- vector(length = mc)
    number.of.trans12 <- length(producing.nodes12)
    number.of.trans21 <- length(producing.nodes21)

    # Identify which nodes are eligible for each direction
    # (nodes not in the "wrong" state)
    node.states <- describe.simmap(anc.state.dt)$states
    anc.cond.nodes12 <- anc.states.cont.trait$ace[
      names(anc.states.cont.trait$ace) %in%
        names(node.states)[node.states != '2']]
    anc.cond.nodes21 <- anc.states.cont.trait$ace[
      names(anc.states.cont.trait$ace) %in%
        names(node.states)[node.states != '1']]

    # Monte Carlo: randomly pick the same number of nodes, take the mean
    for (j in 1:mc) {
      null.orig.val12[j] <- mean(sample(anc.cond.nodes12,
                                        length(producing.nodes12)))
      null.orig.val21[j] <- mean(sample(anc.cond.nodes21,
                                        length(producing.nodes21)))
    }

    # ---- Calculate p-values ----
    # Count how many null values are as extreme as the observed value
    bigger12 <- (sum(null.orig.val12 >= orig.val12) / mc)
    smaller12 <- (sum(null.orig.val12 <= orig.val12) / mc)
    if (!is.null(producing.nodes12)) {
      if (bigger12 <= smaller12) pval12 <- bigger12
      if (smaller12 < bigger12) pval12 <- smaller12
      if (n.tails == 2) pval12 <- 2 * pval12
    } else {
      pval12 <- NA
    }

    bigger21 <- (sum(null.orig.val21 >= orig.val21) / mc)
    smaller21 <- (sum(null.orig.val21 <= orig.val21) / mc)
    if (!is.null(producing.nodes21)) {
      if (bigger21 <= smaller21) pval21 <- bigger21
      if (smaller21 < bigger21) pval21 <- smaller21
      if (n.tails == 2) pval21 <- 2 * pval21
    } else {
      pval21 <- NA
    }

    # ---- Print results ----
    if (message == TRUE) {
      cat(paste("Mean value for the continuous trait at origin of trait 2:",
                round(orig.val12, digits = 4), "\n"))
      cat(paste("Mean value for the continuous trait at origin of trait 1:",
                round(orig.val21, digits = 4), "\n"))
      cat(paste("Number of producing nodes 1->2:",
                round(mean(number.of.trans12), digits = 4), "\n"))
      cat(paste("Number of producing nodes 2->1:",
                round(mean(number.of.trans21), digits = 4), "\n"))
      cat(paste("Mean of null dist 1->2:",
                round(mean(null.orig.val12), digits = 4), "\n"))
      cat(paste("Mean of null dist 2->1:",
                round(mean(null.orig.val21), digits = 4), "\n"))
      cat(paste("SD of null dist 1->2:",
                round(sd(null.orig.val12), digits = 4), "\n"))
      cat(paste("SD of null dist 2->1:",
                round(sd(null.orig.val21), digits = 4), "\n"))
      cat(paste("pvalue 1->2:", round(pval12, digits = 4), "\n"))
      cat(paste("pvalue 2->1:", round(pval21, digits = 4), "\n\n"))
      if (is.null(producing.nodes12)) {
        cat('No 1 -> 2 transitions occured NA and NaN values produced.')
      }
      if (is.null(producing.nodes21)) {
        cat('No 2 -> 1 transitions occured NA and NaN values produced.')
      }
    }

    # Package results into a named list
    results <- list(
      "OriginatingVals1->2" = orig.val12,
      "NTrans1->2"          = number.of.trans12,
      "NullDist1->2"        = null.orig.val12,
      "pval1->2"            = pval12,
      "OriginatingVals2->1" = orig.val21,
      "NTrans2->1"          = number.of.trans21,
      "NullDist2->1"        = null.orig.val21,
      "pval2->1"            = pval21
    )

  } else {
    # ---- Irreversible model (only 1->2 transitions) ----
    wanted_nodes <- gsub(",.*", "", wanted_nodes)
    producing.nodes <- unique(wanted_nodes)

    # Mean continuous trait at transition nodes
    anc.states <- anc.states.cont.trait
    orig.val <- mean(anc.states$ace[names(anc.states$ace) %in%
                                      producing.nodes])

    # Build null distribution
    null.orig.val <- vector(length = mc)
    number.of.trans <- length(producing.nodes)
    node.states <- describe.simmap(anc.state.dt)$states
    anc.cond.nodes <- anc.states.cont.trait$ace[
      names(anc.states.cont.trait$ace) %in%
        names(node.states)[node.states != '2']]

    for (j in 1:mc) {
      null.orig.val[j] <- mean(sample(anc.cond.nodes,
                                      length(producing.nodes)))
    }

    # Calculate p-value
    bigger <- (sum(null.orig.val >= orig.val) / mc)
    smaller <- (sum(null.orig.val <= orig.val) / mc)
    if (bigger <= smaller) pval <- bigger
    if (smaller < bigger) pval <- smaller
    if (n.tails == 2) pval <- 2 * pval

    # Print results
    if (message == TRUE) {
      cat(paste("Mean value for the continuous trait at origin of derived trait:",
                round(orig.val, digits = 4), "\n"))
      cat(paste("Number of producing nodes:",
                round(mean(number.of.trans), digits = 4), "\n"))
      cat(paste("Mean of null dist:",
                round(mean(null.orig.val), digits = 4), "\n"))
      cat(paste("SD of null dist:",
                round(sd(null.orig.val), digits = 4), "\n"))
      cat(paste("pvalue:", round(pval, digits = 4), "\n\n"))
    }

    # Package results
    results <- list(
      OriginatingVals = orig.val,
      NTrans          = number.of.trans,
      NullDist        = null.orig.val,
      pval            = pval
    )
  }

  return(results)
}
