#' Fix Failed Branches in Stochastic Mapping
#'
#' Repairs branches that failed during stochastic character mapping
#' (from \code{\link{make.simmap2}} or \code{phytools::make.simmap}).
#'
#' \strong{What is a "failed" branch?}  Stochastic character mapping
#' works by simulating trait evolution along each branch of a tree.
#' For some branches, the simulation can't find a valid history that
#' starts and ends at the correct states -- like trying to find a
#' path through a maze with too many dead ends.  After many attempts
#' (controlled by \code{rejmax} in \code{\link{make.simmap2}}), the
#' branch is marked as "fail" instead of running forever.
#'
#' This function repairs those failed branches by using graph theory
#' (shortest path algorithm) to find the simplest valid sequence of
#' state transitions between the start and end states, then divides
#' the branch time equally among the required intermediate states.
#'
#' @param hists A \code{simmap} or \code{multiSimmap} object (the output
#'   of \code{\link{make.simmap2}} or \code{\link[phytools]{make.simmap}})
#'   that may contain failed branches.
#' @param tips A two-column data frame where:
#'   \describe{
#'     \item{Column 1}{Species names matching the tree tip labels.}
#'     \item{Column 2}{The known discrete state for each species.}
#'   }
#' @param transition.matrix A symmetric matrix describing which state
#'   transitions are possible.  Non-zero entries indicate allowed
#'   transitions.  This is used to find valid paths between states
#'   using shortest-path algorithms.
#'
#' @return A fixed \code{simmap} or \code{multiSimmap} object with the
#'   failed branches repaired.  The \code{maps} and \code{mapped.edge}
#'   components are updated to reflect the inferred state transitions.
#'
#' @details
#' For each failed branch, the function:
#' \enumerate{
#'   \item Identifies the start state (from the parent node) and end
#'     state (from the data or child node).
#'   \item Finds the shortest valid path between these states using
#'     the transition matrix as a graph.
#'   \item Divides the branch length equally among all states in the
#'     path.
#'   \item Updates both the \code{maps} list and the \code{mapped.edge}
#'     matrix.
#' }
#'
#' @examples
#' \dontrun{
#' # Run stochastic mapping with a rejection limit
#' mapped_tree <- make.simmap2(tree, states, model = "ARD",
#'                              nsim = 10, rejmax = 100000)
#' # Fix any branches that failed
#' fixed_tree <- fix.simmap(mapped_tree, tip_data, trans_matrix)
#' }
#'
#' @importFrom igraph graph_from_data_frame shortest_paths
#' @export
fix.simmap <- function(hists, tips, transition.matrix) {

  # ---- Handle single vs. multiple maps ----
  # If a single simmap, wrap it in a list for uniform processing
  if ("simmap" %in% class(hists)) {
    hist <- hists
    hists <- list()
    hists[[1]] <- hist
  }

  # ---- Catalog all failed edges across all simulations ----
  fail.list <- rep(list(c()), length(hists))

  for (i in 1:length(hists)) {
    for (j in 1:nrow(hists[[i]]$edge)) {
      if ("fail" %in% names(hists[[i]]$maps[[j]])) {
        fail.list[[i]] <- c(fail.list[[i]], j)
      } else {
        names(hists[[i]]$maps[[j]])
      }
    }
  }

  # ---- Fix each simulation ----
  for (i in 1:length(hists)) {

    # Skip simulations with no failures
    if (length(fail.list[[i]]) == 0) {
      next
    }

    # Reset the working table so state never leaks between simulations
    branches <- NULL

    # ---- Catalog tip states from the current simulation ----
    tip.states <- c()
    for (j in 1:length(hists[[i]]$tip.label)) {
      tip.leading.edge <- which(hists[[i]]$edge[, 2] == j)
      tip.state <- names(hists[[i]]$maps[[tip.leading.edge]]
                         [length(hists[[i]]$maps[[tip.leading.edge]])])
      tip.states <- c(tip.states, tip.state)
    }
    names(tip.states) <- hists[[i]]$tip.label

    # ---- Identify failed tips ----
    # These are tips where the simulation couldn't resolve the state
    fail.taxa <- tips[which(tips[, 1] %in%
                              names(which(tip.states == "fail"))), ]

    if (nrow(fail.taxa) != 0) {
      # Add tip numbers to fail.taxa
      for (j in 1:nrow(fail.taxa)) {
        fail.taxa[j, 3] <- which(hists[[i]]$tip.label %in%
                                   fail.taxa[j, 1])
      }

      fail.tips.names <- which(tip.states == "fail")
      fail.tips.numbers <- matrix(hists[[i]]$edge[fail.list[[i]], ],
                                  ncol = 2)
      fail.tips.edges <- fail.list[[i]][which(fail.tips.numbers[, 2]
                                              <= length(hists[[i]]$tip.label))]
      fail.tips.numbers <- matrix(fail.tips.numbers[which(fail.tips.numbers[, 2]
                                                          <= length(hists[[i]]$tip.label)), ],
                                  ncol = 2)
      fail.tips.start <- fail.tips.numbers[, 1]

      # Get the state at the start of each failed tip branch.  Handle
      # parents that are the root (which have no leading edge) so that
      # names stay aligned with start nodes.
      root.node <- length(hists[[i]]$tip.label) + 1
      unique.starts <- unique(fail.tips.start)
      fail.tips.start.states <- character(length(unique.starts))
      for (j in 1:length(unique.starts)) {
        if (unique.starts[j] != root.node) {
          leading.edge <- which(hists[[i]]$edge[, 2] == unique.starts[j])
          leading.map <- hists[[i]]$maps[[leading.edge]]
          fail.tips.start.states[j] <- names(leading.map)[length(leading.map)]
        } else {
          # Root case: use the first state of a non-failed edge
          # descending from the root
          root.edges <- which(hists[[i]]$edge[, 1] == root.node)
          state.found <- NA_character_
          for (re in root.edges) {
            nm <- names(hists[[i]]$maps[[re]])[1]
            if (nm != "fail") {
              state.found <- nm
              break
            }
          }
          fail.tips.start.states[j] <- state.found
        }
      }
      names(fail.tips.start.states) <- unique.starts

      # Get the desired end state from the tip data
      fail.tips.end.states <- fail.taxa[, 2]
      names(fail.tips.end.states) <- fail.taxa[, 3]

      # Get the branch lengths of failed tip branches
      branch.lengths <- c()
      for (k in 1:nrow(fail.tips.numbers)) {
        branch.lengths[k] <- hists[[i]]$edge.length[
          which(fail.tips.numbers[k, 1] == hists[[i]]$edge[, 1]
                & fail.tips.numbers[k, 2] == hists[[i]]$edge[, 2])]
      }
      branches <- cbind(fail.tips.numbers, branch.lengths)
    } else {
      fail.tips.edges <- NA
    }

    # ---- Identify failed internal edges ----
    fail.internal.edges <- fail.list[[i]][which(!fail.list[[i]] %in%
                                                  fail.tips.edges)]

    if (length(fail.internal.edges) != 0) {
      fail.internal <- as.numeric(hists[[i]]$edge[fail.internal.edges, 2])
      fail.internal.start <- hists[[i]]$edge[
        which(hists[[i]]$edge[, 2] %in% fail.internal), 1]

      fail.internal.leading.edge <- c()
      fail.internal.leading.maps <- list()
      for (j in 1:length(fail.internal.start)) {
        if (fail.internal.start[j] != (length(hists[[i]]$tip.label) + 1)) {
          fail.internal.leading.edge <- c(fail.internal.leading.edge,
                                          which(hists[[i]]$edge[, 2] ==
                                                  fail.internal.start[j]))
          fail.internal.leading.maps[j] <- hists[[i]]$maps[which(hists[[i]]$edge[, 2] == fail.internal.start[j])]
        } else {
          # Root edge case: use sibling edge's first state
          fail.internal.leading.edge <- c(fail.internal.leading.edge, 0)
          sibling.leading.edge <- which(hists[[i]]$edge[, 1] == hists[[i]]$edge[fail.internal.edges[j], 1] &
                                          hists[[i]]$edge[, 2] != hists[[i]]$edge[fail.internal.edges[j], 2])
          fail.internal.leading.maps[[j]] <- hists[[i]]$maps[[sibling.leading.edge]][1]
        }
      }

      fail.internal.start.states <- c()
      for (j in 1:length(fail.internal.leading.maps)) {
        fail.internal.start.states[j] <- names(fail.internal.leading.maps[[j]])[
          length(fail.internal.leading.maps[[j]])]
      }
      names(fail.internal.start.states) <- fail.internal.start

      # Remove duplicate entries
      if (TRUE %in% duplicated(names(fail.internal.start.states))) {
        fail.internal.start.states <- fail.internal.start.states[
          -which(duplicated(names(fail.internal.start.states)))]
      }

      fail.internal.end.states <- rep("fail", length(fail.internal))
      names(fail.internal.end.states) <- fail.internal

      fail.internal.numbers <- matrix(hists[[i]]$edge[fail.internal.edges, ],
                                      ncol = 2)
      branches.internal <- cbind(fail.internal.numbers,
                                 hists[[i]]$edge.length[fail.internal.edges])

      if (!is.null(branches) && ncol(branches) == 3) {
        branches <- rbind(branches, branches.internal)
      } else {
        branches <- branches.internal
      }
    }

    # ---- Build the working table for all failed branches ----
    branches <- cbind(branches,
                      rep(NA, nrow(branches)),  # extra edge length
                      rep(NA, nrow(branches)),  # edge id
                      rep(NA, nrow(branches)),  # extra edge id
                      rep(NA, nrow(branches)),  # start state
                      rep(NA, nrow(branches)))  # end state

    colnames(branches) <- c("start", "end", "edge.length",
                             "extra.edge.length", "edge", "extra.edge",
                             "start.state", "end.state")

    # Fill in edge IDs
    if (length(fail.internal.edges) == 0) {
      branches[, 5] <- fail.tips.edges
    } else if (nrow(fail.taxa) == 0) {
      branches[, 5] <- fail.internal.edges
    } else {
      branches[, 5] <- c(fail.tips.edges, fail.internal.edges)
    }

    # ---- Determine start and end states for each failed branch ----
    for (l in 1:nrow(branches)) {

      # -- Start states --
      if (branches[l, 1] %in% branches[, 2]) {
        # Leading edge also failed - use sibling
        sibling.edge <- which(hists[[i]]$edge[, 1] %in% branches[l, 1] &
                                !hists[[i]]$edge[, 2] %in% branches[l, 2])
        sibling.map <- hists[[i]]$maps[sibling.edge]

        if (names(sibling.map[[1]])[1] == "fail") {
          # Sibling also failed - go up another level
          start.1 <- names(fail.internal.start.states[
            names(fail.internal.end.states) == branches[l, 1]])
          extra.length <- hists[[i]]$edge.length[
            which(start.1 == hists[[i]]$edge[, 1]
                  & branches[l, 1] == hists[[i]]$edge[, 2])]
          branches[l, 4] <- extra.length
          branches[l, 6] <- which(start.1 == hists[[i]]$edge[, 1]
                                  & branches[l, 1] == hists[[i]]$edge[, 2])
          branches[l, 7] <- fail.internal.start.states[
            names(fail.internal.end.states) == branches[l, 1]]
        } else {
          branches[l, 7] <- names(sibling.map[[1]])[1]
        }
      } else {
        if (as.numeric(branches[l, 2]) <= length(hists[[i]]$tip.label)) {
          branches[l, 7] <-
            fail.tips.start.states[names(fail.tips.start.states) == branches[l, 1]]
        } else {
          branches[l, 7] <-
            fail.internal.start.states[names(fail.internal.start.states) == branches[l, 1]]
        }
      }

      # -- End states --
      if (as.numeric(branches[l, 2]) <= length(hists[[i]]$tip.label)) {
        branches[l, 8] <-
          fail.tips.end.states[names(fail.tips.end.states) == branches[l, 2]]
      } else {
        if (branches[l, 2] %in% branches[, 1]) {
          descendent.sibling.edge <- which(hists[[i]]$edge[, 1] %in% branches[l, 2] &
                                             !hists[[i]]$edge[, 2] %in% branches[, 2])
          if (length(descendent.sibling.edge) == 0) {
            branches[l, 8] <- NA
          } else {
            descendent.sibling.map <- hists[[i]]$maps[descendent.sibling.edge]
            branches[l, 8] <- names(descendent.sibling.map[[1]])[1]
          }
        } else {
          trailing <- which(hists[[i]]$edge[, 1] %in% branches[l, 2])[1]
          trailing.maps <- hists[[i]]$maps[trailing]
          branches[l, 8] <- names(trailing.maps[[1]][1])
        }
      }
    }

    # Remove branches with NA end states (redundant with other entries)
    if (NA %in% branches[, 8]) {
      branches <- branches[-which(is.na(branches[, 8])), ]
    }
    branches <- as.matrix(branches)

    # ---- Find valid paths and repair maps ----
    # Build a graph from the transition matrix to find shortest paths.
    # Vertices must be named by STATE LABEL (dimnames), not by matrix
    # row/column index, so that states such as "0"/"1" resolve
    allowed <- which(transition.matrix != 0, arr.ind = TRUE)
    if (!is.null(dimnames(transition.matrix))) {
      graph.edges <- data.frame(
        from = rownames(transition.matrix)[allowed[, 1]],
        to = colnames(transition.matrix)[allowed[, 2]])
    } else {
      graph.edges <- allowed
    }
    graph.paths <- graph_from_data_frame(graph.edges)

    for (j in 1:nrow(branches)) {
      # Find the shortest valid path between start and end states
      found.path <- shortest_paths(graph.paths,
                                   branches[j, 7],
                                   branches[j, 8])[[1]][[1]]

      if (is.na(branches[j, 6])) {
        # Simple case: single failed edge
        fail.time <- as.numeric(branches[j, 3]) / length(found.path)
        hists[[i]]$maps[[as.numeric(branches[j, 5])]] <-
          rep(fail.time, length(found.path))
        names(hists[[i]]$maps[[as.numeric(branches[j, 5])]]) <-
          names(found.path)
      } else {
        # Complex case: failed edge connected to another failed edge
        fail.time <- (as.numeric(branches[j, 3]) +
                        as.numeric(branches[j, 4])) / length(found.path)

        leading.multiple <- as.integer(as.numeric(branches[j, 4]) / fail.time)
        trailing.multiple <- as.integer(as.numeric(branches[j, 3]) / fail.time)
        leading.leftover <- as.numeric(branches[j, 4]) - (leading.multiple * fail.time)
        trailing.leftover <- as.numeric(branches[j, 3]) - (trailing.multiple * fail.time)

        # Fix the leading edge map
        hists[[i]]$maps[[as.numeric(branches[j, 6])]] <-
          c(rep(fail.time, leading.multiple), leading.leftover)
        names(hists[[i]]$maps[[as.numeric(branches[j, 6])]]) <-
          names(found.path)[1:(leading.multiple + 1)]

        # Fix the trailing edge map
        hists[[i]]$maps[[as.numeric(branches[j, 5])]] <-
          c(trailing.leftover, rep(fail.time, trailing.multiple))
        names(hists[[i]]$maps[[as.numeric(branches[j, 5])]]) <-
          names(found.path)[(length(found.path) -
                               trailing.multiple):length(found.path)]

        # Handle sibling failures sharing the same leading edge
        # (>= 2 because branches[j, ] always matches itself)
        if (sum(as.numeric(branches[, 6]) == branches[j, 6], na.rm = TRUE) >= 2) {
          repeat.index <- which(branches[, 6] == branches[j, 6] &
                                  branches[, 5] != branches[j, 5])
          branches[repeat.index, 4] <- NA
          branches[repeat.index, 6] <- NA
          branches[repeat.index, 7] <- names(hists[[i]]$maps[[as.numeric(
            branches[j, 5])]])[1]
        }
      }
    }

    # ---- Update mapped.edge matrix to match repaired maps ----
    for (j in 1:nrow(branches)) {
      edge <- as.numeric(branches[j, 5])
      for (k in 1:length(hists[[i]]$maps[[edge]])) {
        state <- names(hists[[i]]$maps[[edge]])[k]
        if (state %in% colnames(hists[[i]]$mapped.edge)) {
          hists[[i]]$mapped.edge[edge, state] <- hists[[i]]$maps[[edge]][k]
        } else {
          # Add a new column for a state not previously seen
          colnames.new <- c(colnames(hists[[i]]$mapped.edge), state)
          hists[[i]]$mapped.edge <- cbind(hists[[i]]$mapped.edge, c(0))
          colnames(hists[[i]]$mapped.edge) <- colnames.new
          hists[[i]]$mapped.edge <- hists[[i]]$mapped.edge[,
            sort(colnames(hists[[i]]$mapped.edge))]
          hists[[i]]$mapped.edge[edge, state] <- hists[[i]]$maps[[edge]][k]
        }
      }
      # Zero out the "fail" column for this edge
      if ("fail" %in% colnames(hists[[i]]$mapped.edge)) {
        hists[[i]]$mapped.edge[edge, "fail"] <- 0
      }
    }
  }

  # ---- Return the fixed object ----
  # Unwrap from list if there was only one simulation
  if (length(hists) == 1) {
    hists <- hists[[1]]
  }
  return(hists)
}
