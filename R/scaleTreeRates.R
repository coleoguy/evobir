#' Analyze Rate Heterogeneity in Discrete Character Evolution
#'
#' Tests whether different branches of a phylogeny have different rates
#' of discrete character evolution.  The function assigns a "scalar"
#' multiplier to each edge: scalar > 1 means faster-than-average
#' evolution, scalar < 1 means slower, and scalar = 1 means the
#' background rate.
#'
#' The method works by trying different scalar values for each edge
#' (via dynamic programming) and keeping whichever value maximizes
#' the likelihood.
#'
#' @param tree A phylogenetic tree of class \code{"phylo"}.
#' @param tip.states A named vector of discrete character states at
#'   the tips.  Names must match \code{tree$tip.label}.
#' @param model Character or matrix.  The model for
#'   \code{\link[phytools]{fitMk}}:
#'   \itemize{
#'     \item \code{"ER"} -- equal rates
#'     \item \code{"SYM"} -- symmetric rates
#'     \item \code{"ARD"} -- all rates different
#'     \item A custom rate index matrix
#'   }
#' @param fixedQ A pre-estimated Q (transition rate) matrix.
#'   If \code{NULL} (default), rates are estimated from the data.
#' @param max.ratio Numeric.  Maximum rate scalar to test (default 2).
#'   The function tests scalars from \code{1/max.ratio} to
#'   \code{max.ratio}.
#' @param nbins Integer.  Number of scalar bins above and below 1
#'   (default 10).  More bins gives finer resolution but takes longer.
#' @param max.transition Integer.  Maximum number of bin steps an edge's
#'   scalar can differ from its parent's scalar (default 1).  This
#'   prevents rate estimates from jumping wildly between adjacent
#'   branches.
#' @param var.start Logical.  If \code{TRUE}, tries multiple starting
#'   scalar values at the root and returns the best tree.  If
#'   \code{FALSE} (default), fixes the root scalar at 1.
#' @param pi Character.  Method for estimating the root state prior.
#'   Passed to \code{\link[phytools]{fitMk}}.  Default is
#'   \code{"fitzjohn"}.
#'
#' @return A tree of class \code{c("phyloscaled", "phylo")} with an
#'   additional \code{$scalar} element: a numeric vector (one value
#'   per edge) of rate multipliers.  Plot with
#'   \code{\link{plot.phyloscaled}}.
#'
#' @details
#' The algorithm proceeds edge-by-edge from root to tips.  For each
#' edge, it tries every allowed scalar value (constrained to be within
#' \code{max.transition} bins of the parent's scalar), fits the Mk
#' model with scaled branch lengths, and picks the scalar that gives
#' the highest likelihood.
#'
#' @examples
#' \dontrun{
#' library(phytools)
#' tree <- rcoal(20)
#' states <- setNames(sample(1:3, 20, replace = TRUE), tree$tip.label)
#' scaled <- scaleTreeRates(tree, states, model = "ARD", nbins = 5)
#' plot(scaled, method = "color")
#' }
#'
#' @importFrom phytools fitMk nodeHeights
#' @importFrom expm expm
#' @importFrom stats optim nlminb reorder optimize
#' @export
scaleTreeRates <- function(tree, tip.states, model, fixedQ = NULL,
                           max.ratio = 2, nbins = 10,
                           max.transition = 1, var.start = FALSE,
                           pi = "fitzjohn") {

  # ---- Internal helper: custom fitMk for likelihood calculations ----
  # This is a local version of phytools::fitMk needed for the
  # edge-by-edge likelihood computations
  fitMkNew <- function(tree, x, model = "SYM", fixedQ = NULL, ...) {
    if (hasArg(output.liks))
      output.liks <- list(...)$output.liks
    else output.liks <- FALSE
    if (hasArg(q.init))
      q.init <- list(...)$q.init
    else q.init <- length(unique(x)) / sum(tree$edge.length)
    if (hasArg(opt.method))
      opt.method <- list(...)$opt.method
    else opt.method <- "nlminb"
    if (hasArg(min.q))
      min.q <- list(...)$min.q
    else min.q <- 1e-12
    if (hasArg(max.q))
      max.q <- list(...)$max.q
    else max.q <- max(nodeHeights(tree)) * 100
    if (hasArg(logscale))
      logscale <- list(...)$logscale
    else logscale <- FALSE
    N <- Ntip(tree)
    M <- tree$Nnode
    if (is.matrix(x)) {
      x <- x[tree$tip.label, ]
      m <- ncol(x)
      states <- colnames(x)
    } else {
      x <- to.matrix(x, sort(unique(x)))
      m <- ncol(x)
      states <- colnames(x)
    }
    if (hasArg(pi))
      pi <- list(...)$pi
    else pi <- "equal"
    if (is.numeric(pi))
      root.prior <- "given"
    if (pi[1] == "equal") {
      pi <- setNames(rep(1/m, m), states)
      root.prior <- "flat"
    } else if (pi[1] == "estimated") {
      pi <- if (!is.null(fixedQ))
        statdist(fixedQ)
      else statdist(summary(fitMk(tree, x, model), quiet = TRUE)$Q)
      cat(paste("Using pi estimated from the stationary",
                "distribution of Q assuming a flat prior.\npi =\n"))
      print(round(pi, 6))
      cat("\n")
      root.prior <- "stationary"
    } else if (pi[1] == "fitzjohn")
      root.prior <- "nuisance"
    if (is.numeric(pi)) {
      pi <- pi / sum(pi)
      if (is.null(names(pi)))
        pi <- setNames(pi, states)
      pi <- pi[states]
    }
    if (is.null(fixedQ)) {
      if (is.character(model)) {
        rate <- matrix(NA, m, m)
        if (model == "ER") {
          k <- rate[] <- 1
          diag(rate) <- NA
        } else if (model == "ARD") {
          k <- m * (m - 1)
          rate[col(rate) != row(rate)] <- 1:k
        } else if (model == "SYM") {
          k <- m * (m - 1) / 2
          ii <- col(rate) < row(rate)
          rate[ii] <- 1:k
          rate <- t(rate)
          rate[ii] <- 1:k
        }
      } else {
        if (ncol(model) != nrow(model))
          stop("model is not a square matrix")
        rate <- model
        m <- ncol(rate)
        states <- as.character(1:ncol(rate))
        k <- max(rate)
        if (length(states) != ncol(x)) {
          missing <- states[which(!states %in% colnames(x))]
          x <- cbind(x, matrix(data = 0, nrow = nrow(x), ncol = length(missing),
                                dimnames = list(rownames(x), missing)))
          x <- x[, states]
        }
      }
      Q <- matrix(0, m, m)
    } else {
      m <- ncol(fixedQ)
      states <- as.character(1:ncol(fixedQ))
      rate <- matrix(NA, m, m)
      k <- m * (m - 1)
      rate[col(rate) != row(rate)] <- 1:k
      Q <- fixedQ
      if (ncol(fixedQ) != ncol(x)) {
        missing <- states[which(!states %in% colnames(x))]
        x <- cbind(x, matrix(data = 0, nrow = nrow(x), ncol = length(missing),
                              dimnames = list(rownames(x), missing)))
      }
    }
    index.matrix <- rate
    tmp <- cbind(1:m, 1:m)
    rate[tmp] <- 0
    rate[rate == 0] <- k + 1
    liks <- rbind(x, matrix(0, M, m, dimnames = list(1:M + N, states)))
    pw <- reorder(tree, "pruningwise")

    # Likelihood function (pruning algorithm)
    lik <- function(Q, output.liks = FALSE, pi, ...) {
      if (hasArg(output.pi))
        output.pi <- list(...)$output.pi
      else output.pi <- FALSE
      if (is.Qmatrix(Q)) Q <- unclass(Q)
      if (any(is.nan(Q)) || any(is.infinite(Q))) return(1e+50)
      comp <- vector(length = N + M, mode = "numeric")
      parents <- unique(pw$edge[, 1])
      root <- min(parents)
      for (i in 1:length(parents)) {
        anc <- parents[i]
        ii <- which(pw$edge[, 1] == parents[i])
        desc <- pw$edge[ii, 2]
        el <- pw$edge.length[ii]
        v <- vector(length = length(desc), mode = "list")
        for (j in 1:length(v)) {
          v[[j]] <- EXPM(Q * el[j]) %*% liks[desc[j], ]
        }
        if (anc == root) {
          if (is.numeric(pi))
            vv <- Reduce("*", v)[, 1] * pi
          else if (pi[1] == "fitzjohn") {
            D <- Reduce("*", v)[, 1]
            pi <- D / sum(D)
            vv <- D * D / sum(D)
          }
        } else {
          vv <- Reduce("*", v)[, 1]
        }
        comp[anc] <- sum(vv)
        liks[anc, ] <- vv / comp[anc]
      }
      if (output.liks) return(liks[1:M + N, , drop = FALSE])
      else if (output.pi) return(pi)
      else {
        logL <- -sum(log(comp[1:M + N]))
        if (is.na(logL)) logL <- Inf
        return(logL)
      }
    }

    if (is.null(fixedQ)) {
      if (length(q.init) != k)
        q.init <- rep(q.init[1], k)
      q.init <- if (logscale) log(q.init) else q.init
      if (opt.method == "optim") {
        fit <- if (logscale)
          optim(q.init, function(p) lik(makeQ(m, exp(p), index.matrix), pi = pi),
                method = "L-BFGS-B",
                lower = rep(log(min.q), k), upper = rep(log(max.q), k))
        else
          optim(q.init, function(p) lik(makeQ(m, p, index.matrix), pi = pi),
                method = "L-BFGS-B",
                lower = rep(min.q, k), upper = rep(max.q, k))
      } else if (opt.method == "none") {
        fit <- list(objective = lik(makeQ(m, q.init, index.matrix), pi = pi),
                    par = q.init)
      } else {
        fit <- if (logscale)
          nlminb(q.init, function(p) lik(makeQ(m, exp(p), index.matrix), pi = pi),
                 lower = rep(log(min.q), k), upper = rep(log(max.q), k))
        else
          nlminb(q.init, function(p) lik(makeQ(m, p, index.matrix), pi = pi),
                 lower = rep(0, k), upper = rep(max.q, k))
      }
      if (logscale) fit$par <- exp(fit$par)
      if (pi[1] == "fitzjohn")
        pi <- setNames(lik(makeQ(m, fit$par, index.matrix),
                           FALSE, pi = pi, output.pi = TRUE), states)
      obj <- list(
        logLik = if (opt.method == "optim") -fit$value else -fit$objective,
        rates = fit$par, index.matrix = index.matrix, states = states,
        pi = pi, method = opt.method, root.prior = root.prior)
      if (output.liks)
        obj$lik.anc <- lik(makeQ(m, obj$rates, index.matrix), TRUE, pi = pi)
    } else {
      fit <- lik(Q, pi = pi)
      if (pi[1] == "fitzjohn")
        pi <- setNames(lik(Q, FALSE, pi = pi, output.pi = TRUE), states)
      obj <- list(
        logLik = -fit,
        rates = Q[sapply(1:k, function(x, y) which(x == y), index.matrix)],
        index.matrix = index.matrix, states = states,
        pi = pi, root.prior = root.prior)
      if (output.liks)
        obj$lik.anc <- lik(makeQ(m, obj$rates, index.matrix), TRUE, pi = pi)
    }
    lik.f <- function(q) -lik(q, output.liks = FALSE,
                               pi = if (root.prior == "nuisance") "fitzjohn" else pi)
    obj$lik <- lik.f
    class(obj) <- "fitMk"
    return(obj)
  }

  # ---- Small helper functions ----
  is.Qmatrix <- function(x) "Qmatrix" %in% class(x)

  makeQ <- function(m, q, index.matrix) {
    Q <- matrix(0, m, m)
    Q[] <- c(0, q)[index.matrix + 1]
    diag(Q) <- 0
    diag(Q) <- -rowSums(Q)
    Q
  }

  # Memoization cache for matrix exponentials.  During the main
  # edge-scaling loop Q is fixed, so the product Q * t (and hence the
  # cache key) is determined by the scaled edge length; the vast
  # majority of EXPM calls are repeats of previously seen edge lengths.
  # The cache is only enabled once Q is fixed (see below) so it is
  # never polluted during the initial rate estimation, where Q varies.
  expm.cache <- new.env(hash = TRUE, parent = emptyenv())
  use.expm.cache <- FALSE

  EXPM <- function(x, ...) {
    if (use.expm.cache) {
      key <- paste(x, collapse = " ")
      e_x <- expm.cache[[key]]
      if (is.null(e_x)) {
        e_x <- if (isSymmetric(x)) matexpo(x) else expm(x, ...)
        expm.cache[[key]] <- e_x
      }
      dimnames(e_x) <- dimnames(x)
      return(e_x)
    }
    e_x <- if (isSymmetric(x)) matexpo(x) else expm(x, ...)
    dimnames(e_x) <- dimnames(x)
    e_x
  }

  statdist <- function(Q) {
    foo <- function(theta, Q) {
      Pi <- c(theta[1:(nrow(Q) - 1)], 1 - sum(theta[1:(nrow(Q) - 1)]))
      sum((Pi %*% Q)^2)
    }
    k <- nrow(Q)
    if (nrow(Q) > 2) {
      fit <- optim(rep(1/k, k - 1), foo, Q = Q,
                   control = list(reltol = 1e-16))
      return(setNames(c(fit$par[1:(k - 1)], 1 - sum(fit$par[1:(k - 1)])),
                      rownames(Q)))
    } else {
      fit <- optimize(foo, interval = c(0, 1), Q = Q)
      return(setNames(c(fit$minimum, 1 - fit$minimum), rownames(Q)))
    }
  }

  # ---- Sort tip states to match tree order ----
  x <- tip.states[order(factor(names(tip.states), levels = tree$tip.label))]

  # ---- Estimate or use supplied Q matrix ----
  if (!is.null(fixedQ)) {
    print("using supplied Q-matrix")
    Q <- fixedQ
  } else {
    print("estimating rates")
    fit <- fitMkNew(tree, x, model = model, pi = pi)
    rates <- fit$rates
    print("building Q-matrix")
    Q <- fit$index.matrix
    Q[is.na(Q)] <- 0
    for (i in 1:length(rates)) {
      Q[Q == i] <- rates[i]
    }
    diag(Q) <- -rowSums(Q)
  }

  # Q is now fixed for the remainder of the search: turn on the
  # matrix-exponential cache (keyed on Q * scaled edge length)
  use.expm.cache <- TRUE

  # ---- Define the scalar bins to test ----
  # Creates a symmetric set of bins around 1
  set.range <- max.transition
  bins <- c(1 / seq(from = max.ratio, to = 1, length.out = nbins + 1),
            seq(from = 1, to = max.ratio, length.out = nbins + 1))
  bins <- unique(bins)

  # Initialize all edges with scalar = 1 (background rate)
  tree$scalar <- rep(1, times = length(tree$edge.length))
  trees <- list()

  # Determine starting scalars to try
  if (var.start == FALSE) {
    start.bins <- 1
  } else {
    start.bins <- bins
  }

  # ---- Main optimization loop ----
  # For each starting scalar, traverse edges and find best scalar per edge
  for (i in 1:length(start.bins)) {
    print(paste0("starting scalar: ", start.bins[i]))
    start.scalar <- start.bins[i]

    for (j in 1:length(tree$edge.length)) {
      print(paste0("starting scalar ", start.bins[i], ", edge ", j))

      # Get parent edge's scalar (or use starting scalar for root edges)
      if (tree$edge[j, 1] == length(tree$tip.label) + 1) {
        prev.scalar <- start.scalar
      } else {
        prev.edge <- which(tree$edge[, 2] == tree$edge[j, 1])
        prev.scalar <- tree$scalar[prev.edge]
      }

      # Determine which scalar bins are reachable from the parent
      parent.bin <- which(bins == prev.scalar)
      low.bound <- max(1, parent.bin - set.range)
      high.bound <- min(length(bins), parent.bin + set.range)
      poss.scalars <- bins[low.bound:high.bound]

      # Test each possible scalar and pick the best likelihood
      liks <- c()
      for (k in 1:length(poss.scalars)) {
        temp.tree <- tree
        temp.tree$scalar[j] <- poss.scalars[k]
        temp.tree$edge.length <- temp.tree$edge.length * temp.tree$scalar
        liks[k] <- fitMkNew(temp.tree, x, model = model,
                            pi = pi, fixedQ = Q)$logLik
      }

      # Select the scalar with highest likelihood
      index.best <- which(liks == max(liks))
      # If tied, keep the parent's scalar (parsimony)
      if (length(index.best) > 1) {
        index.best <- which(poss.scalars == prev.scalar)
      }
      tree$scalar[j] <- poss.scalars[index.best]
    }

    trees[[i]] <- tree
  }

  # ---- Select the best starting scalar (if var.start) ----
  if (var.start == FALSE) {
    final.tree <- trees[[1]]
  } else {
    final.liks <- c()
    for (i in 1:length(trees)) {
      print(paste0("getting final likelihood for starting scalar ", start.bins[i]))
      temp.tree <- trees[[i]]
      temp.tree$edge.length <- temp.tree$edge.length * temp.tree$scalar
      final.liks[i] <- fitMkNew(temp.tree, x, model = model,
                                fixedQ = Q)$logLik
    }
    max.final.lik <- max(final.liks)
    max.index <- which.max(final.liks)
    final.start <- bins[max.index]
    final.tree <- trees[[max.index]]
    print(paste0("the best starting scalar is ", final.start,
                 ", associated logLik is ", max.final.lik))
  }

  # "phyloscaled" must come FIRST so that plot() dispatches to
  # plot.phyloscaled rather than ape::plot.phylo
  class(final.tree) <- c("phyloscaled", "phylo")
  print(paste0("returning tree with scalar"))
  return(final.tree)
}
