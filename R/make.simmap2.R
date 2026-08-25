#' Enhanced Stochastic Character Mapping
#'
#' An enhanced version of \code{\link[phytools]{make.simmap}} that adds
#' configurable rejection limits and monitoring for large or complex
#' datasets.  Stochastic character mapping simulates possible histories
#' of a discrete trait on a phylogeny -- imagine painting each branch of
#' the tree with colors representing different trait states.
#'
#' The key enhancement is the \code{rejmax} parameter: when the Gillespie
#' simulation can't find a valid path on a branch after \code{rejmax}
#' attempts, it marks that branch as "fail" instead of running forever.
#' Failed branches can then be repaired with \code{\link{fix.simmap}}.
#'
#' @param tree A phylogenetic tree of class \code{"phylo"} or
#'   \code{"multiPhylo"}.
#' @param x A named vector of discrete character states at the tips,
#'   or a matrix of state probabilities (rows = tips, columns = states).
#' @param model Character or matrix.  Model of trait evolution:
#'   \itemize{
#'     \item \code{"SYM"} (default) -- symmetric rates
#'     \item \code{"ER"} -- equal rates
#'     \item \code{"ARD"} -- all rates different
#'     \item A custom rate matrix
#'   }
#' @param nsim Integer.  Number of stochastic maps to simulate
#'   (default 1).
#' @param monitor Logical.  If \code{TRUE}, prints progress updates
#'   showing acceptance/rejection of each branch.  Useful for
#'   debugging slow runs.  Default is \code{FALSE}.
#' @param rejmax Integer or \code{NULL}.  Maximum number of rejection
#'   sampling attempts per branch before marking it as "fail."
#'   \code{NULL} (default) means no limit (keeps trying forever).
#' @param rejint Integer.  Interval for printing rejection count
#'   progress messages (default 1000000).  Only used when
#'   \code{monitor = TRUE}.
#' @param ... Additional arguments passed to \code{\link[phytools]{fitMk}}
#'   for rate estimation.  Common options include:
#'   \describe{
#'     \item{Q}{Transition rate matrix: \code{"empirical"} (default,
#'       estimated from the data), \code{"mcmc"} (sampled from posterior),
#'       or a fixed matrix.}
#'     \item{pi}{Root state priors: \code{"equal"}, \code{"estimated"},
#'       \code{"fitzjohn"}, or a numeric vector.}
#'     \item{message}{Logical.  Print the Q matrix? Default \code{TRUE}.}
#'   }
#'
#' @return A tree of class \code{"simmap"} (if \code{nsim = 1}) or a
#'   list of class \code{c("multiSimmap", "multiPhylo")} (if
#'   \code{nsim > 1}).  Each tree has additional components:
#'   \describe{
#'     \item{maps}{A list of named vectors, one per edge, giving the
#'       time spent in each state along that edge.}
#'     \item{mapped.edge}{A matrix summarizing total time in each state
#'       per edge.}
#'     \item{Q}{The transition rate matrix used.}
#'     \item{logL}{The log-likelihood.}
#'     \item{fail}{Integer vector of edge indices that failed
#'       (hit \code{rejmax}).}
#'   }
#'
#' @seealso \code{\link{fix.simmap}} to repair failed branches.
#'
#' @examples
#' \dontrun{
#' library(phytools)
#' # Basic usage (same as make.simmap)
#' mapped <- make.simmap2(tree, states, model = "ARD", nsim = 100)
#'
#' # With rejection limit and monitoring
#' mapped <- make.simmap2(tree, states, model = "ARD", nsim = 10,
#'                         rejmax = 500000, monitor = TRUE)
#' # Fix any failed branches
#' fixed <- fix.simmap(mapped, tip_data, trans_matrix)
#' }
#'
#' @importFrom phytools fitMk matchNodes to.matrix rstate
#' @importFrom ape Ntip matexpo is.binary reorder.phylo multi2di
#' @importFrom expm expm
#' @importFrom methods hasArg
#' @importFrom stats logLik rnorm rgamma dgamma runif rexp
#' @importFrom utils flush.console
#' @export
make.simmap2 <- function(tree, x, model = "SYM", nsim = 1,
                          monitor = FALSE, rejmax = NULL,
                          rejint = 1000000, ...) {

  # ---- Handle multiPhylo input ----
  # If given multiple trees, apply make.simmap2 to each one so that
  # rejmax/rejint/monitor are honored for every tree
  if (inherits(tree, "multiPhylo")) {
    ff <- function(yy, x, model, nsim, ...) {
      zz <- make.simmap2(yy, x, model, nsim, monitor = monitor,
                         rejmax = rejmax, rejint = rejint, ...)
      if (nsim > 1) class(zz) <- NULL
      return(zz)
    }
    if (nsim > 1)
      mtrees <- unlist(sapply(tree, ff, x, model, nsim,
                              ..., simplify = FALSE), recursive = FALSE)
    else
      mtrees <- sapply(tree, ff, x, model, nsim, ...,
                        simplify = FALSE)
    class(mtrees) <- c("multiSimmap", "multiPhylo")
  } else {

    # ---- Parse additional arguments ----
    if (hasArg(pi)) pi <- list(...)$pi else pi <- "equal"
    if (hasArg(message)) pm <- list(...)$message else pm <- TRUE
    if (hasArg(tol)) tol <- list(...)$tol else tol <- 0
    if (hasArg(Q)) Q <- list(...)$Q else Q <- "empirical"
    if (hasArg(burnin)) burnin <- list(...)$burnin else burnin <- 1000
    if (hasArg(samplefreq)) samplefreq <- list(...)$samplefreq else samplefreq <- 100
    if (hasArg(vQ)) vQ <- list(...)$vQ else vQ <- 0.1
    prior <- list(alpha = 1, beta = 1, use.empirical = FALSE)
    if (hasArg(prior)) {
      pr <- list(...)$prior
      prior[names(pr)] <- pr
    }

    # Validate the tree
    if (!inherits(tree, "phylo"))
      stop("tree should be object of class \"phylo\".")

    # Convert tip states to a matrix format
    if (!is.matrix(x))
      xx <- to.matrix(x, sort(unique(x)))
    else
      xx <- x
    xx <- xx[tree$tip.label, ]
    xx <- xx / rowSums(xx)

    # Reorder tree for efficient traversal
    tree <- bt <- reorder.phylo(tree, "cladewise")
    if (!is.binary(bt)) bt <- multi2di(bt, random = FALSE)
    N <- Ntip(tree)
    m <- ncol(xx)
    root <- N + 1
    sim <- 0

    # ---- Method 1: Empirical Q matrix ----
    # Estimate transition rates from the data using ML
    if (is.character(Q) && Q == "empirical") {
      XX <- getPars(bt, xx, model, Q = NULL, tree, tol,
                    m, pi = pi, args = list(...))
      L <- XX$L
      Q <- XX$Q
      logL <- XX$loglik
      pi <- XX$pi
      if (pi[1] == "equal")
        pi <- setNames(rep(1/m, m), colnames(L))
      else if (pi[1] == "estimated")
        pi <- statdist(Q)
      else if (pi[1] == "fitzjohn")
        pi <- "fitzjohn"
      else
        pi <- pi / sum(pi)
      if (pm) printmessage(Q, pi, method = "empirical")

      # Simulate character histories using the estimated Q
      mtrees <- replicate(nsim,
                          smap2(tree, x, N, m, root, L, Q, pi, logL,
                                rejmax, rejint, monitor, sim),
                          simplify = FALSE)
    }

    # ---- Method 2: MCMC sampling of Q ----
    # Sample Q matrices from the posterior distribution
    else if (is.character(Q) && Q == "mcmc") {
      if (prior$use.empirical) {
        qq <- fitMk(bt, xx, model)$rates
        prior$alpha <- qq * prior$beta
      }
      get.stationary <- if (pi[1] == "estimated") TRUE else FALSE
      if (pi[1] %in% c("equal", "estimated"))
        pi <- setNames(rep(1/m, m), colnames(xx))
      else if (pi[1] == "fitzjohn")
        pi <- "fitzjohn"
      else
        pi <- pi / sum(pi)

      XX <- mcmcQ(bt, xx, model, tree, tol, m, burnin,
                  samplefreq, nsim, vQ, prior, pi = pi)
      L <- lapply(XX, function(x) x$L)
      Q <- lapply(XX, function(x) x$Q)
      logL <- lapply(XX, function(x) x$loglik)
      pi <- if (get.stationary)
        lapply(Q, statdist)
      else if (pi[1] == "fitzjohn")
        lapply(XX, function(x) x$pi)
      else
        lapply(1:nsim, function(x, y) y, y = pi)
      if (pm)
        printmessage(Reduce("+", Q) / length(Q),
                     Reduce("+", pi) / length(pi), method = "mcmc")

      mtrees <- if (nsim > 1)
        mapply(smap2, L = L, Q = Q, pi = pi, logL = logL,
               MoreArgs = list(tree = tree, x = x, N = N,
                               m = m, root = root, rejmax = rejmax,
                               rejint = rejint, monitor = monitor,
                               sim = sim), SIMPLIFY = FALSE)
      else
        list(smap2(tree = tree, x = x, N = N, m = m,
                    root = root, rejmax = rejmax, rejint = rejint,
                    monitor = monitor, sim = sim,
                    L = L[[1]], Q = Q[[1]], pi = pi[[1]],
                    logL = logL[[1]]))
    }

    # ---- Method 3: Fixed Q matrix ----
    # User provides their own transition rate matrix
    else if (is.matrix(Q)) {
      XX <- getPars(bt, xx, model, Q = Q, tree, tol, m,
                    pi = pi, args = list(...))
      # state-name indexing downstream requires dimnames on Q; if the
      # user's matrix has none, assume states in the same order as the
      # (sorted) data columns
      if (is.null(dimnames(Q)))
        dimnames(Q) <- list(colnames(XX$L), colnames(XX$L))
      L <- XX$L
      logL <- XX$loglik
      pi <- XX$pi
      if (pi[1] == "equal")
        pi <- setNames(rep(1/m, m), colnames(L))
      else if (pi[1] == "estimated")
        pi <- statdist(Q)
      else if (pi[1] == "fitzjohn")
        pi <- "fitzjohn"
      else
        pi <- pi / sum(pi)
      if (pm) printmessage(Q, pi, method = "fixed")

      # Use a for loop (not replicate) so we can monitor progress
      mtrees <- vector("list", nsim)
      for (sim in 1:nsim) {
        if (monitor == TRUE) {
          print(paste("simulation", sim, "of", nsim, sep = " "))
        }
        mtrees[[sim]] <- smap2(tree, x, N, m, root, L, Q, pi, logL,
                               rejmax, rejint, monitor, sim)
      }
    }

    # ---- Set proper S3 classes on output ----
    if (nsim == 1) {
      mtrees <- mtrees[[1]]
      class(mtrees) <- c("simmap", "phylo")
    } else {
      class(mtrees) <- c("multiSimmap", "multiPhylo")
      for (i in 1:nsim) {
        class(mtrees[[i]]) <- c("simmap", "phylo")
      }
    }
  }

  # Print completion message
  (if (hasArg(message)) list(...)$message else TRUE)
  if ((if (hasArg(message)) list(...)$message else TRUE) &&
      inherits(tree, "phylo"))
    message("Done.")

  return(mtrees)
}


# ---- Core simulation function (internal) ----
# Simulates a single stochastic character map on a tree using the
# Gillespie algorithm with rejection sampling.
#
# @keywords internal
smap2 <- function(tree, x, N, m, root, L, Q, pi, logL,
                  rejmax, rejint, monitor, sim) {

  # Create the output tree with empty maps

  mtree <- tree
  mtree$maps <- list()
  mtree$mapped.edge <- matrix(0, nrow(tree$edge), m,
                              dimnames = list(
                                paste(tree$edge[, 1], ",",
                                      tree$edge[, 2], sep = ""),
                                colnames(L)))

  # ---- Assign states at nodes via pre-order traversal ----
  # Start at the root and work down, drawing states from the
  # conditional probabilities at each node
  NN <- matrix(NA, nrow(tree$edge), 2)
  NN[which(tree$edge[, 1] == root), 1] <-
    rstate(L[as.character(root), ] / sum(L[as.character(root), ]))

  # Track failed branches
  fail <- list()
  fail.count <- 0

  for (j in 1:nrow(tree$edge)) {
    # Draw end state conditioned on start state and data at descendant
    # NOTE: index by state NAME (character), never by as.numeric(state
    # label) -- labels like "0"/"1" or chromosome numbers are not row
    # numbers of Q
    p <- EXPM(Q * tree$edge.length[j])[NN[j, 1], ] *
      L[as.character(tree$edge[j, 2]), ]
    NN[j, 2] <- rstate(p / sum(p))
    NN[which(tree$edge[, 1] == tree$edge[j, 2]), 1] <- NN[j, 2]

    # ---- Simulate the path along this branch ----
    # Use rejection sampling: keep trying until the simulated path
    # ends in the correct state, or we hit rejmax
    accept <- FALSE
    counter <- 0
    while (!accept) {
      map <- sch(NN[j, 1], tree$edge.length[j], Q)
      if (!is.null(rejmax) && counter == rejmax) {
        # Give up on this branch -- mark as failed
        if (monitor == TRUE) {
          print(paste("sim", sim, ": branch", j, "of", nrow(tree$edge),
                      "has exceeded the rejection limit of",
                      format(rejmax, scientific = FALSE),
                      "and will be skipped", sep = " "))
        }
        map <- NA
        fail.count <- fail.count + 1
        fail[[fail.count]] <- j
        accept <- TRUE
      } else if (names(map)[length(map)] == NN[j, 2]) {
        # Path ends in the right state -- accept it
        if (monitor == TRUE) {
          print(paste("sim", sim, ": branch", j, "of", nrow(tree$edge),
                      "ACCEPTED after", format(counter, scientific = FALSE),
                      "total rejections", sep = " "))
        }
        accept <- TRUE
      } else {
        # Wrong ending state -- reject and try again
        counter <- counter + 1
        if ((counter / rejint) %% 1 == 0) {
          if (monitor == TRUE) {
            print(paste("sim", sim, ": branch", j, "of", nrow(tree$edge),
                        "rejected with", format(counter, scientific = FALSE),
                        "total rejections", sep = " "))
          }
        }
      }
    }

    # Store the simulated map
    mtree$maps[[j]] <- map
    for (k in 1:length(mtree$maps[[j]])) {
      mtree$mapped.edge[j, names(mtree$maps[[j]])[k]] <-
        mtree$mapped.edge[j, names(mtree$maps[[j]])[k]] + mtree$maps[[j]][k]
    }
  }

  # Label failed branches
  for (L in 1:length(mtree$maps)) {
    if (is.na(mtree$maps[[L]][1]) == TRUE) {
      named.edge <- mtree$edge.length[L]
      names(named.edge) <- "fail"
      mtree$maps[[L]] <- named.edge
    }
  }

  mtree$Q <- Q
  mtree$logL <- logL
  mtree$fail <- unlist(fail)
  if (!inherits(mtree, "simmap"))
    class(mtree) <- c("simmap", setdiff(class(mtree), "simmap"))
  attr(mtree, "map.order") <- "right-to-left"
  return(mtree)
}


# ---- MCMC Q-matrix sampler (internal) ----
# Estimates the transition rate matrix Q using MCMC with a gamma prior.
#
# @keywords internal
mcmcQ <- function(bt, xx, model, tree, tol, m, burnin,
                  samplefreq, nsim, vQ, prior, pi, args = list()) {
  # Proposal function: perturb rates by small random amounts
  update <- function(x) {
    x <- abs(x + rnorm(n = np, mean = 0, sd = sqrt(vQ)))
    return(x)
  }

  # Build the model matrix (determines which rates are free parameters)
  if (is.character(model)) {
    rate <- matrix(NA, m, m)
    if (model == "ER") {
      np <- rate[] <- 1
      diag(rate) <- NA
    }
    if (model == "ARD") {
      np <- m * (m - 1)
      rate[col(rate) != row(rate)] <- 1:np
    }
    if (model == "SYM") {
      np <- m * (m - 1) / 2
      sel <- col(rate) < row(rate)
      rate[sel] <- 1:np
      rate <- t(rate)
      rate[sel] <- 1:np
    }
  } else {
    if (ncol(model) != nrow(model))
      stop("the matrix given as 'model' is not square")
    if (ncol(model) != m)
      stop("the matrix 'model' must have as many rows as the number of categories in 'x'")
    rate <- model
    np <- max(rate)
  }

  # ---- Burn-in phase ----
  p <- rgamma(np, prior$alpha, prior$beta)
  Q <- matrix(c(0, p)[rate + 1], m, m)
  diag(Q) <- -rowSums(Q, na.rm = TRUE)
  yy <- getPars(bt, xx, model, Q, tree, tol, m, pi = pi, args = args)
  cat("Running MCMC burn-in. Please wait....\n")
  flush.console()

  for (i in 1:burnin) {
    pp <- update(p)
    Qp <- matrix(c(0, pp)[rate + 1], m, m)
    diag(Qp) <- -rowSums(Qp, na.rm = TRUE)
    zz <- getPars(bt, xx, model, Qp, tree, tol, m, FALSE, pi = pi, args)
    # Metropolis-Hastings acceptance ratio
    p.odds <- exp(zz$loglik +
                    sum(dgamma(pp, prior$alpha, prior$beta, log = TRUE)) -
                    yy$loglik -
                    sum(dgamma(p, prior$alpha, prior$beta, log = TRUE)))
    if (p.odds >= runif(n = 1)) {
      yy <- zz
      p <- pp
    }
  }

  # ---- Sampling phase ----
  cat(paste("Running", samplefreq * nsim, "generations of MCMC, sampling every",
            samplefreq, "generations.\nPlease wait....\n\n"))
  flush.console()
  XX <- vector("list", nsim)

  for (i in 1:(samplefreq * nsim)) {
    pp <- update(p)
    Qp <- matrix(c(0, pp)[rate + 1], m, m)
    diag(Qp) <- -rowSums(Qp, na.rm = TRUE)
    zz <- getPars(bt, xx, model, Qp, tree, tol, m, FALSE, pi = pi, args)
    p.odds <- exp(zz$loglik +
                    sum(dgamma(pp, prior$alpha, prior$beta, log = TRUE)) -
                    yy$loglik -
                    sum(dgamma(p, prior$alpha, prior$beta, log = TRUE)))
    if (p.odds >= runif(n = 1)) {
      yy <- zz
      p <- pp
    }
    # Save sampled Q at the specified frequency
    if (i %% samplefreq == 0) {
      Qi <- matrix(c(0, p)[rate + 1], m, m)
      diag(Qi) <- -rowSums(Qi, na.rm = TRUE)
      XX[[i / samplefreq]] <- getPars(bt, xx, model, Qi, tree, tol, m,
                                      TRUE, pi = pi, args)
    }
  }
  return(XX)
}


# ---- Stationary distribution of Q (internal) ----
# Finds the equilibrium frequencies implied by a rate matrix Q.
#
# @keywords internal
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


# ---- Parameter extraction (internal) ----
# Fits the Mk model and extracts Q matrix, likelihoods, and root priors.
#
# @keywords internal
getPars <- function(bt, xx, model, Q, tree, tol, m, liks = TRUE,
                    pi, args = list()) {
  if (!is.null(args$pi)) args$pi <- NULL
  args <- c(list(tree = bt, x = xx, model = model, fixedQ = Q,
                 output.liks = liks, pi = pi), args)
  obj <- do.call(fitMk, args)
  N <- length(bt$tip.label)
  pi <- obj$pi
  II <- obj$index.matrix + 1
  lvls <- obj$states
  if (liks) {
    L <- obj$lik.anc
    rownames(L) <- N + 1:nrow(L)
    if (!is.binary(tree)) {
      ancNames <- matchNodes(tree, bt)
      L <- L[as.character(ancNames[, 2]), ]
      rownames(L) <- ancNames[, 1]
    }
    L <- rbind(xx, L)
    rownames(L)[1:N] <- 1:N
  } else {
    L <- NULL
  }
  Q <- matrix(c(0, obj$rates)[II], m, m, dimnames = list(lvls, lvls))
  if (any(rowSums(Q, na.rm = TRUE) < tol)) {
    message(paste("\nWarning: some rows of Q not numerically distinct from 0; setting to",
                  tol, "\n"))
    ii <- which(rowSums(Q, na.rm = TRUE) < tol)
    for (i in 1:length(ii))
      Q[ii[i], setdiff(1:ncol(Q), ii[i])] <- tol / (ncol(Q) - 1)
  }
  diag(Q) <- -rowSums(Q, na.rm = TRUE)
  return(list(Q = Q, L = L, loglik = logLik(obj), pi = pi))
}


# ---- Console message printer (internal) ----
#
# @keywords internal
printmessage <- function(Q, pi, method) {
  if (method == "empirical" || method == "fixed")
    cat("make.simmap is sampling character histories conditioned on\nthe transition matrix\n\nQ =\n")
  else if (method == "mcmc") {
    cat("make.simmap is simulating with a sample of Q from\nthe posterior distribution\n")
    cat("\nMean Q from the posterior is\nQ =\n")
  }
  print(Q)
  if (method == "empirical") cat("(estimated using likelihood);\n")
  else if (method == "fixed") cat("(specified by the user);\n")
  cat("and (mean) root node prior probabilities\npi =\n")
  if (is.list(pi)) pi <- Reduce("+", pi) / length(pi)
  print(pi)
  flush.console()
}


# ---- Gillespie algorithm path simulator (internal) ----
# Simulates a continuous-time Markov chain path along a branch.
#
# @keywords internal
sch <- function(start, t, Q) {
  tol <- t * 1e-12
  dt <- setNames(0, start)
  while (sum(dt) < (t - tol)) {
    # index Q by state NAME so arbitrary state labels work
    s <- names(dt)[length(dt)]
    # Draw waiting time from exponential distribution
    dt[length(dt)] <- if (-Q[s, s] > 0) rexp(n = 1, rate = -Q[s, s]) else t - sum(dt)
    if (sum(dt) < (t - tol)) {
      # Transition to a new state
      dt <- c(dt, 0)
      if (sum(Q[s, ][-match(s, colnames(Q))]) > 0)
        names(dt)[length(dt)] <- rstate(
          Q[s, ][-match(s, colnames(Q))] /
            sum(Q[s, ][-match(s, colnames(Q))]))
      else
        names(dt)[length(dt)] <- s
    } else {
      # Adjust last segment to fill remaining time exactly
      dt[length(dt)] <- dt[length(dt)] - sum(dt) + t
    }
  }
  return(dt)
}


# ---- Matrix exponentiation wrapper (internal) ----
# Uses the faster matexpo for symmetric matrices, falls back to expm.
#
# @keywords internal
EXPM <- function(x, ...) {
  e_x <- if (isSymmetric(x)) matexpo(x) else expm(x, ...)
  dimnames(e_x) <- dimnames(x)
  e_x
}
