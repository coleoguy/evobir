# ============================================================================
# AncCond: Ancestral Condition Test
#
# Tests whether transitions in a binary discrete trait are associated with
# extreme ancestral values of a continuous trait. Null distribution generated
# by holding the continuous trait fixed and simulating discrete traits under
# the fitted Mk model.
#
# Dependencies: ape, phytools (for fastAnc only)
# ============================================================================


# ── Internal: analytical 2-state transition probabilities ──────────────────
# Returns 2x2 matrix P[from+1, to+1] = P(to | from, t)
# States are 0-indexed but matrix is 1-indexed.
.tp2 <- function(q01, q10, t) {
  if (q01 + q10 < 1e-15) return(diag(2))
  if (q10 < 1e-15) {
    e <- exp(-q01 * t)
    return(matrix(c(e, 1 - e, 0, 1), 2, 2, byrow = TRUE))
  }
  tot <- q01 + q10
  e   <- exp(-tot * t)
  matrix(c(
    (q10 + q01 * e) / tot,  q01 * (1 - e) / tot,
     q10 * (1 - e) / tot,  (q01 + q10 * e) / tot
  ), 2, 2, byrow = TRUE)
}


# ── Internal: pruning algorithm (postorder) with log-rescaling ─────────────
# tree must already be reorder(tree, "postorder").
# y: named 0/1 vector at tips (length n.tips).
# Returns a list with:
#   cl    — matrix cl[node, state+1] of rescaled conditional likelihoods
#           (proportional within each node, suitable for posterior sampling)
#   logL  — log-likelihood at root under flat prior
.pruning2 <- function(tree, y, q01, q10) {
  n.tips  <- length(tree$tip.label)
  n.total <- n.tips + tree$Nnode
  root    <- n.tips + 1L
  edges   <- tree$edge
  el      <- tree$edge.length

  cl <- matrix(1, n.total, 2)
  # tips: indicator at observed state
  cl[1:n.tips, ] <- 0
  for (i in seq_len(n.tips)) cl[i, y[i] + 1L] <- 1

  # log-scaling constant accumulated over internal nodes
  log.scale <- 0

  for (k in seq_len(nrow(edges))) {
    pa <- edges[k, 1L]
    ch <- edges[k, 2L]
    P  <- .tp2(q01, q10, el[k])
    # contribution of child to parent's conditional likelihood
    c0 <- P[1, 1] * cl[ch, 1] + P[1, 2] * cl[ch, 2]
    c1 <- P[2, 1] * cl[ch, 1] + P[2, 2] * cl[ch, 2]
    cl[pa, 1] <- cl[pa, 1] * c0
    cl[pa, 2] <- cl[pa, 2] * c1

    # rescale at parent to prevent underflow
    mx <- max(cl[pa, 1], cl[pa, 2])
    if (mx > 0 && is.finite(mx)) {
      cl[pa, ] <- cl[pa, ] / mx
      log.scale <- log.scale + log(mx)
    }
  }

  # log-likelihood under flat root prior
  logL <- log(sum(cl[root, ] * 0.5)) + log.scale

  list(cl = cl, logL = logL)
}


# ── Internal: fit 2-state Mk model ────────────────────────────────────────
# Fits bidirectional (ARD) and unidirectional (q10=0) models, returns the
# AIC-preferred fit.  ARD uses L-BFGS-B with a rate floor to prevent
# degenerate solutions where one rate collapses to ~0.
.fitMk2 <- function(tree, y, q01, q10) {
  n.tips <- length(tree$tip.label)
  root   <- n.tips + 1L
  rate.floor <- 1e-6   # minimum rate (keep > 0 for numerical stability)

  # log-likelihood given rates (uses log-rescaled pruning)
  ll <- function(q01, q10) {
    .pruning2(tree, y, q01, q10)$logL
  }

  # ARD: optimise q01 and q10 jointly (bounded)
  # Use multiple starting points to avoid local optima
  nll.ard <- function(par) -ll(par[1], par[2])
  starts <- list(c(0.5, 0.5), c(0.01, 0.01), c(0.001, 0.001), c(0.001, 0.01))
  best.ard <- list(value = Inf)
  for (s in starts) {
    fit <- tryCatch(
      optim(s, nll.ard, method = "L-BFGS-B",
            lower = c(rate.floor, rate.floor),
            upper = c(100, 100),
            control = list(maxit = 5000)),
      error = function(e) list(value = Inf)
    )
    if (fit$value < best.ard$value) best.ard <- fit
  }
  fit.ard <- best.ard

  # Unidirectional: only q01
  nll.uni <- function(par) -ll(par, 0)
  fit.uni <- optimize(nll.uni, interval = c(rate.floor, 100))

  aic.ard <- 2 * fit.ard$value + 4  # 2k where k=2
  aic.uni <- 2 * fit.uni$objective + 2  # k=1

  # prefer simpler unless ARD wins by >=4 (conservative: unidirectional
  # models with q10=0 can produce degenerate stochastic maps)
  if (aic.ard + 4 <= aic.uni) {
    list(q01 = fit.ard$par[1], q10 = fit.ard$par[2],
         loglik = -fit.ard$value, AIC = aic.ard, model = "ARD")
  } else {
    list(q01 = fit.uni$minimum, q10 = 0,
         loglik = -fit.uni$objective, AIC = aic.uni, model = "unidirectional")
  }
}


# ── Internal: simulate discrete tip states (batch) ─────────────────────────
# Simulates n.sims independent discrete trait datasets under Mk(q01,q10).
# tree must be in cladewise (preorder) order.
# Returns n.tips × n.sims integer matrix of 0/1 tip states.
.simMk2 <- function(tree, q01, q10, n.sims) {
  n.tips  <- length(tree$tip.label)
  n.total <- n.tips + tree$Nnode
  root    <- n.tips + 1L
  edges   <- tree$edge
  el      <- tree$edge.length

  states <- matrix(NA_integer_, n.total, n.sims)

  # root from equilibrium (or state 0 if unidirectional)
  if (q10 < 1e-15) {
    states[root, ] <- 0L
  } else {
    pi1 <- q01 / (q01 + q10)
    states[root, ] <- rbinom(n.sims, 1L, pi1)
  }

  for (k in seq_len(nrow(edges))) {
    pa <- edges[k, 1L]
    ch <- edges[k, 2L]
    P  <- .tp2(q01, q10, el[k])
    p1.g0 <- P[1, 2]
    p1.g1 <- P[2, 2]
    prob1 <- ifelse(states[pa, ] == 0L, p1.g0, p1.g1)
    states[ch, ] <- rbinom(n.sims, 1L, prob1)
  }

  out <- states[1:n.tips, , drop = FALSE]
  rownames(out) <- tree$tip.label
  out
}


# ── Internal: sample ancestral states + compute test statistics ────────────
# Samples n.maps stochastic ancestral-state configurations via preorder
# traversal, identifies transition edges in both directions, and returns
# the mean interpolated continuous trait value at transition points,
# averaged across maps.
#
# Returns a list with:
#   T.01 — mean interpolated continuous value at 0→1 transitions
#   T.10 — mean interpolated continuous value at 1→0 transitions
#
# tree must be in cladewise order.
# anc.x: numeric vector indexed by node number (1:n.total), continuous trait
#         values at tips and ML estimates at internal nodes.
# clik:  rescaled conditional likelihood matrix from .pruning2()$cl
#         (n.total × 2). Only ratios matter for sampling, so rescaling is fine.
.anccond.T <- function(tree, anc.x, clik, q01, q10, n.maps) {
  n.tips  <- length(tree$tip.label)
  n.total <- n.tips + tree$Nnode
  root    <- n.tips + 1L
  edges   <- tree$edge
  el      <- tree$edge.length
  n.edges <- nrow(edges)

  # ── Sample root (vectorised across maps) ──
  rp <- clik[root, ] * 0.5
  rp <- rp / sum(rp)
  states <- matrix(NA_integer_, n.total, n.maps)
  states[root, ] <- rbinom(n.maps, 1L, rp[2])

  # ── Preorder traversal: sample child states ──
  p1g0 <- numeric(n.edges)
  p1g1 <- numeric(n.edges)

  for (k in seq_len(n.edges)) {
    ch <- edges[k, 2L]
    P  <- .tp2(q01, q10, el[k])
    d0 <- P[1, 1] * clik[ch, 1] + P[1, 2] * clik[ch, 2]
    d1 <- P[2, 1] * clik[ch, 1] + P[2, 2] * clik[ch, 2]
    p1g0[k] <- if (d0 > 0) P[1, 2] * clik[ch, 2] / d0 else 0.5
    p1g1[k] <- if (d1 > 0) P[2, 2] * clik[ch, 2] / d1 else 0.5
  }

  for (k in seq_len(n.edges)) {
    pa <- edges[k, 1L]
    ch <- edges[k, 2L]
    prob1 <- ifelse(states[pa, ] == 0L, p1g0[k], p1g1[k])
    states[ch, ] <- rbinom(n.maps, 1L, prob1)
  }

  # ── Identify transitions on each edge × map ──
  pa.st <- states[edges[, 1], , drop = FALSE]  # n.edges × n.maps
  ch.st <- states[edges[, 2], , drop = FALSE]
  is.01 <- (pa.st == 0L) & (ch.st == 1L)       # logical matrix
  is.10 <- (pa.st == 1L) & (ch.st == 0L)

  # ── Interpolated values at transition points ──
  # For 0→1: transition position drawn from truncated Exp(q01)
  # For 1→0: transition position drawn from truncated Exp(q10)
  u.mat <- matrix(runif(n.edges * n.maps), n.edges, n.maps)

  x.pa <- anc.x[edges[, 1]]   # n.edges
  x.ch <- anc.x[edges[, 2]]

  # 0→1 interpolation
  if (q01 > 1e-15) {
    d.mat.01 <- -log(1 - u.mat * (1 - exp(-q01 * el))) / (q01 * el)
  } else {
    d.mat.01 <- matrix(0.5, n.edges, n.maps)
  }
  d.mat.01[d.mat.01 < 0] <- 0; d.mat.01[d.mat.01 > 1] <- 1
  interp.01 <- x.pa * (1 - d.mat.01) + x.ch * d.mat.01

  # 1→0 interpolation
  if (q10 > 1e-15) {
    d.mat.10 <- -log(1 - u.mat * (1 - exp(-q10 * el))) / (q10 * el)
  } else {
    d.mat.10 <- matrix(0.5, n.edges, n.maps)
  }
  d.mat.10[d.mat.10 < 0] <- 0; d.mat.10[d.mat.10 > 1] <- 1
  interp.10 <- x.pa * (1 - d.mat.10) + x.ch * d.mat.10

  # ── Mean ancestral condition at transition points, per map ──
  T.01 <- numeric(n.maps)
  T.10 <- numeric(n.maps)
  for (m in seq_len(n.maps)) {
    v01 <- interp.01[is.01[, m], m]
    v10 <- interp.10[is.10[, m], m]
    T.01[m] <- if (length(v01) >= 1L) mean(v01) else NA_real_
    T.10[m] <- if (length(v10) >= 1L) mean(v10) else NA_real_
  }
  list(T.01 = mean(T.01, na.rm = TRUE),
       T.10 = mean(T.10, na.rm = TRUE))
}


# ── Internal: map pruned-tree ancestral estimates to full tree ─────────────
.map.pruned <- function(full.tree, pruned.tree, pruned.anc, obs.x) {
  n.f <- length(full.tree$tip.label)
  n.total.f <- n.f + full.tree$Nnode
  n.p <- length(pruned.tree$tip.label)

  anc.x <- numeric(n.total.f)
  anc.x[1:n.f] <- obs.x[full.tree$tip.label]
  anc.x[(n.f + 1):n.total.f] <- NA_real_

  # map each pruned-tree internal node to the corresponding full-tree node
  for (nd in (n.p + 1):(n.p + pruned.tree$Nnode)) {
    desc <- ape::extract.clade(pruned.tree, nd)$tip.label
    if (length(desc) >= 2L) {
      fnode <- ape::getMRCA(full.tree, desc)
    } else {
      fnode <- full.tree$edge[
        which(full.tree$edge[, 2] == match(desc, full.tree$tip.label)), 1]
    }
    nm <- as.character(nd)
    if (nm %in% names(pruned.anc)) anc.x[fnode] <- pruned.anc[nm]
  }

  # fallback: any unmapped internal nodes get the full-tree estimate
  na.nodes <- which(is.na(anc.x))
  if (length(na.nodes) > 0L) {
    full.anc <- phytools::fastAnc(full.tree, obs.x[full.tree$tip.label])
    for (nd in na.nodes) {
      nm <- as.character(nd)
      if (nm %in% names(full.anc)) anc.x[nd] <- full.anc[nm]
    }
  }
  anc.x
}


# ============================================================================
# AncCond — main function
# ============================================================================
#' Ancestral Condition Test
#'
#' Tests whether transitions in a binary discrete trait are associated with
#' extreme ancestral values of a continuous trait.
#'
#' @param tree   a phylo object (ape)
#' @param x      named numeric vector — continuous trait at tips
#' @param y      named integer vector — binary discrete trait (0/1) at tips
#' @param n.maps number of stochastic maps per test-statistic computation
#' @param n.sims number of null simulations
#' @param prune  if TRUE, estimate continuous trait using only state-0 lineages
#' @param Q      optional 2x2 rate matrix (rows/cols labelled "0","1"). If
#'               supplied, the Mk model is NOT fitted — the provided rates
#'               are used directly for stochastic mapping and null simulation.
#'               Useful when the user has a preferred model or wants to avoid
#'               automatic model selection.
#' @param hypothesis direction of the alternative hypothesis. One of:
#'               "none"   — two-sided test (default): p = min(1, 2 * min(p.upper, p.lower))
#'               "higher" — one-sided: transitions occur at HIGH continuous values (p = p.upper)
#'               "lower"  — one-sided: transitions occur at LOW continuous values (p = p.lower)
#'
#' @return A list of class "anccond" with elements:
#'   T.obs.01   — observed mean continuous value at 0→1 transitions
#'   T.obs.10   — observed mean continuous value at 1→0 transitions
#'   p.01       — p-value for 0→1 transitions (adjusted per hypothesis)
#'   p.10       — p-value for 1→0 transitions (adjusted per hypothesis)
#'   p.combined — Bonferroni-corrected p-value across both transition
#'                directions: min(1, 2 * min(p.01, p.10)). Use this when
#'                interpreting both p.01 and p.10 as independent tests.
#'   p.upper.01 — one-sided P(T_null >= T_obs) for 0→1
#'   p.lower.01 — one-sided P(T_null <= T_obs) for 0→1
#'   p.upper.10 — one-sided P(T_null >= T_obs) for 1→0
#'   p.lower.10 — one-sided P(T_null <= T_obs) for 1→0
#'   T.null.01  — null distribution for 0→1 (length n.sims)
#'   T.null.10  — null distribution for 1→0 (length n.sims)
#'   Q          — fitted rate matrix
#'   model      — "ARD", "unidirectional", or "user-supplied"
#'   hypothesis — hypothesis direction used
#'
#' @examples
#' \dontrun{
#'   library(ape); library(phytools)
#'   tree <- pbtree(n = 100, scale = 1)
#'   x <- fastBM(tree)
#'   y <- rbinom(100, 1, 0.3); names(y) <- tree$tip.label
#'   res <- AncCond(tree, x, y, n.maps = 50, n.sims = 500)
#'   res$p.upper.01   # 0->1 transitions at high values
#'   res$p.upper.10   # 1->0 transitions at high values
#' }
AncCond <- function(tree, x, y,
                    n.maps     = 100,
                    n.sims     = 1000,
                    prune      = FALSE,
                    Q          = NULL,
                    hypothesis = "none") {

  if (!inherits(tree, "phylo")) stop("tree must be a phylo object")
  hypothesis <- match.arg(hypothesis, c("none", "lower", "higher"))
  n.tips  <- length(tree$tip.label)
  n.total <- n.tips + tree$Nnode

  # ── Match data to tree ──
  if (is.null(names(x)) || is.null(names(y)))
    stop("x and y must be named vectors matching tree$tip.label")
  x <- x[tree$tip.label]
  y <- y[tree$tip.label]
  if (any(is.na(x)) || any(is.na(y)))
    stop("All tips must have non-missing values for x and y")
  y <- as.integer(y)
  if (!all(y %in% 0:1)) stop("y must be binary (0 and 1)")

  # ── Step 1: Continuous trait ancestral reconstruction ──
  if (prune) {
    keep <- tree$tip.label[y == 0L]
    if (length(keep) < 3L) stop("Pruning leaves fewer than 3 tips in state 0")
    pruned.tree <- ape::drop.tip(tree, setdiff(tree$tip.label, keep))
    pruned.anc  <- phytools::fastAnc(pruned.tree, x[keep])
    anc.x <- .map.pruned(tree, pruned.tree, pruned.anc, x)
  } else {
    anc.est <- phytools::fastAnc(tree, x)
    anc.x <- numeric(n.total)
    anc.x[1:n.tips] <- x
    anc.x[(n.tips + 1):n.total] <- anc.est[as.character((n.tips + 1):n.total)]
  }

  # ── Step 2: Mk rates (fit or user-supplied) ──
  tree.po <- reorder(tree, "postorder")
  if (!is.null(Q)) {
    # user-supplied Q matrix: extract rates, skip model selection
    if (!is.matrix(Q) || !all(dim(Q) == c(2, 2)))
      stop("Q must be a 2x2 rate matrix")
    q01 <- Q[1, 2]
    q10 <- Q[2, 1]
    if (q01 < 0 || q10 < 0) stop("Q rates must be non-negative")
    mk <- list(q01 = q01, q10 = q10, loglik = NA, AIC = NA,
               model = "user-supplied")
  } else {
    mk <- .fitMk2(tree.po, y, q01 = NA, q10 = NA)
    q01 <- mk$q01
    q10 <- mk$q10
  }

  # ── Step 3–5: Empirical test statistics (both directions) ──
  tree.pre <- reorder(tree, "cladewise")
  prun.obs <- .pruning2(tree.po, y, q01, q10)
  T.obs    <- .anccond.T(tree.pre, anc.x, prun.obs$cl, q01, q10, n.maps)

  T.obs.01 <- T.obs$T.01
  T.obs.10 <- T.obs$T.10

  # ── Guard: no transitions in observed data ──
  if (!is.finite(T.obs.01) && !is.finite(T.obs.10)) {
    warning("No transitions detected in stochastic maps of observed data; ",
            "returning NA p-values. Consider whether the binary trait has ",
            "enough transitions for a meaningful test.")
    Qmat <- matrix(c(-q01, q01, q10, -q10), 2, 2,
                   dimnames = list(c("0","1"), c("0","1")))
    out <- list(
      T.obs.01   = NA_real_, T.obs.10   = NA_real_,
      p.01       = NA_real_, p.10       = NA_real_,
      p.combined = NA_real_,
      p.upper.01 = NA_real_, p.lower.01 = NA_real_,
      p.upper.10 = NA_real_, p.lower.10 = NA_real_,
      T.null.01  = rep(NA_real_, n.sims),
      T.null.10  = rep(NA_real_, n.sims),
      Q = Qmat, model = mk$model, n.maps = n.maps, n.sims = n.sims,
      prune = prune, hypothesis = hypothesis
    )
    class(out) <- "anccond"
    return(out)
  }

  # ── Null distribution ──
  y.sims <- .simMk2(tree.pre, q01, q10, n.sims)

  T.null.01 <- numeric(n.sims)
  T.null.10 <- numeric(n.sims)
  for (b in seq_len(n.sims)) {
    y.b   <- y.sims[, b]
    prun.b <- .pruning2(tree.po, y.b, q01, q10)
    T.b   <- .anccond.T(tree.pre, anc.x, prun.b$cl, q01, q10, n.maps)
    T.null.01[b] <- T.b$T.01
    T.null.10[b] <- T.b$T.10
  }

  # ── p-values (with continuity correction) ──
  .pvals <- function(T.obs, T.null) {
    T.valid <- T.null[is.finite(T.null)]
    nv <- length(T.valid)
    if (nv < 10L) {
      warning("Fewer than 10 valid null replicates (", nv,
              "); p-values may be unreliable.", call. = FALSE)
    }
    if (nv == 0L || !is.finite(T.obs)) {
      return(list(p.upper = NA_real_, p.lower = NA_real_))
    }
    list(
      p.upper = (sum(T.valid >= T.obs) + 1) / (nv + 1),
      p.lower = (sum(T.valid <= T.obs) + 1) / (nv + 1)
    )
  }

  pv.01 <- .pvals(T.obs.01, T.null.01)
  pv.10 <- .pvals(T.obs.10, T.null.10)

  # ── Compute final p-values based on hypothesis direction ──
  .final.p <- function(p.upper, p.lower, hypothesis) {
    switch(hypothesis,
      higher = p.upper,
      lower  = p.lower,
      none   = min(1, 2 * min(p.upper, p.lower))
    )
  }
  p.01 <- .final.p(pv.01$p.upper, pv.01$p.lower, hypothesis)
  p.10 <- .final.p(pv.10$p.upper, pv.10$p.lower, hypothesis)

  # Bonferroni correction for testing both transition directions
  p.combined <- min(1, 2 * min(p.01, p.10, na.rm = TRUE))
  if (!is.finite(p.combined)) p.combined <- NA_real_

  Qmat <- matrix(c(-q01, q01, q10, -q10), 2, 2,
                 dimnames = list(c("0","1"), c("0","1")))
  out <- list(
    T.obs.01   = T.obs.01,
    T.obs.10   = T.obs.10,
    p.01       = p.01,
    p.10       = p.10,
    p.combined = p.combined,
    p.upper.01 = pv.01$p.upper,
    p.lower.01 = pv.01$p.lower,
    p.upper.10 = pv.10$p.upper,
    p.lower.10 = pv.10$p.lower,
    T.null.01  = T.null.01,
    T.null.10  = T.null.10,
    Q          = Qmat,
    model      = mk$model,
    n.maps     = n.maps,
    n.sims     = n.sims,
    prune      = prune,
    hypothesis = hypothesis
  )
  class(out) <- "anccond"
  out
}


# ── print method ───────────────────────────────────────────────────────────
print.anccond <- function(x, ...) {
  hyp.label <- switch(x$hypothesis,
    none   = "two-sided (2 * min of one-sided p-values)",
    higher = "one-sided (transitions at HIGH values)",
    lower  = "one-sided (transitions at LOW values)"
  )
  cat("\nAncestral Condition Test\n")
  cat("─────────────────────────────────────\n")
  cat("Mk model:    ", x$model, "\n")
  cat("  q01 =", formatC(x$Q[1,2], digits = 4, format = "f"), "\n")
  cat("  q10 =", formatC(x$Q[2,1], digits = 4, format = "f"), "\n")
  cat("Maps/test:   ", x$n.maps, "\n")
  cat("Null sims:   ", x$n.sims, "\n")
  cat("Pruned:      ", x$prune, "\n")
  cat("Hypothesis:  ", hyp.label, "\n\n")
  cat("0→1 transitions (gains):\n")
  cat("  T.obs =", formatC(x$T.obs.01, digits = 4, format = "f"), "\n")
  cat("  p     =", formatC(x$p.01, digits = 4, format = "f"), "\n")
  cat("  p (upper) =", formatC(x$p.upper.01, digits = 4, format = "f"),
      "  p (lower) =", formatC(x$p.lower.01, digits = 4, format = "f"), "\n\n")
  cat("1→0 transitions (losses):\n")
  cat("  T.obs =", formatC(x$T.obs.10, digits = 4, format = "f"), "\n")
  cat("  p     =", formatC(x$p.10, digits = 4, format = "f"), "\n")
  cat("  p (upper) =", formatC(x$p.upper.10, digits = 4, format = "f"),
      "  p (lower) =", formatC(x$p.lower.10, digits = 4, format = "f"), "\n\n")
  cat("Combined (Bonferroni across directions):\n")
  cat("  p.combined =", formatC(x$p.combined, digits = 4, format = "f"), "\n")
  cat("  Note: use p.combined if interpreting both p.01\n")
  cat("  and p.10 as independent tests from the same data.\n")
  cat("─────────────────────────────────────\n")
  invisible(x)
}
