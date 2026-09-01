\name{AncCond}
\alias{AncCond}
\title{Calculate the mean of a continuous character at transitions in a binary charachter}
\description{Tests whether transitions in a binary discrete trait are associated with
extreme ancestral values of a continuous trait. Null distribution generated
by holding the continuous trait fixed and simulating discrete traits under
the fitted Mk model. The continuous and discrete traits may each be observed
on any subset of the tips, including disjoint subsets.}
\usage{
AncCond(tree,
        x,
        y,
        n.maps = 100,
        n.sims = 1000,
        prune = FALSE,
        Q = NULL,
        hypothesis = "none")
}
\arguments{
  \item{tree}{tree of class phylo}
  \item{x}{named numeric vector: continuous trait at tips. Tips may be absent
from \code{x} or coded \code{NA}; at least 3 tips must carry data.}
    \item{y}{named integer vector: binary discrete trait (0/1) at tips. Tips may
be absent from \code{y} or coded \code{NA}; both states must be observed in at
least one tip. \code{x} and \code{y} need not be observed on the same tips.}
    \item{n.maps}{number of stochastic maps per test-statistic computation}
    \item{n.sims}{number of null simulations}
    \item{prune}{if TRUE, estimate the continuous trait using only lineages
observed in state 0 (tips with unknown \code{y} are excluded as well)}
    \item{Q}{optional 2x2 rate matrix (rows/cols labelled "0","1"). If supplied,
the Mk model is NOT fitted; the provided rates are used directly for stochastic mapping
and null simulation. Useful when the user has a preferred model or wants to avoid
automatic model selection.}
    \item{hypothesis}{direction of the alternative hypothesis. One of:
"none": two-sided test (default): p = min(1, 2 * min(p.upper, p.lower))
"higher": one-sided: transitions occur at HIGH continuous values (p = p.upper)
"lower": one-sided: transitions occur at LOW continuous values (p = p.lower)}

}
\details{This function uses ancestral state estimates to determine if the transitions in the binary trait are associated with extreme values of the continuous trait.  This test can incorporate the possibility that the derived state of a binary character may lead to correlated selection in the continuous trait.  If this is desired then \code{prune = TRUE} restricts the ancestral state estimation of the continuous trait to lineages observed in state 0.

Missing data are handled as follows. Tips without \code{y} are marginalised
out of the Mk likelihood and of the stochastic maps (their state is sampled
from the posterior given the rest of the data), and the null simulations are
masked to the same missingness pattern so that observed and null test
statistics are computed from equivalent information. Tips and internal nodes
that lack a direct continuous-trait estimate receive their Brownian motion
conditional expectation given the tips that carry \code{x}: nodes lying on an
edge of the tree pruned to the \code{x}-bearing tips are linearly interpolated
between the estimates at the ends of that edge (weighted by branch length), and
nodes in subtrees with no \code{x} data take the value at the point where that
subtree attaches. Transitions inferred in lineages without \code{x} therefore
contribute an imputed value, and power falls as the overlap between the
\code{x}- and \code{y}-bearing tips shrinks because the signal then comes only
through phylogenetic signal in \code{x}.}

\value{Returns a list of class \code{anccond} with elements:
  \item{T.obs.01}{ observed mean continuous value at 0 to 1 transitions}
  \item{T.obs.10}{ observed mean continuous value at 1 to 0 transitions}
  \item{p.01}{ p-value for 0 to 1 transitions (adjusted per hypothesis)}
  \item{p.10}{ p-value for 1 to 0 transitions (adjusted per hypothesis)}
  \item{p.combined}{ Bonferroni-corrected p-value across both transition directions: min(1, 2 * min(p.01, p.10)). Use this when interpreting both p.01 and p.10 as independent tests.}
  \item{p.upper.01}{ one-sided P(T_null >= T_obs) for 0 to 1}
  \item{p.lower.01}{ one-sided P(T_null <= T_obs) for 0 to 1}
  \item{p.upper.10}{ one-sided P(T_null >= T_obs) for 1 to 0}
  \item{p.lower.10}{ one-sided P(T_null <= T_obs) for 1 to 0}
  \item{T.null.01}{ null distribution for 0 to 1 (length n.sims)}
  \item{T.null.10}{ null distribution for 1 to 0 (length n.sims)}
  \item{Q}{ fitted rate matrix}
  \item{model}{ "ARD", "unidirectional", or "user-supplied"}
  \item{hypothesis}{ hypothesis direction used}
  \item{n.obs}{ named vector: number of tips in the tree, with \code{x} observed, and with \code{y} observed}

}
\author{Meghann McConnell, Nathan Anderson, Richard H. Adams, and Heath Blackmon}
\references{http://coleoguy.github.io/}
\examples{
\dontrun{
# Simulate a toy dataset
set.seed(1)
tree <- pbtree(n = 100, scale = 1)
x <- fastBM(tree)                                   # continuous trait
y <- rbinom(100, 1, 0.3); names(y) <- tree$tip.label # binary trait

# Run the test (small settings for a quick demo)
res <- AncCond(tree, x, y, n.maps = 50, n.sims = 500)

# Continuous and discrete traits measured on different species
res2 <- AncCond(tree, x[1:60], y[41:100], n.maps = 50, n.sims = 500)
}
}
\keyword{ comparative method }
\keyword{ discrete trait }
\keyword{ continuous trait }
