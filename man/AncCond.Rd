\name{AncCond}
\alias{AncCond}
\title{Calculate the mean of a continuous character at transitions in a binary charachter}
\description{Tests whether transitions in a binary discrete trait are associated with
extreme ancestral values of a continuous trait. Null distribution generated
by holding the continuous trait fixed and simulating discrete traits under
the fitted Mk model.}
\usage{
AncCond(tree, 
        x, 
        y,
        n.maps = 100, 
        n.sims = 500, 
        prune = FALSE,
        Q = c(0,2,1,0),
        hypothesis = "none")
}
\arguments{
  \item{tree}{tree of class phylo}
  \item{x}{named numeric vector — continuous trait at tips}
    \item{y}{named integer vector — binary discrete trait (0/1) at tips}
    \item{n.maps}{number of stochastic maps per test-statistic computation}
    \item{n.sims}{number of null simulations}
    \item{prune}{if TRUE, estimate continuous trait using only state-0 lineages}
    \item{Q}{optional 2x2 rate matrix (rows/cols labelled "0","1"). If supplied, 
the Mk model is NOT fitted — the provided rates are used directly for stochastic mapping 
and null simulation. Useful when the user has a preferred model or wants to avoid
automatic model selection.}
    \item{hypothesis}{direction of the alternative hypothesis. One of:
"none"   — two-sided test (default): p = min(1, 2 * min(p.upper, p.lower))
"higher" — one-sided: transitions occur at HIGH continuous values (p = p.upper)
"lower"  — one-sided: transitions occur at LOW continuous values (p = p.lower)}

}
\details{This function uses ancestral state estimates to determine if the transitions in the binary trait are associated with extreme values of the continuous trait.  This test can incorporate the possibility that the derived state of a binary character may lead to correlated selection in the continuous trait.  If this is desired then the drop.state argument should be used to specify the derived state of the binary character that should not be used in the ancestral state estimation of the continuous trait.}

\value{Returns a list of 14 AncCond elements:
  \item{T.obs.01}{ observed mean continuous value at 0→1 transitions}
  \item{T.obs.10}{ observed mean continuous value at 1→0 transitions}
  \item{p.01}{ p-value for 0→1 transitions (adjusted per hypothesis)}
  \item{p.10}{ p-value for 1→0 transitions (adjusted per hypothesis)}
  \item{p.combined}{ Bonferroni-corrected p-value across both transition directions: min(1, 2 * min(p.01, p.10)). Use this when interpreting both p.01 and p.10 as independent tests.}
  \item{p.upper.01}{ one-sided P(T_null >= T_obs) for 0→1}
  \item{p.lower.01}{ one-sided P(T_null <= T_obs) for 0→1}
  \item{p.upper.10}{ one-sided P(T_null >= T_obs) for 1→0}
  \item{p.lower.10}{ one-sided P(T_null <= T_obs) for 1→0}
  \item{T.null.01}{ null distribution for 0→1 (length n.sims)}
  \item{T.null.10}{ null distribution for 1→0 (length n.sims)}
  \item{Q}{ fitted rate matrix}
  \item{model}{ "ARD", "unidirectional", or "user-supplied"}
  \item{hypothesis}{ hypothesis direction used}

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
}
}
\keyword{ comparative method }
\keyword{ discrete trait }
\keyword{ continuous trait }