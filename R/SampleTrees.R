#' Sample Trees from Bayesian MCMC Output
#'
#' Reads a collection of phylogenetic trees from a NEXUS file (typically
#' output from a Bayesian analysis like MrBayes or BEAST), removes the
#' burn-in, and randomly samples a specified number of post-burn-in trees.
#'
#' \strong{Why do you need this?}  Bayesian phylogenetic programs like
#' MrBayes produce thousands (sometimes millions) of trees as they
#' explore "tree space."  The first chunk of trees are unreliable
#' because the program was still searching for good trees -- this
#' early period is called the "burn-in" and should be discarded.
#' After removing the burn-in, you typically don't need ALL the
#' remaining trees (they would be too many to work with), so you
#' randomly subsample a manageable number that still represents the
#' full range of uncertainty in your phylogeny.
#'
#' @param trees Character.  Path to a NEXUS format file containing
#'   phylogenetic trees (e.g., the \code{.t} file from MrBayes or
#'   the \code{.trees} file from BEAST).
#' @param burnin Numeric between 0 and 1.  Proportion of trees to
#'   discard as burn-in.  For example, \code{0.25} removes the first
#'   25\% of trees.  A common starting point is 0.25, but you should
#'   check convergence in Tracer or a similar tool to choose wisely.
#' @param final.number Integer.  Number of trees to randomly sample from
#'   the post-burn-in trees.  Common choices are 100 or 1000.
#' @param format Character.  Output format:
#'   \itemize{
#'     \item \code{"new"} -- Newick format (.nwk file)
#'     \item \code{"nex"} -- Nexus format (.nex file)
#'   }
#' @param prefix Character.  Prefix for the output file name.  The
#'   appropriate extension (.nwk or .nex) is added automatically.
#'
#' @return Writes the sampled trees to a file.  No value is returned
#'   (called for its side effect of writing a file).
#'
#' @examples
#' \dontrun{
#' # Your MrBayes run produced "my_analysis.nex.t" with 10,000 trees.
#' # Remove the first 25% as burn-in, then randomly pick 100 trees:
#' SampleTrees("my_analysis.nex.t", burnin = 0.25,
#'             final.number = 100, format = "new",
#'             prefix = "sampled_trees")
#' # This creates "sampled_trees.nwk" with 100 trees
#' }
#'
#' @importFrom ape read.nexus write.tree write.nexus
#' @export
SampleTrees <- function(trees, burnin, final.number, format, prefix) {
  # Read the full collection of trees from the NEXUS file
  trees <- read.nexus(trees)
  original.number <- length(trees)

  # Validate burnin
  if (!is.numeric(burnin) || burnin < 0 || burnin >= 1) {
    stop("burnin must be a proportion in [0, 1)")
  }

  # Remove the burn-in portion (first floor(burnin*N) trees are
  # discarded; the sample is drawn from the remaining trees)
  post.burnin.trees <- trees[(floor(burnin * original.number) + 1):original.number]

  # Validate the requested sample size
  if (final.number > length(post.burnin.trees)) {
    stop("final.number (", final.number,
         ") exceeds the number of post-burn-in trees (",
         length(post.burnin.trees), ")")
  }

  # Randomly sample the desired number of post-burn-in trees
  final.trees <- sample(post.burnin.trees, final.number)

  # Save in the requested format
  if (format == "new") {
    write.tree(final.trees, file = paste0(prefix, ".nwk"))
    print("Your trees were saved in the Newick format")
  }
  if (format == "nex") {
    write.nexus(final.trees, file = paste0(prefix, ".nex"))
    print("Your trees were saved in the Nexus format")
  }
}
