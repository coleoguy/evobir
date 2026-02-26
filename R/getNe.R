#' Effective Population Size from Census Numbers
#'
#' Calculates the effective population size (Ne) from the number of
#' breeding males and females.  Ne is almost always smaller than the
#' actual census size because not all individuals contribute equally
#' to the next generation.  Unequal sex ratios reduce Ne further.
#'
#' Two formulas are available:
#' \itemize{
#'   \item \strong{Autosomal} (default): \eqn{N_e = 4MF / (M + F)}
#'   \item \strong{X-linked}: \eqn{N_e = 9MF / (4M + 2F)}
#' }
#' where M = number of breeding males, F = number of breeding females.
#'
#' @param males Integer.  Number of breeding males in the population.
#' @param females Integer.  Number of breeding females in the population.
#' @param locus Character.  Which locus type to calculate Ne for:
#'   \itemize{
#'     \item \code{"A"} (default) -- autosomal locus.
#'     \item \code{"X"} -- X-linked locus.
#'   }
#'
#' @return Numeric.  The effective population size.
#'
#' @details
#' When M = F, autosomal Ne = 2N (same as census).  As the sex ratio
#' becomes more skewed, Ne drops.  X-linked Ne differs because X
#' chromosomes spend 2/3 of their time in females and 1/3 in males.
#'
#' @examples
#' # Equal sex ratio: 50 males, 50 females
#' getNe(50, 50)
#' # Returns 100 (same as census)
#'
#' # Skewed: 10 males, 90 females (like elephant seals)
#' getNe(10, 90)
#' # Returns 36 -- much less than the census of 100!
#'
#' # X-linked Ne with equal sex ratio
#' getNe(50, 50, locus = "X")
#'
#' @export
getNe <- function(males, females, locus = "A") {
  if (locus == "A") {
    # Autosomal: classic formula from Wright (1931)
    ne <- (4 * males * females) / (males + females)
  }
  if (locus == "X") {
    # X-linked: accounts for hemizygosity in males
    ne <- (9 * males * females) / (4 * males + 2 * females)
  }
  return(ne)
}
