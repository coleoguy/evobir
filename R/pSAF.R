#' Probability of Sex-Autosome Fusion
#'
#' Computes the expected proportion of chromosomal fusions that join a sex
#' chromosome to an autosome (SA-fusions), as opposed to autosome-autosome
#' (AA) or sex-sex (SS) fusions. Two null models are available:
#'
#' \strong{Occurrence model} (\code{weighted = FALSE}): The original null from
#' Anderson et al. (2020). Every pair of non-homologous chromosomes is equally
#' likely to fuse. The probability of each fusion class (AA, SA, SS) depends
#' only on how many chromosomes of each type are present in each sex, weighted
#' by the proportion of fusions originating in each sex (\code{mu}).
#'
#' \strong{Fixation-weighted model} (\code{weighted = TRUE}, default): Extends
#' the occurrence model by accounting for the fact that a new fusion's
#' probability of drifting to fixation under neutrality depends on how many
#' copies of that chromosome exist in the population (Kimura 1962). Autosomes
#' are present in 2N copies, the shared sex chromosome (X or Z) in 3N/2
#' copies, and the limited sex chromosome (Y or W) in only N/2 copies. Because
#' fewer copies means higher fixation probability, Y-autosome (or W-autosome)
#' fusions are disproportionately likely to reach fixation relative to their
#' occurrence rate. This model therefore predicts a higher baseline proportion
#' of SA-fusions than the occurrence model alone.
#'
#' @param n_auto Integer >= 1. The haploid autosome count — i.e., the number
#'   of autosome PAIRS. For example, humans have 22 autosome pairs, so
#'   \code{n_auto = 22}. The diploid autosome count used internally is
#'   \code{2 * n_auto}. Can be a vector to compute probabilities across a
#'   range of autosome counts simultaneously.
#'
#' @param n_x Integer >= 0. Number of X chromosomes present in the
#'   heterogametic sex (males in XY systems). For a standard XY system,
#'   \code{n_x = 1}. For an XXY system, \code{n_x = 2}. For an XO system,
#'   \code{n_x = 1} and \code{n_y = 0}. Must not be used together with
#'   \code{n_z} or \code{n_w}.
#'
#' @param n_y Integer >= 0. Number of Y chromosomes present in the
#'   heterogametic sex. For a standard XY system, \code{n_y = 1}. For an
#'   XYY system, \code{n_y = 2}. Set to 0 for XO systems. Must not be used
#'   together with \code{n_z} or \code{n_w}.
#'
#' @param n_z Integer >= 0. Number of Z chromosomes present in the
#'   heterogametic sex (females in ZW systems). Analogous to \code{n_x} but
#'   for ZW systems. Must not be used together with \code{n_x} or \code{n_y}.
#'
#' @param n_w Integer >= 0. Number of W chromosomes present in the
#'   heterogametic sex. Analogous to \code{n_y} but for ZW systems. Must not
#'   be used together with \code{n_x} or \code{n_y}.
#'
#' @param mu Numeric between 0 and 1. The proportion of fusion mutations that
#'   originate in \strong{females}, regardless of which sex is homogametic.
#'   Default is 0.5, meaning both sexes contribute equally. Setting
#'   \code{mu = 0.7} always means 70\% of fusions arise in females — in an
#'   XY system that is the homogametic sex, but in a ZW system it is the
#'   heterogametic sex. The function handles this mapping internally so the
#'   user can always think in terms of phenotypic sex.
#'
#' @param weighted Logical. If \code{TRUE} (default), returns the fixation-
#'   weighted probability of SA-fusions (Schmalz & Blackmon). If \code{FALSE},
#'   returns the occurrence probability only (Anderson et al. 2020).
#'
#' @return A numeric value (or vector, if \code{n_auto} is a vector) giving
#'   the probability that a random fusion is a sex-autosome fusion under the
#'   chosen null model.
#'
#' @details
#' \strong{Terminology:}
#' \itemize{
#'   \item \emph{Shared sex chromosome}: the sex chromosome present in both
#'     sexes — X in XY systems (females are XX, males have at least one X) or
#'     Z in ZW systems (males are ZZ, females have at least one Z).
#'   \item \emph{Limited sex chromosome}: the sex chromosome present only in
#'     the heterogametic sex — Y in XY systems or W in ZW systems.
#'   \item \emph{Homogametic sex}: the sex with two copies of the shared sex
#'     chromosome (e.g., XX females in XY systems, ZZ males in ZW systems).
#'   \item \emph{Heterogametic sex}: the sex carrying both shared and limited
#'     sex chromosomes (e.g., XY males, ZW females).
#'   \item \emph{Non-homologous fusion}: a fusion between two chromosomes that
#'     are NOT homologs of each other. Fusions between homologs (e.g., both
#'     copies of chromosome 3) produce unbalanced gametes and are excluded.
#' }
#'
#' \strong{How the math works:}
#'
#' Within each sex, the probability that a random non-homologous fusion falls
#' into each class (AA, SA, SS) is determined by the number of chromosomes of
#' each type. The overall probability is the weighted average across the two
#' sexes (using \code{mu}).
#'
#' For the fixation-weighted model, the SA class is further split into fusions
#' involving the shared sex chromosome (XA or ZA) versus the limited sex
#' chromosome (YA or WA). Each class is multiplied by a fixation weight
#' proportional to 1/(population copy number), then renormalized so
#' probabilities sum to 1.
#'
#' \strong{Fixation weights (relative to autosomes):}
#' \itemize{
#'   \item Autosome: weight = 1 (present in 2N copies)
#'   \item Shared sex chrom (X/Z): weight = 4/3 (present in 3N/2 copies)
#'   \item Limited sex chrom (Y/W): weight = 4 (present in N/2 copies)
#' }
#'
#' @references
#' Anderson, N.W., Hjelmen, C.E., & Blackmon, H. (2020). The probability of
#' fusions joining sex chromosomes and autosomes. \emph{Biology Letters}, 16,
#' 20200648. \doi{10.1098/rsbl.2020.0648}
#'
#' Schmalz, S. & Blackmon, H. A fixation-weighted null model for the
#' proportion of sex chromosome-autosome fusions. \emph{In prep}.
#'
#' Kimura, M. (1962). On the probability of fixation of mutant genes in a
#' population. \emph{Genetics}, 47, 713-719.
#'
#' @examples
#' # ---- Basic usage ----
#'
#' # Standard XY system with 5 autosome pairs (Da = 10), fixation-weighted
#' pSAF(n_auto = 5, n_x = 1, n_y = 1)
#'
#' # Same system, occurrence probability only (Anderson et al. 2020)
#' pSAF(n_auto = 5, n_x = 1, n_y = 1, weighted = FALSE)
#'
#' # ---- Different sex chromosome systems ----
#'
#' # XO system (no Y chromosome): 13 autosome pairs
#' pSAF(n_auto = 13, n_x = 1, n_y = 0)
#'
#' # XXY system: 2 X chromosomes and 1 Y in the heterogametic sex
#' pSAF(n_auto = 5, n_x = 2, n_y = 1)
#'
#' # XYY system: 1 X and 2 Y chromosomes in the heterogametic sex
#' pSAF(n_auto = 5, n_x = 1, n_y = 2)
#'
#' # ---- ZW systems ----
#'
#' # Standard ZW (e.g., birds): 10 autosome pairs
#' pSAF(n_auto = 10, n_z = 1, n_w = 1)
#'
#' # ---- Vectorized over autosome count ----
#'
#' # Compute across a range of autosome counts for plotting
#' n <- 1:30
#' p_fix <- pSAF(n_auto = n, n_x = 1, n_y = 1)
#' p_occ <- pSAF(n_auto = n, n_x = 1, n_y = 1, weighted = FALSE)
#' plot(n, p_fix, type = "l", col = "red",
#'      xlab = "Haploid autosome count", ylab = "P(SA fusion)")
#' lines(n, p_occ, col = "blue", lty = 2)
#' legend("topright", c("Fixation-weighted", "Occurrence"),
#'        col = c("red", "blue"), lty = 1:2)
#'
#' # ---- Habronattus example from the paper (XXO, Da = 26) ----
#' pSAF(n_auto = 13, n_x = 2, n_y = 0)                  # fixation: 0.243
#' pSAF(n_auto = 13, n_x = 2, n_y = 0, weighted = FALSE) # occurrence: 0.194
#'
#' @export
pSAF <- function(n_auto,
                 n_x = 0, n_y = 0,
                 n_z = 0, n_w = 0,
                 mu = 0.5, weighted = TRUE) {

  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------

  # Determine which sex chromosome system the user specified. Exactly one
  # system (XY-type or ZW-type) must be provided — never both.
  xy <- n_x > 0 || n_y > 0   # TRUE if user specified any X or Y chromosomes
  zw <- n_z > 0 || n_w > 0   # TRUE if user specified any Z or W chromosomes
  if (xy && zw)   stop("cannot mix XY and ZW systems")
  if (!xy && !zw) stop("must specify at least one sex chromosome")

  # Internally we use generic names:
  #   s = count of the SHARED sex chromosome (X in XY systems, Z in ZW systems)
  #       This is the one present in BOTH sexes. The homogametic sex carries
  #       2 * s copies; the heterogametic sex carries s copies.
  #   l = count of the LIMITED sex chromosome (Y in XY systems, W in ZW systems)
  #       This is the one present ONLY in the heterogametic sex.
  #   is_xy = TRUE for XY-type systems, FALSE for ZW-type systems.
  #       We need this flag solely to correctly map the user-facing mu
  #       parameter (which is always "proportion from females") onto the
  #       internal m_hom / m_het split (see below).
  is_xy <- xy
  s <- if (xy) n_x else n_z
  l <- if (xy) n_y else n_w

  # Every sex chromosome system must have at least one shared chromosome.
  # (You can have zero limited chromosomes — that's an XO or ZO system.)
  if (s < 1) stop("need at least 1 X (or Z) chromosome")

  # The haploid autosome count must be at least 1 (one pair of autosomes).
  if (any(n_auto < 1)) stop("n_auto must be >= 1")

  # ---------------------------------------------------------------------------
  # CHROMOSOME COUNTS
  # ---------------------------------------------------------------------------

  # Convert haploid autosome count to diploid (the total number of individual
  # autosome chromosomes in a cell). Example: humans have n_auto = 22, so
  # Da = 44 autosome chromosomes per cell.
  Da <- 2 * n_auto

  # Total diploid chromosome number in the HOMOGAMETIC sex.
  # The homogametic sex carries all the autosomes PLUS two copies of the
  # shared sex chromosome (e.g., an XX female has Da + 2 chromosomes).
  Dhom <- Da + 2 * s

  # Total diploid chromosome number in the HETEROGAMETIC sex.
  # The heterogametic sex carries all the autosomes PLUS the shared sex
  # chromosomes PLUS the limited sex chromosomes (e.g., an XY male has
  # Da + 1 + 1 = Da + 2 chromosomes; an XXY male has Da + 2 + 1 = Da + 3).
  Dhet <- Da + s + l

  # ---------------------------------------------------------------------------
  # OCCURRENCE PROBABILITIES WITHIN EACH SEX
  #
  # Following Anderson et al. (2020), we compute the probability that a
  # random fusion between two non-homologous chromosomes falls into each
  # class: autosome-autosome (AA), sex-sex (SS), or sex-autosome (SA).
  #
  # The logic: imagine picking two chromosomes at random from a cell. The
  # first pick is any of the D chromosomes. The second must be non-homologous
  # to the first (it cannot be the other copy of the same chromosome),
  # leaving (D - 2) choices if the first chromosome was one of a homologous
  # pair, or (D - 1) if it was unpaired (hemizygous).
  # ---------------------------------------------------------------------------

  # ---- Homogametic sex ----
  # In the homogametic sex, ALL chromosomes are paired (autosomes in pairs,
  # and the shared sex chromosomes in pairs — e.g., XX). So picking any
  # chromosome always removes itself AND its homolog from the pool.

  # P(AA | homogametic): both chromosomes are autosomes.
  # Numerator: Da ways to pick the first autosome, times (Da - 2) remaining
  #   autosomes that are non-homologous (we subtract 2 because the first
  #   pick and its homolog are removed).
  # Denominator: Dhom total chromosomes for the first pick, times (Dhom - 2)
  #   non-homologous chromosomes for the second pick.
  AA_hom <- Da * (Da - 2) / (Dhom * (Dhom - 2))

  # P(SS | homogametic): both chromosomes are sex chromosomes.
  # There are 2*s sex chromosomes in the homogametic sex (all are the shared
  # type, in s homologous pairs). Same logic as above:
  # Numerator: 2s choices for the first, times (2s - 2) non-homologous sex
  #   chromosomes remaining.
  # Denominator: same Dhom * (Dhom - 2).
  SS_hom <- 2 * s * (2 * s - 2) / (Dhom * (Dhom - 2))

  # P(SA | homogametic): one autosome and one sex chromosome.
  # Since AA + SA + SS = 1, this is just the complement.
  SA_hom <- 1 - AA_hom - SS_hom

  # ---- Heterogametic sex ----
  # In the heterogametic sex, autosomes are still paired, but sex chromosomes
  # are NOT paired with each other — X is not homologous to Y, and if there
  # are multiple X's, each is distinct. This changes the denominator when
  # computing SS probabilities, because picking an unpaired chromosome
  # removes only itself (not a homolog) from the pool.

  # P(AA | heterogametic): both chromosomes are autosomes.
  # Same logic as the homogametic case — autosomes are paired in both sexes.
  AA_het <- Da * (Da - 2) / (Dhet * (Dhet - 2))

  # P(shared-shared | heterogametic): both chromosomes are the shared type
  # (e.g., two X's in an XXY male). This is only nonzero when s >= 2.
  # Since each shared sex chromosome in the heterogametic sex is hemizygous
  # (unpaired — there's no homolog to remove), picking the first shared
  # chromosome removes only that one from the pool, leaving (Dhet - 1)
  # non-homologous chromosomes total. BUT we also need to account for the
  # fact that the Da/2 autosome pairs contribute an additional Da/2
  # "forbidden" homologous matches relative to a naive (Dhet - 1).
  # Anderson et al. (2020) derive:
  #   Denominator = Dhet * (Dhet + s - 1)
  # This accounts for the interaction between paired autosomes and unpaired
  # sex chromosomes. See Anderson et al. Eq 2.2 for the full derivation.
  ss_het <- s * (s - 1) / (Dhet * (Dhet + s - 1))

  # P(limited-limited | heterogametic): both chromosomes are the limited type
  # (e.g., two Y's in an XYY male). Same logic as above but for the limited
  # sex chromosome. Only nonzero when l >= 2.
  ll_het <- l * (l - 1) / (Dhet * (Dhet + l - 1))

  # P(SA | heterogametic): one autosome and one sex chromosome (of either
  # type). Again computed as the complement of AA + SS.
  SA_het <- 1 - AA_het - ss_het - ll_het

  # ---------------------------------------------------------------------------
  # MAP mu (PROPORTION FROM FEMALES) TO HOMOGAMETIC / HETEROGAMETIC
  #
  # The user supplies mu as "proportion of fusions from females."  Internally
  # the math is written in terms of the homogametic sex (m_hom) and the
  # heterogametic sex (m_het).  We need to translate:
  #
  #   XY systems: females ARE the homogametic sex (XX).
  #               m_hom = mu           (female = homogametic)
  #               m_het = 1 - mu       (male   = heterogametic)
  #
  #   ZW systems: females ARE the HETEROGAMETIC sex (ZW).
  #               m_hom = 1 - mu       (male   = homogametic)
  #               m_het = mu           (female = heterogametic)
  #
  # Without this flip, setting mu = 0.7 in a ZW system would incorrectly
  # assign 70% of fusions to the homogametic MALES instead of the intended
  # FEMALES.
  # ---------------------------------------------------------------------------

  if (is_xy) {
    m_hom <- mu         # XY: females are homogametic
    m_het <- 1 - mu     # XY: males are heterogametic
  } else {
    m_hom <- 1 - mu     # ZW: males are homogametic
    m_het <- mu          # ZW: females are heterogametic
  }

  # ---------------------------------------------------------------------------
  # OVERALL OCCURRENCE PROBABILITIES
  #
  # The population-level probability is the weighted average across the two
  # sexes. m_hom is the proportion of fusions originating in the homogametic
  # sex; m_het is the proportion from the heterogametic sex.
  # ---------------------------------------------------------------------------

  P_AA <- m_hom * AA_hom + m_het * AA_het
  P_SA <- m_hom * SA_hom + m_het * SA_het

  # If the user only wants the occurrence model, return P(SA) and stop here.
  if (!weighted) return(P_SA)

  # ---------------------------------------------------------------------------
  # FIXATION WEIGHTS
  #
  # Under neutrality, a new mutation's fixation probability is approximately
  # 1 / (number of copies of that chromosome in the population) (Kimura 1962).
  #
  # For a population of N diploid individuals with equal sex ratio:
  #   - Each autosome exists in 2N copies (N males + N females, 2 each)
  #   - Each shared sex chrom (X or Z) exists in 3N/2 copies
  #     (2 per homogametic individual + 1 per heterogametic = 3 per pair)
  #   - Each limited sex chrom (Y or W) exists in N/2 copies
  #     (1 per heterogametic individual only)
  #
  # The fixation probability is proportional to the INVERSE of these counts.
  # We express weights RELATIVE to the autosomal weight (set to 1):
  #   wA = 1          (baseline: autosome, 2N copies)
  #   wX = 4/3 ≈ 1.33 (shared sex chrom: 2N / (3N/2) = 4/3 times more
  #                      likely to fix than an autosomal mutation)
  #   wY = 4          (limited sex chrom: 2N / (N/2) = 4 times more likely
  #                      to fix than an autosomal mutation)
  #
  # N cancels out when we compute proportions, so these weights are universal.
  # ---------------------------------------------------------------------------

  wA <- 1      # autosomal fixation weight (baseline)
  wX <- 4 / 3  # shared sex chromosome fixation weight
  wY <- 4      # limited sex chromosome fixation weight

  # ---------------------------------------------------------------------------
  # SPLIT SA AND SS FUSIONS BY SEX CHROMOSOME TYPE
  #
  # To apply the fixation weights, we need to know what fraction of SA fusions
  # involve the shared vs. limited sex chromosome, because they have different
  # fixation weights. Similarly, SS fusions can be shared-shared (e.g., X-X)
  # or limited-limited (e.g., Y-Y), and these also get different weights.
  # (Note: mixed SS fusions like X-Y are already counted in SA_het via the
  # complement — they aren't a separate class in the Anderson framework.)
  # ---------------------------------------------------------------------------

  # Total number of sex chromosomes in the heterogametic sex
  stot <- s + l

  # P(shared-autosome): XA or ZA fusions.
  # In the homogametic sex, ALL sex chromosomes are the shared type, so every
  # SA fusion from that sex is a shared-autosome fusion.
  # In the heterogametic sex, SA fusions are split between shared-A and
  # limited-A in proportion to the number of each sex chromosome type present.
  P_sA <- m_hom * SA_hom + m_het * SA_het * s / stot

  # P(limited-autosome): YA or WA fusions.
  # These can ONLY arise in the heterogametic sex (the homogametic sex has
  # no limited sex chromosomes). The fraction is proportional to l / stot.
  P_lA <- m_het * SA_het * l / stot

  # P(shared-shared): e.g., X-X fusions across the population.
  # From the homogametic sex (where all sex chroms are shared) plus the
  # shared-shared component of the heterogametic sex.
  P_ss <- m_hom * SS_hom + m_het * ss_het

  # P(limited-limited): e.g., Y-Y fusions. Only from the heterogametic sex.
  P_ll <- m_het * ll_het

  # ---------------------------------------------------------------------------
  # FIXATION-WEIGHTED PROPORTION OF SA FUSIONS
  #
  # We multiply each fusion class by its fixation weight, then renormalize so
  # the weighted proportions sum to 1. The fixation-weighted P(SA) is:
  #
  #            wX * P(shared-A) + wY * P(limited-A)
  #   P_fix = -----------------------------------------------
  #           wA*P(AA) + wX*[P(sA)+P(ss)] + wY*[P(lA)+P(ll)]
  #
  # The numerator is the fixation-weighted SA component.
  # The denominator is the fixation-weighted total across all fusion classes.
  # ---------------------------------------------------------------------------

  # Numerator: weighted sum of the two SA sub-types
  num <- wX * P_sA + wY * P_lA

  # Denominator: weighted sum of ALL fusion classes
  # AA fusions get weight wA (autosomal).
  # Shared sex-chrom fusions (both SA and SS involving X/Z) get weight wX.
  # Limited sex-chrom fusions (both SA and SS involving Y/W) get weight wY.
  den <- wA * P_AA + wX * (P_sA + P_ss) + wY * (P_lA + P_ll)

  # Return the fixation-weighted proportion of SA fusions
  num / den
}
