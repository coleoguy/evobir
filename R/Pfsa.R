#' Proportion of Sex-Autosome Fusions
#'
#' Calculates the expected proportion of chromosomal fusions that join a sex
#' chromosome and an autosome. Two modes are available: \code{"fixation"}
#' weights each fusion class by its neutral fixation probability (accounting
#' for the different effective copy numbers of autosomes, X, and Y
#' chromosomes), while \code{"occurrence"} returns the raw occurrence
#' probability from Anderson et al. (2020).
#'
#' @param Da Integer. The diploid autosome count (must be >= 4).
#' @param scs Character string describing the sex chromosome system in the
#'   heterogametic sex. Built from any combination of X's and Y's (XY-type)
#'   or Z's and W's (ZW-type). Use \code{"O"} for the absent chromosome in
#'   X0/Z0 systems. Common examples: \code{"XY"}, \code{"XXY"}, \code{"XYY"},
#'   \code{"XXXXY"}, \code{"ZW"}, \code{"ZZW"}, \code{"XO"}, \code{"ZO"}.
#'   The string must contain at least one sex chromosome and only characters
#'   from one system (do not mix X/Y with Z/W).
#' @param mud Numeric between 0 and 1. The proportion of fusions that
#'   originate in the homogametic sex (females in XY systems, males in ZW
#'   systems). Default is 0.5 (equal contribution from both sexes).
#' @param type Character. Either \code{"fixation"} (default) to return the
#'   fixation-weighted proportion, or \code{"occurrence"} to return the raw
#'   occurrence probability from Anderson et al. (2020).
#' @param wA Numeric. Relative neutral fixation weight for autosome fusions.
#'   Default is 1. Only used when \code{type = "fixation"}.
#' @param wX Numeric. Relative neutral fixation weight for X-linked (or
#'   Z-linked) fusions. Default is 4/3. Only used when
#'   \code{type = "fixation"}.
#' @param wY Numeric. Relative neutral fixation weight for Y-linked (or
#'   W-linked) fusions. Default is 4. Only used when
#'   \code{type = "fixation"}.
#'
#' @return A named list with class \code{"PropSA"} containing:
#'   \describe{
#'     \item{P_SA}{The proportion of SA fusions (fixation-weighted or
#'       occurrence, depending on \code{type}).}
#'     \item{P_AA}{Proportion of autosome-autosome fusions.}
#'     \item{P_SS}{Proportion of sex-sex fusions.}
#'     \item{type}{Which mode was used.}
#'     \item{scs}{The sex chromosome system.}
#'     \item{Da}{The diploid autosome count.}
#'   }
#'
#' @details
#' The occurrence model follows Anderson et al. (2020, Biology Letters).
#' Under a null model where any chromosome is equally likely to fuse with
#' any other non-homologous chromosome, the probability of each fusion type
#' depends on the chromosome counts in each sex, weighted by the proportion
#' of fusions originating in each sex (\code{mud}).
#'
#' The fixation model extends this by recognizing that a new fusion's
#' probability of reaching fixation under neutrality is approximately
#' 1/(effective number of copies). Autosomes exist in 2N copies, X
#' chromosomes in ~3N/2 copies (two in females, one in males), and Y
#' chromosomes in ~N/2 copies (one in males only). The fixation weights
#' (\code{wA}, \code{wX}, \code{wY}) are the inverses of these copy
#' numbers, so Y-linked fusions fix more readily than autosomal ones.
#'
#' @references
#' Anderson, N.W., Hjelmen, C.E., & Blackmon, H. (2020). The probability
#' of fusions joining sex chromosomes and autosomes. \emph{Biology Letters},
#' 16, 20200648.
#'
#' @examples
#' # Standard XY system with 20 autosomes, fixation-weighted (default)
#' Prop_SA(Da = 20, scs = "XY")
#'
#' # Same system, occurrence probability only
#' Prop_SA(Da = 20, scs = "XY", type = "occurrence")
#'
#' # ZW system with female-biased fusion origin
#' Prop_SA(Da = 10, scs = "ZW", mud = 0.7)
#'
#' # Multi-X system
#' Prop_SA(Da = 26, scs = "XXO", type = "fixation")
#'
#' # Arbitrary multi-sex-chromosome system
#' Prop_SA(Da = 20, scs = "XXXXY")
#'
#' @export
Prop_SA <- function(Da,
                    scs,
                    mud  = 0.5,
                    type = "fixation",
                    wA   = 1,
                    wX   = 4 / 3,
                    wY   = 4) {

  # ---- Input validation ----
  # Make sure the user gave us sensible inputs before we do any math

  # type must be one of our two options
  type <- match.arg(type, choices = c("fixation", "occurrence"))

  # Da must be a positive even-ish integer >= 4 so denominators don't blow up
  if (!is.numeric(Da) || Da < 4) {
    stop("Da (diploid autosome count) must be a number >= 4.")
  }

  # mud is the proportion of fusions from the homogametic sex: must be 0-1
  if (!is.numeric(mud) || mud < 0 || mud > 1) {
    stop("mud must be a number between 0 and 1.")
  }

  # ---- Parse the sex chromosome string ----
  # Split "XXY" into c("X", "X", "Y") so we can count each type.
  # The string can be ANY length, the math is fully general.
  sexchroms <- strsplit(scs, split = "")[[1]]

  # Figure out which system we're in by looking at what letters are present.
  # XY-type systems contain X's and optionally Y's (plus O for XO systems).
  # ZW-type systems contain Z's and optionally W's (plus O for ZO systems).
  # We don't allow mixing (e.g., "XZ" makes no biological sense).
  has_X <- "X" %in% sexchroms
  has_Y <- "Y" %in% sexchroms
  has_Z <- "Z" %in% sexchroms
  has_W <- "W" %in% sexchroms
  has_O <- "O" %in% sexchroms

  # Check that only one system type is represented
  is_xy_type <- has_X || has_Y
  is_zw_type <- has_Z || has_W

  if (is_xy_type && is_zw_type) {
    stop("scs mixes XY and ZW characters. Use only X/Y/O or Z/W/O.")
  }
  if (!is_xy_type && !is_zw_type) {
    stop("scs must contain at least one sex chromosome (X, Y, Z, or W).")
  }

  # Make sure there are no unexpected characters in the string
  allowed <- if (is_xy_type) c("X", "Y", "O") else c("Z", "W", "O")
  bad_chars <- setdiff(sexchroms, allowed)
  if (length(bad_chars) > 0) {
    stop("scs contains invalid characters: ",
         paste(bad_chars, collapse = ", "),
         ". Only ", paste(allowed, collapse = "/"), " are allowed.")
  }

  # We use a unified naming convention internally:
  #   n1 = count of the "shared" sex chromosome (X in XY, Z in ZW)
  #        This is the one present in BOTH sexes (homogametic sex has 2*n1)
  #   n2 = count of the "limited" sex chromosome (Y in XY, W in ZW)
  #        This is the one present ONLY in the heterogametic sex
  #   is_xy = TRUE for XY-type, FALSE for ZW-type
  #   "O" is just a placeholder for the absent chromosome it adds 0
  if (is_xy_type) {
    is_xy <- TRUE
    n1 <- sum(sexchroms == "X")  # number of X chromosomes in males
    n2 <- sum(sexchroms == "Y")  # number of Y chromosomes in males
  } else {
    is_xy <- FALSE
    n1 <- sum(sexchroms == "Z")  # number of Z chromosomes in females
    n2 <- sum(sexchroms == "W")  # number of W chromosomes in females
  }

  # Need at least one sex chromosome (n1 must be >= 1)
  if (n1 < 1) {
    stop("scs must contain at least one X (for XY-type) or Z (for ZW-type).")
  }

  # ---- Compute chromosome totals for each sex ----
  #
  # In XY systems:
  #   Homogametic sex (female) has Da autosomes + 2*Xs X chromosomes
  #   Heterogametic sex (male) has Da autosomes + Xs X + Y Y chromosomes
  #
  # In ZW systems it's flipped:
  #   Homogametic sex (male) has Da autosomes + 2*Zd Z chromosomes
  #   Heterogametic sex (female) has Da autosomes + Zd Z + W W chromosomes
  #
  # We call these D_hom (homogametic diploid number) and
  #                D_het (heterogametic diploid number)

  D_hom <- Da + 2 * n1             # e.g., for XY: Da + 2*Xs (the XX sex)
  D_het <- Da + n1 + n2            # e.g., for XY: Da + Xs + Y (the XY sex)

  # ---- Proportion of fusions from each sex ----
  #
  # mud = proportion from the homogametic sex (dam in XY, sire in ZW)
  # For ZW systems the paper swaps the roles: the "dam" is heterogametic,
  # so we flip mud to get the proportion from each sex correctly.

  if (is_xy) {
    # XY: mud = fraction from females (homogametic)
    m_hom <- mud           # fraction of fusions from the homogametic sex
    m_het <- 1 - mud       # fraction from the heterogametic sex
  } else {
    # ZW: mud = fraction from females, but females are HETEROGAMETIC
    # so the homogametic sex (males) contributes (1 - mud) of fusions
    m_hom <- 1 - mud       # fraction from homogametic sex (males in ZW)
    m_het <- mud           # fraction from heterogametic sex (females in ZW)
  }

  # ---- Helper: safe denominator ----
  # Avoid division by zero for edge cases (e.g., very small chromosome counts)
  den_safe <- function(x) {
    if (abs(x) < 1e-12) return(NA_real_)
    return(x)
  }

  # ====================================================================
  # OCCURRENCE PROBABILITIES
  # From Anderson et al. (2020), Equations 2.2 - 2.4
  #
  # The key idea: when two chromosomes fuse, we're picking 2 chromosome
  # arms. The model excludes fusions between homologs (which would make
  # unbalanced gametes).
  #
  # For the HOMOGAMETIC sex (e.g., XX females in XY systems):
  #   - All chromosomes are paired, so picking one arm means the other
  #     arm of the same chromosome is excluded. That's why we see
  #     Da*(Da-2) for autosomes: Da arms total, minus the chosen arm
  #     and its partner on the same chromosome = (Da-2) options.
  #   - Same logic for sex chromosomes: 2*n1 total sex chrom arms,
  #     and 2*n1*(2*n1-2) ordered pairs excluding same-chromosome.
  #
  # For the HETEROGAMETIC sex (e.g., XY males):
  #   - Autosomes are still paired, so Da*(Da-2) / [D_het*(D_het-2)]
  #   - But X and Y chromosomes are NOT paired with each other.
  #     X chromosomes pair with other X's (if multiple), and Y's with Y's.
  #     So each type gets its own term with a type-specific denominator:
  #       X-X fusions: n1*(n1-1) / [D_het*(D_het + n1 - 1)]
  #       Y-Y fusions: n2*(n2-1) / [D_het*(D_het + n2 - 1)]
  #     (See Anderson et al. Eq 2.2 for the derivation of these
  #      denominators.)
  # ====================================================================

  # -- P(AA): probability both fused chromosomes are autosomes --
  #
  # Homogametic sex: Da autosome arms, D_hom total arms
  P_AA_hom <- (Da * (Da - 2)) / den_safe(D_hom * (D_hom - 2))

  # Heterogametic sex: same numerator, different total
  P_AA_het <- (Da * (Da - 2)) / den_safe(D_het * (D_het - 2))

  # Weighted average across the two sexes
  P_AA <- m_hom * P_AA_hom + m_het * P_AA_het

  # -- P(SS): probability both fused chromosomes are sex chromosomes --
  #
  # Homogametic sex: 2*n1 sex chromosome arms, all paired
  # Number of ordered non-homolog pairs = 2*n1 * (2*n1 - 2)
  P_SS_hom <- (2 * n1 * (2 * n1 - 2)) / den_safe(D_hom * (D_hom - 2))

  # Heterogametic sex: separate terms for each sex chromosome type
  # n1-type (X or Z) fusing with another n1-type:
  P_n1n1_het <- (n1 * (n1 - 1)) / den_safe(D_het * (D_het + n1 - 1))

  # n2-type (Y or W) fusing with another n2-type:
  P_n2n2_het <- (n2 * (n2 - 1)) / den_safe(D_het * (D_het + n2 - 1))

  # Total SS in heterogametic sex is the sum of both types
  P_SS_het <- P_n1n1_het + P_n2n2_het

  # Weighted average
  P_SS <- m_hom * P_SS_hom + m_het * P_SS_het

  # -- P(SA): everything that isn't AA or SS --
  P_SA <- 1 - P_AA - P_SS

  # ====================================================================
  # If the user only wants occurrence probabilities, we're done!
  # ====================================================================
  if (type == "occurrence") {
    result <- list(
      P_SA  = P_SA,
      P_AA  = P_AA,
      P_SS  = P_SS,
      type  = type,
      scs   = scs,
      Da    = Da
    )
    class(result) <- "PropSA"
    return(result)
  }

  # ====================================================================
  # FIXATION-WEIGHTED PROPORTIONS
  #
  # Under neutrality, a new mutation's fixation probability is roughly
  # 1 / (number of copies in the population). For a population of N
  # diploid individuals with equal sex ratio:
  #   - Each autosome exists in 2N copies  -> weight wA ~ 1/2N
  #   - Each X exists in 3N/2 copies       -> weight wX ~ 2/3N
  #   - Each Y exists in N/2 copies        -> weight wY ~ 2/N
  #
  # We only care about RELATIVE weights (N cancels out in proportions),
  # so the defaults are wA=1, wX=4/3, wY=4.
  #
  # To apply these weights, we need to split P(SA) into the part
  # involving n1-type sex chromosomes (X or Z) and the part involving
  # n2-type (Y or W), because they have different fixation weights.
  # ====================================================================

  # -- Split P(SA) into n1-Autosome and n2-Autosome components --
  #
  # In the homogametic sex, all sex chromosomes are n1-type (XX or ZZ),
  # so ALL SA fusions from that sex are n1-A fusions.
  P_SA_hom <- 1 - P_AA_hom - P_SS_hom   # total SA from homogametic sex
  P_SA_het <- 1 - P_AA_het - P_SS_het    # total SA from heterogametic sex

  # In the heterogametic sex, SA fusions could involve either n1 or n2
  # chromosomes. We split in proportion to how many of each type are
  # present in that sex.
  n_total_het <- n1 + n2  # total sex chromosomes in heterogametic sex

  if (n_total_het > 0) {
    # n1-A fusions: all from homogametic + proportional share from heterogametic
    P_n1A <- m_hom * P_SA_hom +
      m_het * P_SA_het * (n1 / n_total_het)

    # n2-A fusions: only from heterogametic sex
    P_n2A <- m_het * P_SA_het * (n2 / n_total_het)
  } else {
    # Edge case: no sex chromosomes at all (shouldn't happen, but be safe)
    P_n1A <- m_hom * P_SA_hom
    P_n2A <- 0
  }

  # -- Apply fixation weights and renormalize --
  #
  # Each class of fusion gets multiplied by its weight, then we
  # divide by the total to get the proportion among FIXED fusions.

  # Numerator: fixation-weighted SA fusions
  fix_num <- wX * P_n1A + wY * P_n2A

  # Denominator: fixation-weighted total of all fusion types
  # Note: SS fusions also split into n1-n1 and n2-n2 for weighting.
  # For n1-n1 (e.g., X-X) fusions, both chromosomes are n1-type -> use wX
  # For n2-n2 (e.g., Y-Y) fusions, both chromosomes are n2-type -> use wY
  P_n1n1 <- m_hom * P_SS_hom + m_het * P_n1n1_het
  P_n2n2 <- m_het * P_n2n2_het

  fix_den <- wA * P_AA + wX * (P_n1A + P_n1n1) + wY * (P_n2A + P_n2n2)

  # The fixation-weighted proportion of SA fusions
  P_SA_fix <- fix_num / den_safe(fix_den)

  # Also compute fixation-weighted AA and SS for completeness
  P_AA_fix <- (wA * P_AA) / den_safe(fix_den)
  P_SS_fix <- (wX * P_n1n1 + wY * P_n2n2) / den_safe(fix_den)

  result <- list(
    P_SA  = P_SA_fix,
    P_AA  = P_AA_fix,
    P_SS  = P_SS_fix,
    type  = type,
    scs   = scs,
    Da    = Da
  )
  class(result) <- "PropSA"
  return(result)
}


#' Print method for PropSA objects
#'
#' @param x A \code{PropSA} object returned by \code{\link{Prop_SA}}.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.PropSA <- function(x, ...) {
  cat(sprintf("Sex chromosome system: %s | Diploid autosomes: %d\n",
              x$scs, x$Da))
  cat(sprintf("Mode: %s\n\n", x$type))
  cat(sprintf("  P(SA) = %.5f\n", x$P_SA))
  cat(sprintf("  P(AA) = %.5f\n", x$P_AA))
  cat(sprintf("  P(SS) = %.5f\n", x$P_SS))
  invisible(x)
}
