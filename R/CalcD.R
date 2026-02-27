#' Statistical Mode
#'
#' Returns the most frequently occurring value in a vector.  R does not have
#' a built-in mode function for data (its \code{mode()} returns storage type),
#' so this simple helper fills that gap.
#'
#' @param x A vector of any type.
#'
#' @return The most common value in \code{x}.  If there is a tie, the first
#'   value (in order of appearance) wins.
#'
#' @keywords internal
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


#' Patterson's D-Statistic (ABBA-BABA Test)
#'
#' Calculates Patterson's D-statistic to test for introgression (gene flow)
#' between species using genomic sequence data.  The test looks at two
#' patterns of shared alleles -- called "ABBA" and "BABA" -- that should
#' be equally common if there is no introgression.  A significant departure
#' from D = 0 suggests gene flow between non-sister taxa.
#'
#' The input alignment must contain exactly 4 sequences in this order:
#' \enumerate{
#'   \item \strong{P1} -- first ingroup population
#'   \item \strong{P2} -- second ingroup population (suspected introgression recipient)
#'   \item \strong{P3} -- third population (suspected introgression donor)
#'   \item \strong{O}  -- outgroup
#' }
#'
#' @param alignment Character string.  Path to the alignment file
#'   (default \code{"alignment.fasta"}).
#' @param sig.test Character.  Type of significance test:
#'   \itemize{
#'     \item \code{"N"} (default) -- no significance test, just calculate D.
#'     \item \code{"B"} -- bootstrap: resample sites with replacement.
#'     \item \code{"J"} -- jackknife: drop blocks of sites.  Good when
#'       nearby SNPs may be in linkage disequilibrium.
#'   }
#' @param ambig Character.  How to handle ambiguous bases (e.g., R, Y, S):
#'   \itemize{
#'     \item \code{"D"} (default) -- drop sites with any ambiguity.
#'     \item \code{"R"} -- randomly resolve each ambiguous base.
#'     \item \code{"I"} -- ignore ambiguities (keep them as-is).
#'   }
#' @param block.size Integer.  Number of contiguous sites to drop in each
#'   jackknife replicate (default 1000).  Only used when \code{sig.test = "J"}.
#' @param replicate Integer.  Number of bootstrap or jackknife replicates
#'   (default 1000).
#' @param align.format Character.  Format of the alignment file
#'   (default \code{"fasta"}).  Any format accepted by
#'   \code{\link[seqinr]{read.alignment}}.
#'
#' @return The D-statistic value (numeric).  Also prints summary statistics
#'   to the console including site counts and (if requested) p-values.
#'
#' @details
#' D is calculated as (ABBA - BABA) / (ABBA + BABA).
#' \itemize{
#'   \item D = 0 means no evidence of introgression.
#'   \item D > 0 suggests gene flow between P2 and P3.
#'   \item D < 0 suggests gene flow between P1 and P3.
#' }
#' Only biallelic sites are used.  Significance is assessed by computing
#' a Z-score (|D| / SD) from the resampled distribution, then converting
#' to a two-tailed p-value.
#'
#' @references
#' Durand, E.Y., Patterson, N., Reich, D., & Slatkin, M. (2011). Testing
#' for ancient admixture between closely related populations.
#' \emph{Molecular Biology and Evolution}, 28(8), 2239-2252.
#'
#' @examples
#' \dontrun{
#' # Basic D-statistic with no significance test
#' CalcD("my_alignment.fasta")
#'
#' # With bootstrap significance test
#' CalcD("my_alignment.fasta", sig.test = "B", replicate = 500)
#'
#' # With jackknife (good for whole-genome data with LD)
#' CalcD("my_alignment.fasta", sig.test = "J", block.size = 5000)
#' }
#'
#' @importFrom seqinr read.alignment
#' @importFrom stats sd pnorm lm
#' @export
CalcD <- function(alignment = "alignment.fasta",
                  sig.test = "N",
                  ambig = "D",
                  block.size = 1000,
                  replicate = 1000,
                  align.format = 'fasta') {

  # ---- Inner function: core D calculation (vectorized) ----
  # Takes a 4-row alignment matrix and counts ABBA/BABA patterns.
  # All comparisons are done column-wise with no R-level for-loop.
  d.calc <- function(alignment) {
    # Extract each taxon's row as a vector for fast element-wise comparison
    p1 <- alignment[1, ]
    p2 <- alignment[2, ]
    p3 <- alignment[3, ]
    o  <- alignment[4, ]

    # Biallelic filter: exactly 2 unique alleles at the site.
    # A site is biallelic when not all four are the same AND at most
    # two distinct values exist.  Fastest check: all pairwise combos.
    n_alleles <- (p1 != p2) | (p1 != p3) | (p1 != o)   # at least 2
    too_many  <- (p1 != p2) & (p1 != p3) & (p2 != p3)  # proxy for 3+
    # Refine: sites where 3+ alleles exist among P1,P2,P3 or outgroup
    # differs from all three ingroup.  Full exact check:
    biallelic <- n_alleles & !( (p1 != p2) & (p1 != p3) & (p2 != p3) &
                                 !((o == p1) | (o == p2) | (o == p3)) )
    # Simpler exact method: count unique alleles per column
    biallelic <- vapply(seq_len(ncol(alignment)), function(i) {
      length(unique(alignment[, i])) == 2L
    }, logical(1))

    # Informative sites: biallelic, P1 != P2, and outgroup != P3
    informative <- biallelic & (p1 != p2) & (o != p3)

    # ABBA: P3 matches P2 (not P1) at informative sites
    abba <- sum(informative & (p2 == p3))

    # BABA: P3 matches P1 (not P2) at informative sites
    baba <- sum(informative & (p3 == p1))

    d <- (abba - baba) / (abba + baba)
    return(list(d, abba, baba))
  }

  # ---- Read and parse the alignment ----
  alignment <- read.alignment(alignment, format = align.format,
                              forceToLower = TRUE)
  # Vectorized parsing: split all sequences at once, bind into matrix
  alignment.matrix <- do.call(rbind, strsplit(unlist(alignment$seq), ""))

  # ---- Handle ambiguous bases ----
  # IUPAC ambiguity codes: R=A/G, Y=C/T, S=G/C, W=A/T, K=G/T, M=A/C
  if (ambig == "D") {
    # Drop: remove any site containing a non-standard base
    target <- c("a", "c", "g", "t")
    keep <- vapply(1:ncol(alignment.matrix), function(i) {
      all(alignment.matrix[, i] %in% target)
    }, logical(1))
    alignment.matrix <- alignment.matrix[, keep]
  }

  if (ambig == "R") {
    # Random resolution: randomly pick one of the two possible bases
    target <- c("a", "c", "g", "t", "r", "y", "s", "w", "k", "m")
    keep <- vapply(1:ncol(alignment.matrix), function(i) {
      all(alignment.matrix[, i] %in% target)
    }, logical(1))
    alignment.matrix <- alignment.matrix[, keep]

    # Resolve each ambiguous base to one of its two possibilities
    resolver <- function(x) {
      switch(x,
             "r" = sample(c("a", "g"), 1),
             "y" = sample(c("c", "t"), 1),
             "s" = sample(c("g", "c"), 1),
             "w" = sample(c("a", "t"), 1),
             "k" = sample(c("g", "t"), 1),
             "m" = sample(c("a", "c"), 1),
             x)  # if it's already a/c/g/t, keep it
    }
    alignment.matrix <- apply(alignment.matrix, c(1, 2), resolver)
  }

  # ---- Calculate D on the real data ----
  results <- d.calc(alignment.matrix)
  d <- results[[1]]
  if (is.nan(d)) d <- 0  # no informative sites -> D is zero
  abba <- results[[2]]
  baba <- results[[3]]

  # ---- Bootstrap significance test ----
  # Resample columns (sites) WITH replacement to make new datasets
  # of equal size, then calculate D for each replicate
  if (sig.test == "B") {
    sim.d <- numeric(replicate)  # pre-allocate for speed
    foo <- ncol(alignment.matrix)
    cat("\nperforming bootstrap")
    for (k in 1:replicate) {
      if (k %% 100 == 0) cat(".")
      sim.matrix <- alignment.matrix[1:4, sample.int(foo, foo, replace = TRUE)]
      sim.d[k] <- d.calc(sim.matrix)[[1]]
    }
    sim.d[is.nan(sim.d)] <- 0
    # Z-score: how many SDs away from zero is our observed D?
    z <- abs(d / sd(sim.d))
    new.pval <- 2 * (1 - pnorm(z))  # two-tailed p-value

    cat("\nSites in alignment =", ncol(alignment.matrix))
    cat("\nNumber of sites with ABBA pattern =", abba)
    cat("\nNumber of sites with BABA pattern =", baba)
    cat("\n\nD raw statistic / Z-score = ", d, " / ", z)
    cat("\n\nResults from ", replicate, "bootstraps")
    cat("\nSD D statistic =", sd(sim.d))
    cat("\nP-value (that D=0) = ", new.pval, "\n\n")
  }

  # ---- Jackknife significance test ----
  # Drop a contiguous block of sites and recalculate D
  # Better than bootstrap when SNPs are in linkage disequilibrium
  if (sig.test == "J") {
    max.rep <- ncol(alignment.matrix) - block.size
    if (block.size >= (ncol(alignment.matrix) / 2)) {
      stop(call. = FALSE,
           paste("\nWith a block size of", block.size,
                 "and an alignment of", ncol(alignment.matrix),
                 "sites \nsome sites would never be included in the \nanalysis",
                 "\n\nThe maximum block size is 1/2 the alignment length"))
    }
    if (max.rep < replicate) {
      stop(call. = FALSE,
           paste("\nWith a block size of", block.size,
                 "and an alignment of", ncol(alignment.matrix),
                 "sites", replicate, "replicates\nare not possible"))
    }

    # Evenly space the drop positions across the alignment
    drop.pos <- seq.int(from = 1, to = (max.rep - 1), length.out = replicate)
    replicate2 <- replicate

    sim.d <- numeric(replicate2)  # pre-allocate for speed
    foo <- ncol(alignment.matrix)
    cat("\nperforming jackknife")
    for (k in 1:replicate2) {
      if (k / 2 == round(k / 2)) cat(".")
      drop_range <- drop.pos[k]:(drop.pos[k] + block.size)
      sim.matrix <- alignment.matrix[1:4, -drop_range]
      sim.d[k] <- d.calc(sim.matrix)[[1]]
    }
    sim.d[is.nan(sim.d)] <- 0
    z <- abs(d / sd(sim.d))
    new.pval <- 2 * (1 - pnorm(z))

    cat("\nSites in alignment =", ncol(alignment.matrix))
    cat("\nNumber of sites with ABBA pattern =", abba)
    cat("\nNumber of sites with BABA pattern =", baba)
    cat("\nD raw statistic", d)
    cat("\nZ-score = ", z)
    cat("\n\nResults from", replicate2, "jackknifes with block size of", block.size)
    cat("\nSD D statistic =", sd(sim.d))
    cat("\nP-value (that D=0) = ", new.pval, "\n\n")
  }

  # ---- No significance test: just print the raw D ----
  if (sig.test == "N") {
    cat("\nSites in alignment =", ncol(alignment.matrix))
    cat("\nNumber of sites with ABBA pattern =", abba)
    cat("\nNumber of sites with BABA pattern =", baba)
    cat("\n\nD raw statistic = ", d, "\n\n")
  }

  return(d)
}


#' Population-Level Patterson's D-Statistic
#'
#' Like \code{\link{CalcD}}, but designed for population-level data where
#' you have multiple individuals per population.  Instead of counting
#' discrete ABBA/BABA patterns, it uses allele frequencies to compute
#' a frequency-based D-statistic (Equation 2 from Durand et al. 2011).
#'
#' \strong{When should you use this instead of CalcD?}  Use
#' \code{CalcPopD} when you have sequenced multiple individuals from
#' each population (e.g., 10 humans, 10 chimps, 10 gorillas, 10
#' orangutans).  The frequency-based approach extracts more information
#' from your data because it considers how common each allele is within
#' each population, rather than just looking at a single representative
#' sequence.
#'
#' The alignment must contain sequences grouped by population name.
#' Sequences from the same population should share the same name in the
#' FASTA header.  Populations must appear in this order: P1, P2, P3,
#' Outgroup.
#'
#' @inheritParams CalcD
#'
#' @return A list with:
#'   \describe{
#'     \item{d.stat}{The D-statistic value.}
#'     \item{pval}{P-value (only if \code{sig.test != "N"}).}
#'     \item{align.length}{Total number of sites in the alignment.}
#'     \item{useful.sites}{Number of biallelic sites with ABBA/BABA signal.}
#'     \item{seg.sites}{Number of useful sites still segregating within
#'       at least one population.}
#'   }
#'
#' @details
#' For each biallelic site, the function calculates the frequency of the
#' derived allele (B) in each population:
#' \deqn{D = \sum[(1-\hat{p}_1)\hat{p}_2\hat{p}_3(1-\hat{p}_4) -
#'   \hat{p}_1(1-\hat{p}_2)\hat{p}_3(1-\hat{p}_4)] /
#'   \sum[(1-\hat{p}_1)\hat{p}_2\hat{p}_3(1-\hat{p}_4) +
#'   \hat{p}_1(1-\hat{p}_2)\hat{p}_3(1-\hat{p}_4)]}
#'
#' This frequency-based approach is more powerful than the single-sequence
#' version when population samples are available because it uses all the
#' information about allele frequencies rather than just presence/absence.
#'
#' @references
#' Durand, E.Y., Patterson, N., Reich, D., & Slatkin, M. (2011). Testing
#' for ancient admixture between closely related populations.
#' \emph{Molecular Biology and Evolution}, 28(8), 2239-2252.
#'
#' @examples
#' \dontrun{
#' # Population-level D-statistic with bootstrap
#' result <- CalcPopD("population_alignment.fasta", sig.test = "B")
#' result$d.stat
#' result$pval
#' }
#'
#' @export
CalcPopD <- function(alignment = "alignment.fasta",
                     sig.test = "N",
                     ambig = "D",
                     block.size = 1000,
                     replicate = 1000,
                     align.format = 'fasta') {

  # ---- Read alignment ----
  # The first column will hold population names, the rest hold bases
  alignment <- read.alignment(alignment, format = align.format,
                              forceToLower = TRUE)
  # Vectorized parsing: split all sequences at once, bind into matrix
  seq.matrix <- do.call(rbind, strsplit(unlist(alignment$seq), ""))
  # Prepend population names as column 1
  alignment.matrix <- cbind(unlist(alignment$nam), seq.matrix)

  # ---- Handle ambiguous bases ----
  if (ambig == "D") {
    target <- c("a", "c", "g", "t")
    # Vectorized column filtering (skip column 1 = names)
    keep <- c(TRUE, vapply(2:ncol(alignment.matrix), function(i) {
      all(alignment.matrix[, i] %in% target)
    }, logical(1)))
    alignment.matrix <- alignment.matrix[, keep]
  }

  if (ambig == "R") {
    target <- c("a", "c", "g", "t", "r", "y", "s", "w", "k", "m")
    keep <- c(TRUE, vapply(2:ncol(alignment.matrix), function(i) {
      all(alignment.matrix[, i] %in% target)
    }, logical(1)))
    alignment.matrix <- alignment.matrix[, keep]

    resolver <- function(x) {
      switch(x,
             "r" = sample(c("a", "g"), 1),
             "y" = sample(c("c", "t"), 1),
             "s" = sample(c("g", "c"), 1),
             "w" = sample(c("a", "t"), 1),
             "k" = sample(c("g", "t"), 1),
             "m" = sample(c("a", "c"), 1),
             x)
    }
    temp.mat <- alignment.matrix[, -1]
    alignment.matrix[, 2:ncol(alignment.matrix)] <-
      apply(temp.mat, c(1, 2), resolver)
  }

  # Get the 4 unique population names (in order: P1, P2, P3, Outgroup)
  groups <- unique(alignment$nam)

  # ---- Inner function: frequency-based D calculation ----
  dpop.calc <- function(alignment.matrix) {
    numerator <- 0
    denominator <- 0
    useful <- 0      # sites with ABBA/BABA signal
    segregating <- 0 # useful sites still polymorphic in >= 1 population

    for (i in 2:ncol(alignment.matrix)) {
      seg.pos <- FALSE

      # Only biallelic sites
      if (length(unique(alignment.matrix[, i])) == 2) {
        # "A" = most common allele in outgroup, "B" = the other
        A <- Mode(alignment.matrix[alignment.matrix[, 1] == groups[4], i])
        B <- unique(alignment.matrix[, i])[unique(alignment.matrix[, i]) != A]

        # P3 must carry allele B
        if (B %in% unique(alignment.matrix[alignment.matrix[, 1] == groups[3], i])) {
          # P1 and P2 must differ
          if (length(unique(alignment.matrix[alignment.matrix[, 1] %in% groups[1:2], i])) == 2) {
            useful <- useful + 1

            # Check if site is still segregating within any population
            if (length(unique(alignment.matrix[alignment.matrix[, 1] == groups[1], i])) == 2) seg.pos <- TRUE
            if (length(unique(alignment.matrix[alignment.matrix[, 1] == groups[2], i])) == 2) seg.pos <- TRUE
            if (length(unique(alignment.matrix[alignment.matrix[, 1] == groups[3], i])) == 2) seg.pos <- TRUE
            if (length(unique(alignment.matrix[alignment.matrix[, 1] == groups[4], i])) == 2) seg.pos <- TRUE
            if (seg.pos) segregating <- segregating + 1

            # Calculate allele frequencies (proportion of allele B)
            p1 <- 1 - (sum(alignment.matrix[alignment.matrix[, 1] == groups[1], i] == A) /
                          length(alignment.matrix[alignment.matrix[, 1] == groups[1], i]))
            p2 <- 1 - (sum(alignment.matrix[alignment.matrix[, 1] == groups[2], i] == A) /
                          length(alignment.matrix[alignment.matrix[, 1] == groups[2], i]))
            p3 <- 1 - (sum(alignment.matrix[alignment.matrix[, 1] == groups[3], i] == A) /
                          length(alignment.matrix[alignment.matrix[, 1] == groups[3], i]))
            p4 <- 1 - (sum(alignment.matrix[alignment.matrix[, 1] == groups[4], i] == A) /
                          length(alignment.matrix[alignment.matrix[, 1] == groups[4], i]))

            # Durand et al. Equation 2: frequency-based ABBA-BABA
            numerator <- ((1 - p1) * p2 * p3 * (1 - p4)) -
              (p1 * (1 - p2) * p3 * (1 - p4)) + numerator
            denominator <- ((1 - p1) * p2 * p3 * (1 - p4)) +
              (p1 * (1 - p2) * p3 * (1 - p4)) + denominator
          }
        }
      }
    }
    d <- numerator / denominator
    return(list(d, useful, segregating))
  }

  # ---- Calculate D on real data ----
  results <- dpop.calc(alignment.matrix)
  d <- results[[1]]
  if (is.nan(d)) d <- 0

  # ---- Bootstrap ----
  if (sig.test == "B") {
    sim.d <- numeric(replicate)  # pre-allocate for speed
    foo <- ncol(alignment.matrix)
    cat("\nperforming bootstrap")
    for (k in 1:replicate) {
      if (k %% 100 == 0) cat(".")
      # Resample columns (keeping column 1 = names)
      sim.matrix <- alignment.matrix[, c(1, sample.int(foo - 1, replace = TRUE) + 1L)]
      sim.d[k] <- dpop.calc(sim.matrix)[[1]]
    }
    sim.d[is.nan(sim.d)] <- 0
    z <- abs(d / sd(sim.d))
    new.pval <- 2 * (1 - pnorm(z))

    cat("\nSites in alignment =", ncol(alignment.matrix))
    cat("\n\nD raw statistic / Z-score = ", d, " / ", z)
    cat("\n\nResults from ", replicate, "bootstraps")
    cat("\nSD D statistic =", sd(sim.d))
    cat("\nP-value (that D=0) = ", new.pval, "\n\n")
  }

  # ---- Jackknife ----
  if (sig.test == "J") {
    max.rep <- ncol(alignment.matrix) - block.size
    if (block.size >= (ncol(alignment.matrix) / 2)) {
      stop(call. = FALSE,
           paste("\nWith a block size of", block.size,
                 "and an alignment of", ncol(alignment.matrix),
                 "sites \nsome sites would never be included in the \nanalysis",
                 "\n\nThe maximum block size is 1/2 the alignment length"))
    }
    if (max.rep < replicate) {
      stop(call. = FALSE,
           paste("\nWith a block size of", block.size,
                 "and an alignment of", ncol(alignment.matrix),
                 "sites", replicate, "replicates\nare not possible"))
    }
    drop.pos <- seq.int(from = 2, to = (max.rep - 1), length.out = replicate)
    replicate2 <- replicate

    sim.d <- numeric(replicate2)  # pre-allocate for speed
    foo <- ncol(alignment.matrix)
    cat("\nperforming jackknife")
    for (k in 1:replicate2) {
      if (k / 2 == round(k / 2)) cat(".")
      drop_range <- drop.pos[k]:(drop.pos[k] + block.size)
      sim.matrix <- alignment.matrix[, -drop_range]
      sim.d[k] <- dpop.calc(sim.matrix)[[1]]
    }
    sim.d[is.nan(sim.d)] <- 0
    z <- abs(d / sd(sim.d))
    new.pval <- 2 * (1 - pnorm(z))

    cat("\nSites in alignment =", ncol(alignment.matrix))
    cat("\nD raw statistic", d)
    cat("\nZ-score = ", z)
    cat("\n\nResults from", replicate2, "jackknifes with block size of", block.size)
    cat("\nSD D statistic =", sd(sim.d))
    cat("\nP-value (that D=0) = ", new.pval, "\n\n")
  }

  # ---- No significance test ----
  if (sig.test == "N") {
    print(paste("Sites in alignment =", ncol(alignment.matrix) - 1))
    print(paste("Number of sites with ABBA or BABA patterns =", results[[2]]))
    print(paste("Number of ABBA or BABA sites that are still segregating in at least one population =", results[[3]]))
    print(paste("D statistic =", d))
  }

  # Return a list with all relevant results
  user.result <- list()
  user.result$d.stat <- d
  if (sig.test != "N") user.result$pval <- new.pval
  user.result$align.length <- ncol(alignment.matrix) - 1
  user.result$useful.sites <- results[[2]]
  user.result$seg.sites <- results[[3]]
  return(user.result)
}
