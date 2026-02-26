#' Windowed Patterson's D-Statistic
#'
#' Calculates Patterson's D-statistic in sliding windows across a
#' sequence alignment.  This lets you see how introgression signal
#' varies along the genome -- some regions may show strong introgression
#' while others don't.
#'
#' Like \code{\link{CalcD}}, the alignment must have exactly 4 sequences
#' in the order P1, P2, P3, Outgroup.
#'
#' @param alignment Character.  Path to the FASTA alignment file
#'   (default \code{"alignment.fasta"}).
#' @param win.size Integer.  Width of each window in base pairs
#'   (default 100).
#' @param step.size Integer.  How far the window slides between
#'   calculations (default 50).  Setting \code{step.size < win.size}
#'   gives overlapping windows for smoother results.
#' @param boot Logical.  If \code{TRUE}, performs bootstrap significance
#'   testing within each window.  Default is \code{FALSE}.
#' @param replicate Integer.  Number of bootstrap replicates per window
#'   (default 1000).  Only used when \code{boot = TRUE}.
#'
#' @return A data frame with one row per window and columns:
#'   \describe{
#'     \item{range}{Genomic coordinates of the window (e.g., "1:100").}
#'     \item{abba}{Number of ABBA sites in the window.}
#'     \item{baba}{Number of BABA sites in the window.}
#'     \item{d}{D-statistic for the window.}
#'     \item{Z}{Z-score (only if \code{boot = TRUE}).}
#'     \item{pval}{P-value (only if \code{boot = TRUE}).}
#'   }
#'
#' @details
#' D is calculated as (ABBA - BABA) / (ABBA + BABA) within each window.
#' Windows with no informative sites return D = 0.  When bootstrapping,
#' sites within each window are resampled with replacement to generate
#' a null distribution.
#'
#' @examples
#' \dontrun{
#' # Calculate D in 500-bp windows sliding every 250 bp
#' results <- WinCalcD("my_alignment.fasta",
#'                     win.size = 500, step.size = 250)
#'
#' # Plot D across the genome
#' plot(1:nrow(results), as.numeric(results$d), type = "l",
#'      xlab = "Window", ylab = "D-statistic")
#' abline(h = 0, lty = 2, col = "red")
#' }
#'
#' @importFrom seqinr read.alignment
#' @export
WinCalcD <- function(alignment = "alignment.fasta", win.size = 100,
                     step.size = 50, boot = FALSE, replicate = 1000) {

  # ---- Read and parse the alignment ----
  alignment <- read.alignment(alignment, format = "fasta")
  alignment.matrix <- matrix(, length(alignment$nam),
                             nchar(alignment$seq[[1]]))
  for (i in 1:length(alignment$nam)) {
    alignment.matrix[i, ] <- unlist(strsplit(alignment$seq[[i]], ""))
  }
  full.align <- alignment.matrix

  # ---- Set up sliding window positions ----
  total <- ncol(full.align)
  spots <- seq(from = 1, to = (total - win.size), by = step.size)

  # Prepare results container
  results.matrix <- as.data.frame(matrix(, 1, 6))
  colnames(results.matrix) <- c("range", "abba", "baba", "d", "Z", "pval")

  # ---- Loop through each window ----
  for (q in 1:length(spots)) {
    # Extract the current window from the full alignment
    alignment.matrix <- full.align[, spots[q]:(spots[q] + win.size - 1)]
    starting <- spots[q]
    ending <- spots[q] + win.size - 1

    # Count ABBA and BABA patterns in this window
    abba <- 0
    baba <- 0
    for (i in 1:ncol(alignment.matrix)) {
      # Only biallelic sites
      if (length(unique(alignment.matrix[, i])) == 2) {
        if (alignment.matrix[1, i] != alignment.matrix[2, i]) {
          if (alignment.matrix[4, i] != alignment.matrix[3, i]) {
            if (alignment.matrix[3, i] == alignment.matrix[1, i]) baba <- baba + 1
            if (alignment.matrix[2, i] == alignment.matrix[3, i]) abba <- abba + 1
          }
        }
      }
    }
    d <- (abba - baba) / (abba + baba)
    if (is.nan(d)) d <- 0  # no informative sites -> D = 0

    # ---- Bootstrap within this window ----
    if (boot == TRUE) {
      sim.d <- vector()
      foo <- ncol(alignment.matrix)
      sim.matrix <- matrix(, 4, foo)
      for (k in 1:replicate) {
        sim.matrix[1:4, 1:foo] <- alignment.matrix[1:4, sample(1:foo, replace = TRUE)]
        t.abba <- 0
        t.baba <- 0
        for (i in 1:ncol(sim.matrix)) {
          if (length(unique(sim.matrix[, i])) == 2) {
            if (sim.matrix[1, i] != sim.matrix[2, i]) {
              if (sim.matrix[4, i] != sim.matrix[3, i]) {
                if (sim.matrix[3, i] == sim.matrix[1, i]) t.baba <- t.baba + 1
                if (sim.matrix[2, i] == sim.matrix[3, i]) t.abba <- t.abba + 1
              }
            }
          }
        }
        sim.d[k] <- (t.abba - t.baba) / (t.abba + t.baba)
      }
      sim.d[is.nan(sim.d)] <- 0
      z <- abs(d / sd(sim.d, na.rm = TRUE))
      z[is.nan(z)] <- 0
      new.pval <- 2 * (1 - pnorm(z))

      cat("\nSites in alignment =", ncol(alignment.matrix))
      cat("\nNumber of sites with ABBA pattern =", abba)
      cat("\nNumber of sites with BABA pattern =", baba)
      cat("\n\nD raw statistic / Z-score = ", d, " / ", z)
      cat("\n\nResults from ", replicate, "bootstraps")
      cat("\nSD D statistic =", sd(sim.d, na.rm = FALSE))
      cat("\nP-value (that D=0) = ", new.pval)
      results.matrix[q, 1:6] <- c(paste(starting, ":", ending, sep = ""),
                                   abba, baba, d, z, new.pval)
    }

    # ---- No bootstrap: just report raw D ----
    if (boot == FALSE) {
      cat("\nSites in alignment =", ncol(alignment.matrix))
      cat("\nNumber of sites with ABBA pattern =", abba)
      cat("\nNumber of sites with BABA pattern =", baba)
      cat("\n\nD raw statistic = ", d)
      results.matrix[q, 1:4] <- c(paste(starting, ":", ending, sep = ""),
                                   abba, baba, d)
    }
  }

  return(results.matrix)
}
