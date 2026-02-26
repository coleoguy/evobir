#' Reorder Sequences in a FASTA Alignment
#'
#' Reads a FASTA alignment and writes a new file with the sequences
#' reordered by their starting position (i.e., where the first non-gap
#' character appears).
#'
#' \strong{When would you use this?}  When you have an alignment where
#' different sequences cover different regions -- for example, amplicon
#' sequencing data where each primer pair targets a different part of
#' a gene, or partial genome sequences that start at different
#' positions.  Sorting by start position makes the alignment much
#' easier to visualize and inspect in an alignment viewer.
#'
#' Optionally, one or more reference sequences can be pinned to the
#' top of the output file regardless of their start position.
#'
#' @param file Character.  Path to the input FASTA alignment file.
#' @param newfile Character.  Path for the output (reordered) FASTA file.
#' @param ref Optional.  Reference sequence(s) to place first in the
#'   output.  Can be specified as:
#'   \itemize{
#'     \item A character vector of sequence names.
#'     \item A numeric vector of sequence indices (positions in the file).
#'     \item \code{NULL} (default) -- no pinned references.
#'   }
#'
#' @return Writes the reordered alignment to \code{newfile}.  No value is
#'   returned (called for its side effect of writing a file).
#'
#' @examples
#' \dontrun{
#' # Reorder by start position
#' ReOrderAlignment("my_alignment.fasta", "reordered.fasta")
#'
#' # Pin "reference_seq" to the top, then sort the rest
#' ReOrderAlignment("my_alignment.fasta", "reordered.fasta",
#'                  ref = "reference_seq")
#' }
#'
#' @importFrom seqinr read.fasta write.fasta
#' @export
ReOrderAlignment <- function(file, newfile, ref = NULL) {
  # Read the FASTA file into a list of sequences
  dat <- read.fasta(file)

  # Validate reference sequences if provided
  if (!is.null(ref)) {
    # Convert numeric indices to names
    if (is.numeric(ref[1])) {
      ref <- names(dat)[ref]
    }
    if (!sum(ref %in% names(dat)) == length(ref)) {
      stop("One of your reference sequences was not found in your file")
    }
  }

  # ---- Find where each sequence starts (first non-gap position) ----
  startspots <- c()
  for (i in 1:length(dat)) {
    if (dat[[i]][1] != "-") {
      # Sequence starts at position 1 (no leading gaps)
      startspots[i] <- 1
    } else {
      # Find the run of leading gaps and measure its length
      bar <- which(dat[[i]] == "-")
      vec <- bar[-1] - bar[-length(bar)]
      if (length(vec) == sum(vec)) {
        # All gaps are contiguous from the start
        startspots[i] <- length(vec)
      } else {
        # Leading gaps end where the first non-consecutive gap appears
        startspots[i] <- min(which(vec != 1))
      }
    }
  }

  # Sort sequences by their start position
  seqorder <- order(startspots)

  # Move reference sequences to the front
  starts <- which(names(dat) %in% ref)
  seqorder <- c(starts, seqorder[!seqorder %in% starts])

  # Write the reordered alignment
  newlist <- dat[seqorder]
  write.fasta(newlist, names = names(newlist), file.out = newfile)
  print("Done!")
}
