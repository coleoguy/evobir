#' Build a Supermatrix from Multiple Sequence Alignments
#'
#' Reads multiple sequence alignment files and either concatenates them
#' into a single supermatrix or pads each one to include all taxa.
#'
#' \strong{What is a supermatrix?}  Think of it like merging
#' spreadsheets.  You have separate DNA sequence files for different
#' genes (e.g., gene1.fasta with 80 species, gene2.fasta with 95
#' species).  This function pastes all the gene columns together
#' side-by-side into one big alignment.  Species that are missing
#' from a particular gene get filled in with the missing data
#' character (default \code{"-"}).
#'
#' A partition file is also saved so you know which columns in the
#' final alignment came from which gene -- this is important for
#' partitioned phylogenetic analyses where each gene gets its own
#' substitution model.
#'
#' @param missing Character.  Symbol to use for missing data
#'   (default \code{"-"}).
#' @param prefix Character.  Prefix for output file names
#'   (default \code{"concatenated"}).
#' @param save Logical.  If \code{TRUE} (default), saves the supermatrix
#'   FASTA and partition CSV to disk.
#' @param input Character.  A pattern to match alignment file names in
#'   the current directory.  If \code{""} (default), all files are used.
#'   For example, \code{input = ".fasta"} would only read files with
#'   ".fasta" in their name.
#' @param format.in Character.  Input alignment format (default \code{"f"}
#'   for FASTA).  Passed to \code{\link[ape]{read.dna}}.
#' @param format.out Character.  Output alignment format (default \code{"f"}
#'   for FASTA).  Passed to \code{\link[ape]{write.dna}}.
#' @param concatenate Logical.  If \code{TRUE} (default), concatenates all
#'   genes into a single supermatrix.  If \code{FALSE}, produces separate
#'   padded alignments (each containing all taxa).
#'
#' @return When \code{concatenate = TRUE}: a list with two elements:
#'   \describe{
#'     \item{partitions}{A data frame with columns: gene file name,
#'       start position, stop position.}
#'     \item{supermatrix}{A character matrix of the concatenated
#'       alignment.}
#'   }
#'   When \code{concatenate = FALSE}: a list of padded alignment matrices,
#'   one per input file.
#'
#' @examples
#' \dontrun{
#' # Concatenate all .fasta files in the current directory
#' setwd("path/to/alignments")
#' result <- SuperMatrix(input = ".fasta")
#' result$partitions  # see which positions belong to which gene
#' }
#'
#' @importFrom ape read.dna write.dna
#' @importFrom utils write.csv
#' @export
SuperMatrix <- function(missing = "-",
                        prefix = "concatenated",
                        save = TRUE,
                        input = "",
                        format.in = "f",
                        format.out = "f",
                        concatenate = TRUE) {

  # ---- Find and read alignment files ----
  file.names <- list.files()
  if (input != "") {
    # fixed = TRUE: treat input as a literal substring, not a regex
    file.names <- grep(input, file.names, value = TRUE, fixed = TRUE)
  }
  # Never read this function's own output files (e.g., from a previous
  # run in the same directory): the supermatrix, the partition CSV, and
  # padded per-gene outputs all begin with the output prefix
  file.names <- file.names[!startsWith(file.names, prefix)]
  if (length(file.names) == 0) {
    stop("no alignment files found matching input = \"", input, "\"")
  }

  # Read each alignment file into a list
  DNA <- list()
  for (i in 1:length(file.names)) {
    print(paste("Reading alignment", i))
    DNA[[i]] <- read.dna(file = file.names[i],
                         format = format.in,
                         as.character = TRUE)
  }

  # ---- Collect all unique taxon names across all alignments ----
  # Vectorized: unlist all rownames, then take unique (preserves first-seen order)
  taxa <- unique(unlist(lapply(DNA, rownames)))

  # ---- Concatenation mode ----
  if (concatenate == TRUE) {
    # Calculate total alignment length (sum of all gene lengths)
    total.length <- sum(vapply(DNA, ncol, integer(1)))

    # Create an empty matrix filled with the missing character
    seqmatrix <- matrix(as.character(missing), length(taxa), total.length)
    rownames(seqmatrix) <- taxa

    # Track where each gene starts and stops (for partitioned analyses)
    partitions <- as.data.frame(matrix(, length(DNA), 3))
    colnames(partitions) <- c("part", "start", "stop")

    # Fill in the supermatrix gene by gene
    c.col <- 0
    print("Creating supermatrix")
    for (i in 1:length(DNA)) {
      print(paste("Processing alignment", i))
      gene <- DNA[[i]]
      # Vectorized: match all taxon names at once instead of one-by-one
      c.rows <- match(rownames(gene), rownames(seqmatrix))
      col_range <- (c.col + 1):(c.col + ncol(gene))
      seqmatrix[c.rows, col_range] <- gene
      partitions[i, 1:3] <- c(file.names[i], c.col + 1, c.col + ncol(gene))
      c.col <- c.col + ncol(gene)
    }

    results <- list()
    results[[1]] <- partitions
    results[[2]] <- seqmatrix

    # Save files if requested
    if (save == TRUE) {
      print("saving files")
      write.dna(seqmatrix,
                file = paste(prefix, ".fasta", sep = ""),
                format = format.out)
      write.csv(partitions,
                row.names = FALSE,
                file = paste(prefix, ".partitions.csv", sep = ""))
    }
  }

  # ---- Non-concatenation mode (pad each alignment separately) ----
  if (concatenate == FALSE) {
    results <- list()
    for (i in 1:length(DNA)) {
      # Create a taxon-complete matrix for this gene
      seqmatrix <- matrix(as.character(missing), length(taxa), ncol(DNA[[i]]))
      rownames(seqmatrix) <- taxa
      # Vectorized: match all taxon names at once
      hits <- match(rownames(DNA[[i]]), rownames(seqmatrix))
      seqmatrix[hits, ] <- DNA[[i]]
      results[[i]] <- seqmatrix
      if (save == TRUE) {
        write.dna(seqmatrix,
                  file = paste(prefix, file.names[i], ".fasta", sep = ""),
                  format = format.out)
      }
    }
  }

  return(results)
}
