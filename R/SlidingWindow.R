#' Sliding Window Analysis
#'
#' Applies a function across a sliding window over a vector or matrix.
#'
#' \strong{What is a sliding window?}  Imagine looking at a long DNA
#' sequence through a small window that you slide along one step at a
#' time.  At each position, you measure something inside the window
#' (e.g., the average GC content).  This gives you a profile showing
#' how that measurement changes along the sequence.  For example, you
#' could calculate the mean GC content in 100-bp windows sliding every
#' 10 bp across a genome, or compute local diversity statistics across
#' a landscape grid.
#'
#' For vector data the window slides along the length; for matrix data
#' the window slides in both dimensions (rows and columns).
#'
#' @param FUN Function to apply within each window.  Must accept a
#'   vector (for vector data) or matrix (for matrix data) and return
#'   a single value.  Can be a function name as a string (e.g.,
#'   \code{"mean"}) or a function object.
#' @param data A numeric vector or matrix to analyze.
#' @param window Integer.  Size of the sliding window (number of
#'   elements for vectors, or number of rows/columns for matrices).
#' @param step Integer.  How far the window moves between calculations.
#'   A step of 1 means windows overlap maximally; a step equal to
#'   \code{window} means no overlap.
#' @param strict Logical.  If \code{TRUE}, validates that \code{data}
#'   is numeric and that \code{window}/\code{step} are sensible sizes.
#'   Default is \code{FALSE}.
#'
#' @return For vector input: a numeric vector of results, one per window
#'   position.
#'
#'   For matrix input: a matrix of results where each cell is the
#'   function applied to the corresponding window.
#'
#' @examples
#' # Calculate mean in windows of 10, stepping by 5
#' x <- rnorm(100)
#' SlidingWindow("mean", x, window = 10, step = 5)
#'
#' # GC content in sliding windows (simplified example)
#' seq_data <- sample(c(0, 1), 500, replace = TRUE)
#' gc <- SlidingWindow("mean", seq_data, window = 50, step = 10)
#'
#' @export
SlidingWindow <- function(FUN, data, window, step, strict = FALSE) {

  # ---- Optional strict validation ----
  if (strict) {
    if (!is.numeric(data)) stop("Please supply numeric data")

    if (sum(is.vector(data), is.matrix(data)) == 0) {
      stop("You must supply data as a vector or matrix")
    }

    if (is.vector(data)) {
      if (window > length(data)) stop("Window size is too large")
      if ((step + window) > length(data)) {
        stop("Window and step size should be small enough that multiple windows can be examined")
      }
    }

    if (is.matrix(data)) {
      if (window > nrow(data)) stop("Window size is too large for number of rows")
      if (window > ncol(data)) stop("Window size is too large for number of cols")
      if ((step + window) > nrow(data)) {
        stop("Window and step size should be small enough that multiple windows can be examined")
      }
    }
  }

  # ---- Vector sliding window (vectorized with vapply) ----
  if (is.vector(data)) {
    total <- length(data)
    # Calculate the starting positions for each window
    spots <- seq(from = 1, to = (total - window + 1), by = step)
    # Resolve function once outside the loop to avoid repeated lookup
    fun <- match.fun(FUN)
    result <- vapply(spots, function(s) fun(data[s:(s + window - 1)]),
                     numeric(1))
  }

  # ---- Matrix sliding window ----
  if (is.matrix(data)) {
    total.x <- ncol(data)
    spots.x <- seq(from = 1, to = (total.x - window + 1), by = step)
    total.y <- nrow(data)
    spots.y <- seq(from = 1, to = (total.y - window + 1), by = step)
    # Resolve function once outside the loop
    fun <- match.fun(FUN)
    result <- matrix(, length(spots.y), length(spots.x))
    for (i in seq_along(spots.y)) {
      row_range <- spots.y[i]:(spots.y[i] + window - 1)
      for (j in seq_along(spots.x)) {
        result[i, j] <- fun(data[row_range,
                                  spots.x[j]:(spots.x[j] + window - 1)])
      }
    }
  }

  if (!exists("result")) stop("Hmmm unknown error... Sorry")

  return(result)
}
