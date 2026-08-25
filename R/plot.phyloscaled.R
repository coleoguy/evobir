#' Plot a Rate-Scaled Phylogeny
#'
#' S3 plot method for \code{"phyloscaled"} objects produced by
#' \code{\link{scaleTreeRates}}.  Visualizes how evolutionary rates
#' vary across the tree using either branch length scaling or color
#' coding.
#'
#' @param x A tree of class \code{"phyloscaled"} (output of
#'   \code{\link{scaleTreeRates}}).  Contains a \code{$scalar} element
#'   indicating the rate multiplier for each edge.
#' @param method Character.  How to display rate variation:
#'   \itemize{
#'     \item \code{"multiply"} (default) -- multiply each edge length
#'       by its scalar, so fast branches appear longer and slow branches
#'       appear shorter.
#'     \item \code{"color"} -- color edges by their scalar value.
#'       Red/cool colors = slow evolution (scalar < 1), green/warm
#'       colors = fast evolution (scalar > 1).
#'   }
#' @param palette Character.  Color palette for the \code{"color"} method.
#'   Default is \code{"RdYlGn"} (from RColorBrewer).  Automatically
#'   switches to \code{"viridis"} if there are more than 10 categories.
#' @param edge.width Numeric.  Width of plotted edges (default 1).
#' @param cex Numeric.  Character expansion factor for tip labels
#'   (default 1).
#' @param show.tip.label Logical.  Whether to display tip labels
#'   (default \code{TRUE}).
#' @param ... Additional arguments passed to
#'   \code{\link[ape]{plot.phylo}}.
#'
#' @return Produces a plot.  Returns \code{NULL} invisibly.
#'
#' @examples
#' \dontrun{
#' # After running scaleTreeRates:
#' scaled_tree <- scaleTreeRates(tree, tip.states, model = "ARD")
#'
#' # Method 1: scale branch lengths
#' plot(scaled_tree, method = "multiply")
#'
#' # Method 2: color branches by rate
#' plot(scaled_tree, method = "color")
#' }
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom viridis viridis
#' @importFrom ape plot.phylo
#' @exportS3Method plot phyloscaled
plot.phyloscaled <- function(x, method = "multiply", palette = "RdYlGn",
                             edge.width = 1, cex = 1,
                             show.tip.label = TRUE, ...) {
  tree <- x

  # Validate that this is actually a phyloscaled object
  if (!"phyloscaled" %in% class(tree)) {
    print("tree is not of class phyloscaled")
    stop()
  }

  # ---- Method: multiply branch lengths by scalar ----
  if (method == "multiply") {
    tree$edge.length <- tree$edge.length * tree$scalar
    plot.phylo(tree, edge.width = edge.width, cex = cex,
               show.tip.label = show.tip.label)

  # ---- Method: color branches by scalar value ----
  } else if (method == "color") {
    # Get the unique scalar values and how many there are
    states <- sort(unique(tree$scalar))
    states.numeric <- as.numeric(as.factor(tree$scalar))
    n <- length(states)

    # If too many categories for brewer.pal, switch to viridis
    if (n >= 11) {
      print("number of states exceeds maximum for brewer.pal, using viridis")
      palette <- "viridis"
    }

    # Count how many scalars are below 1 (slow) and above 1 (fast)
    slow <- sum(states < 1)
    fast <- sum(states > 1)

    if (slow >= 5 || fast >= 5) {
      print("number of states exceeds maximum for brewer.pal, using viridis")
      palette <- "viridis"
    }

    # ---- Generate color palette ----
    if (palette == "RdYlGn") {
      # Use diverging palette centered on 1 (no change)
      if (1 %in% states) {
        cols <- c(brewer.pal(n = 2 * slow + 1, name = palette)[1:slow],
                  brewer.pal(n = 3, name = palette)[2],
                  brewer.pal(n = 2 * fast + 1, name = palette)[(fast + 1):
                                                                  (2 * fast + 1)])
      } else {
        cols <- c(brewer.pal(n = 2 * slow + 1, name = palette)[1:slow],
                  brewer.pal(n = 2 * fast + 1, name = palette)[(fast + 2):
                                                                  (2 * fast + 1)])
      }
    } else {
      # Viridis: continuous palette split around the midpoint
      if (1 %in% states) {
        cols <- c(viridis(n = slow, begin = 0, end = 0.49),
                  viridis(n = 1, begin = 0.5, end = 0.5),
                  viridis(n = fast, begin = 0.51, end = 1))
      } else {
        cols <- c(viridis(n = slow, begin = 0, end = 0.49),
                  viridis(n = fast, begin = 0.51, end = 1))
      }
    }

    # Map each edge to its color
    edge.cols <- cols[states.numeric]

    plot.phylo(tree, edge.color = edge.cols, edge.width = edge.width,
               cex = cex, show.tip.label = show.tip.label)

  } else {
    stop("unrecognized method")
  }
}
