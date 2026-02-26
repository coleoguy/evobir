#' Residual Selection for Artificial Breeding Experiments
#'
#' Identifies individuals with the most extreme residual values from a
#' linear regression between two traits.  This is a common method in
#' artificial selection experiments.
#'
#' \strong{What are residuals?}  Imagine you plot body size (x-axis)
#' against wing length (y-axis) and draw a best-fit line through the
#' data.  A "residual" is how far above or below that line each
#' individual falls.  Individuals far ABOVE the line have unusually
#' long wings for their body size (positive residuals).  Individuals
#' far BELOW the line have unusually short wings for their body size
#' (negative residuals).
#'
#' The function selects the most extreme individuals from each end:
#' the "high line" (furthest above the line) and the "low line"
#' (furthest below the line) for use as breeding stock.  It also
#' produces a scatter plot showing the regression line (orange),
#' high-line individuals (blue), and low-line individuals (red).
#'
#' @param data A data frame containing at least an identifier column and
#'   the two trait columns.
#' @param traits A numeric vector of length 2 giving the column numbers
#'   of the two traits.  The first is the predictor (x-axis), the second
#'   is the response (y-axis).  For example, \code{c(2, 3)} means
#'   column 2 predicts column 3.
#' @param percent Numeric.  The percentage of individuals to select from
#'   each tail of the residual distribution (default 10, meaning the
#'   top 10\% and bottom 10\%).  Higher values select more individuals
#'   but with less extreme phenotypes.
#' @param identifier Integer.  Column number containing individual IDs
#'   (default 1).
#' @param model Character.  Type of regression model (default \code{"linear"}).
#'   Currently only linear models are supported.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{high line}{Identifiers of individuals in the upper tail
#'       (positive residuals -- their y-trait is bigger than expected
#'       given their x-trait value).}
#'     \item{low line}{Identifiers of individuals in the lower tail
#'       (negative residuals -- their y-trait is smaller than expected
#'       given their x-trait value).}
#'   }
#'
#' @examples
#' \dontrun{
#' # Simulate some data: 100 beetles with body size and wing length
#' df <- data.frame(
#'   id = paste0("beetle_", 1:100),
#'   body_size = rnorm(100, 50, 10),
#'   wing_length = rnorm(100, 30, 5)
#' )
#' # Select top and bottom 15% by residual wing length
#' # (controlling for body size)
#' selected <- ResSel(df, traits = c(2, 3), percent = 15)
#' selected$`high line`  # beetles with unusually LONG wings for their size
#' selected$`low line`   # beetles with unusually SHORT wings for their size
#' }
#'
#' @importFrom stats lm
#' @importFrom graphics plot points abline
#' @export
ResSel <- function(data, traits, percent = 10, identifier = 1,
                   model = "linear") {
  # Fit a linear model: trait2 ~ trait1
  lm.norm <- lm(data[, traits[2]] ~ data[, traits[1]])

  # Rank the residuals from smallest (most negative) to largest (most positive)
  foo <- rank(lm.norm$residuals)

  # Identify individuals in the upper and lower tails
  # "high line" = individuals whose trait2 is higher than predicted
  # "low line" = individuals whose trait2 is lower than predicted
  high <- data[foo > (nrow(data) - (nrow(data) * 0.01 * percent)), identifier]
  low  <- data[foo < (nrow(data) * 0.01 * percent), identifier]

  # ---- Visualization ----
  # Gray dots = all individuals, blue = high line, red = low line
  plot(data[, traits[1]:traits[2]], pch = 19, cex = 0.5,
       main = paste("Top and Bottom ", percent, "%", sep = ""), col = "gray")
  points(data[foo > (nrow(data) - (nrow(data) * 0.01 * percent)),
              traits[1]:traits[2]],
         col = "blue", pch = 19, cex = 1.1)
  points(data[foo < (nrow(data) * 0.01 * percent),
              traits[1]:traits[2]],
         col = "red", pch = 19, cex = 1.1)
  abline(lm.norm, col = 'orange', lwd = 5)

  # Return the selected individuals
  results <- list()
  results[[1]] <- high
  results[[2]] <- low
  names(results) <- c('high line', 'low line')
  return(results)
}
