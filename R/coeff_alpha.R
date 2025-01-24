#' @title Calculate Cronbach's Alpha
#'
#' @description
#'   This function calculates the Cronbach's Alpha for a set of measurments
#'
#' @param RT numeric. A vector of reaction times.
#' @param ACC numeric. A vector of accuracies to determine the proportion correct.
#' @param s numeric. The diffusion constant. Typically set to 1 or 0.1. The default value is 1.
#' @param robust logical. Should robust statistics, that is median and Interquartile range, be used to calculate the ezDM parameters. Default is true.
#' @param use_RTs character. Should only correct RTs ("correct") or all RTs (any other character) be used to calculate mean and sd of RTs.
#'
#' @export
coeff_alpha <- function(measurements) {
  # Get number of items (N) and individuals
  n_items <- ncol(measurements)
  n_persons <- nrow(measurements)
  # Get individual total scores
  x <- rowSums(measurements)
  # Get observed-score variance of whole test (X)
  var_x <- var(x) * (n_persons - 1) / n_persons
  # Get observed-score variance of each item (Y_j)
  var_y <- numeric(n_items)
  for (j in 1:n_items) {
    var_y[j] <- var(measurements[, j]) * (n_persons - 1) / n_persons
  }
  # Apply the alpha formula
  alpha <- (n_items / (n_items - 1)) * (1 - sum(var_y) / var_x)

  alpha
}
