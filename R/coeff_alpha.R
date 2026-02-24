#' @title Calculate Cronbach's Alpha
#'
#' @description
#'   This function calculates the Cronbach's Alpha for a set of measurments
#'
#' @param measurements numeric matrix or data frame. Rows represent individuals,
#'   columns represent items/measurements.
#'
#' @importFrom stats var
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
