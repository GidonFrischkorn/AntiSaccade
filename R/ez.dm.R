#' @title Calculate ezDM parameters
#'
#' @description
#'   This functions calculates the ezDM parameters as proposed by Wagenmakers et al. (2007).
#'   For this it takes a vector of reaction times and accuracies, calculates the mean and sd
#'   of the reaction times and the proportion of correct responses and returns the ezDM
#'   parameters.
#'
#' @param RT numeric. A vector of reaction times  A `formula` or `brmsformula` object that is specifies the model to be analyzed.
#' @param ACC numeric. A vector of accuracies to determine the proportion correct.
#' @param s numeric. The diffusion constant. Typically set to 1 or 0.1. The default value is 1.
#' @param robust logical. Should robust statistics, that is median and Interquartile range, be used to calculate the ezDM parameters. Default is true.
#' @param use_RTs character. Should only be correct RTs ("correct") or all RTs (any other character) be used to calculate mean and sd of RTs.
#'
#' @export
ez.dm <- function(RT, ACC, s = 1, robust = TRUE, use_RTs = "correct"){
  # The default value for the scaling parameter s equals 0.1
  s2 <- s^2

  # compute descriptives
  N = length(ACC)
  Pc = mean(ACC)

  if (use_RTs == "correct") {
    useRTs <- RT[which(ACC == 1)]
  } else {
    useRTs <- RT
  }

  if (robust == TRUE) {
    MRT = median(useRTs)
    VRT = (IQR(useRTs)/(qnorm(p = .75)*2))^2
  } else {
    MRT = mean(useRTs)
    VRT = var(useRTs)
  }


  # edge correction is required.
  if (Pc == 0 | Pc == 0.5) {
    Pc = Pc + (1/(2*N))
  } else if (Pc == 1) {
    Pc = Pc - (1/(2*N))
  }

  # The function "qlogis" calculates the logit.
  L <- qlogis(Pc)

  # This gives drift rate.
  x <- L*(L*Pc^2 - L*Pc + Pc - .5)/VRT
  v <- sign(Pc - .5)*s*x^(1/4)

  # This gives boundary separation.
  a <- s2*L/v

  # This gives nondecision time.
  y <- -v*a/s2
  MDT <- (a/(2*v)) * (1 - exp(y))/(1 + exp(y))
  t0 <- MRT - MDT

  # Return the calculated parameters
  parms <- c(v,a,t0)
  names(parms) <- c("v","a","t0")
  return(parms)
}
