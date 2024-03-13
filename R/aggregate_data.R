#' @title Aggregate data for analyses
#'
#' @description
#'   This function aggregates data for analyses using brms. In most cases this will be done
#'   for analysis of count data that can be modeled with a `binomial` distribution instead
#'   of a `bernoulli` distribution on trial level data.
#'
#' @param formula A `formula` or `brmsformula` object that is specifies the model to be analyzed.
#' @param data A data frame that should be aggregated to be optimized for estimating the specified model
#' @param family The distribution family of the model to be estimated. This will determine how data is aggregated specifically
#'
#' @export
aggregate_data <- function(formula, data, family) {
  if ("brmsformula" %in% class(formula)) {
    vars <- all.vars(formula$formula)
  } else if ("formula" %in% class(formula)) {
    vars <- all.vars(formula)
  }

  resp_var <- dplyr::sym(vars[1])
  vars <- vars[-1]

  data_vars <- dplyr::syms(vars[vars %in% names(data)])

  if (family$family == "binomial") {
    agg_data <- data %>%
      dplyr::group_by(dplyr::pick(!!! data_vars)) %>%
      dplyr::summarise(
        "{{resp_var}}" := sum( !! resp_var ),
        nTrials = dplyr::n(),
        .groups = "drop"
      )
  } else {
    agg_data <- data %>%
      dplyr::group_by(dplyr::pick(!!! data_vars)) %>%
      dplyr::summarise(
        "mean_{{resp_var}}" := mean( !! resp_var ),
        .groups = "drop"
      )
    return(data)
  }

  return(agg_data)
}
