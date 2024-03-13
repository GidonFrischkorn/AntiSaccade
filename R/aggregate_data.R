aggregate_data <- function(formula, data, family) {
  if ("brmsformula" %in% class(formula)) {
    vars <- all.vars(formula$formula)
  } else if ("formula" %in% class(formula)) {
    vars <- all.vars(formula)
  }

  resp_var <- sym(vars[1])
  vars <- vars[-1]

  data_vars <- dplyr::syms(vars[vars %in% names(data)])

  agg_data <- data %>%
    dplyr::group_by(dplyr::pick(!!! data_vars)) %>%
    dplyr::summarise(
      "{{resp_var}}" := sum( !! resp_var ),
      nTrials = dplyr::n(),
      .groups = "drop"
    )

  return(agg_data)
}
