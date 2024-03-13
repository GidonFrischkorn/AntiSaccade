#' @title Plot descriptive statistics
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
plot_descriptives <- function(data, formula, family ,idVar, file = NULL, clean = TRUE, ...) {
  if ("brmsformula" %in% class(bf_binomial)) {
    vars <- all.vars(bf_binomial$formula)
  } else if ("formula" %in% class(bf_binomial)) {
    vars <- all.vars(bf_binomial)
  }

  resp_var <- syms(vars[1:2])
  vars <- vars[-(1:2)]

  vars <- vars[-which(idVar == vars)]

  data_vars <- dplyr::syms(vars[vars %in% names(data)])
  if (length(data_vars) > 1) {
    color_var <- data_vars[[2]]
  }

  plot <- ggplot2::ggplot(data = agg_data,
                  aes(x = !! data_vars[[1]], y = !!resp_var[[1]]/ !!resp_var[[2]],
                      color = !! color_var, fill = !! color_var, shape_var = !! color_var)) +
    stat_summary(position = position_dodge(0.25)) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = .25),
                alpha = 0.5)

  if (clean) {
    plot <- clean_plot(plot)
  }

}

#' @title Adapt the ggplot theme to clean plots up
#'
#' @param plot A `ggplot2` object that should be clean up
#' @param ... additional cleaning settings that should be applied
#'
#' @export
clean_plot <- function(plot, ...) {
  plot <- plot +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line.x = element_line(color = 'black'),
          axis.line.y = element_line(color = 'black'),
          legend.key = element_rect(fill = 'white'),
          text = element_text(size = 18),
          line = element_line(linewidth = 1),
          axis.ticks = element_line(linewidth = 1),
          ...)
}
