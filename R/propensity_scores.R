#' Get and Plot Propensity Scores
#'
#' @param likelihood a fitted likelihood object
#' @param tmle_task a tmle_task data object
#' @param node a character specifing which node to use
#' @rdname propensity_scores
#' @export
density_formula <- function(tmle_task, node = "A") {
  tmle_node <- tmle_task$npsem[[node]]
  if (tmle_node$variable_type$type == "continuous") {
    operator <- "E"
  } else {
    operator <- "p"
  }
  node_parents <- tmle_node$parents
  dens_form <- sprintf("%s(%s|%s)", operator, node, paste(node_parents, collapse = ", "))

  return(dens_form)
}

#' @rdname propensity_scores
#' @export
get_propensity_scores <- function(likelihood, tmle_task, node = "A") {
  stop("this needs to be re-factored to consider different intervention types. Currently disabled")
  tmle_node <- tmle_task$npsem[[node]]
  # kludge for Rcmd::check with data.table:
  # see https://github.com/Rdatatable/data.table/issues/850
  value <- NULL
  node_values <- tmle_node$variable_type$levels
  node_parents <- tmle_node$parents
  scoremat <- likelihood$get_initial_likelihoods(tmle_task, node, only_observed = FALSE)
  colnames(scoremat) <- node_values
  propensity_scores <- melt.data.table(scoremat, measure.vars = node_values)
  setnames(propensity_scores, c("value", "likelihood"))
  propensity_scores[, value := factor(value)]
  return(propensity_scores)
}

#' @rdname propensity_scores
#' @import ggplot2
#' @export
propensity_score_plot <- function(likelihood, tmle_task, node = "A") {
  propensity_scores <- get_propensity_scores(likelihood, tmle_task, node)
  dens_form <- density_formula(tmle_task, node)
  p <- ggplot(propensity_scores, aes_(x = ~likelihood, fill = ~value)) +
    geom_histogram(binwidth = 0.05) +
    xlab(dens_form) +
    ylab("Density") +
    scale_fill_discrete(name = node) +
    facet_grid(value ~ .) +
    theme_bw() +
    theme(strip.text = element_blank())

  return(p)
}

#' @rdname propensity_scores
#' @export
propensity_score_table <- function(likelihood, tmle_task, node = "A") {
  # kludge for Rcmd::check with data.table:
  # see https://github.com/Rdatatable/data.table/issues/850
  value <- NULL
  propensity_scores <- get_propensity_scores(likelihood, tmle_task, node)
  quants <- propensity_scores[, as.list(quantile(likelihood)), by = list(value)]
  setnames(quants, "value", "A")
  return(quants)
}
