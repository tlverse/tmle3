#' Helper functions to get and plot propensity scores
#' 
#' @param likelihood a fitted likelihood object
#' @param tmle_task a tmle_task data object
#' @param node a character specifing which node to use
#' @rdname propensity_scores
#' @export
density_formula <- function(tmle_task, node = "A"){
  tmle_node <- tmle_task$tmle_nodes[[node]]
  node_parents <- tmle_node$parents
  dens_form <- sprintf("p(%s|%s)", node, paste(node_parents, collapse=", "))
  
  return(dens_form)
}

#' @rdname propensity_scores
#' @export
get_propensity_scores <- function(likelihood, tmle_task, node="A"){
  tmle_node <- tmle_task$tmle_nodes[[node]]
  node_values <- tmle_node$variable_type$levels
  node_parents <- tmle_node$parents
  lf <- likelihood$get_factor(node)
  scoremat <- lf$get_likelihood(tmle_task, only_observed = FALSE)
  colnames(scoremat) = node_values
  propensity_scores <- melt(scoremat)
  setDT(propensity_scores)
  
  setnames(propensity_scores,c("junk", "value", "likelihood"))
  propensity_scores[ , junk:=NULL]
  propensity_scores[ , value:=factor(value)]
  return(propensity_scores)
}

#' @rdname propensity_scores
#' @export
propensity_score_plot <- function(likelihood, tmle_task, node="A"){
  propensity_scores <- get_propensity_scores(likelihood, tmle_task, node)
  dens_form <- density_formula(tmle_task, node)
  p <- ggplot(propensity_scores,aes(x=likelihood, fill=value)) + geom_histogram(binwidth = 0.05) +
    xlab(dens_form) + ylab("Density") + scale_fill_discrete(name=node) + 
    facet_grid(value~.) + theme_bw() + theme(strip.text = element_blank())
  
  return(p)
}

#' @rdname propensity_scores
#' @export
propensity_score_table <- function(likelihood, tmle_task, node="A"){
  propensity_scores <- get_propensity_scores(likelihood, tmle_task, node)
  quants <- propensity_scores[, as.list(quantile(likelihood)), by=list(value)]
  setnames(quants, "value", "A")
  return(quants)  
}
