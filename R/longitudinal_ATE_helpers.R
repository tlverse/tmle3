
#' @export
ipw_late <- function(task, lik, ipw_args, fold_number){

  cf_likelihood_control = ipw_args$cf_likelihood_control
  cf_likelihood_treatment = ipw_args$cf_likelihood_treatment
  Y <- task$get_tmle_node("Y", format = T)[[1]]
  A_nodes <- grep("A", names(task$npsem), value = T)


  nodes <- c(A_nodes)

  #TODO should do divisions component-wise for large t?
  g <- as.vector(apply(lik$get_likelihoods(task, A_nodes, fold_number = fold_number), 1, prod))
  cf_g_trt <- as.vector(apply(cf_likelihood_treatment$get_likelihoods(task, A_nodes, fold_number = fold_number), 1, prod))
  cf_g_control <- as.vector(apply(cf_likelihood_control$get_likelihoods(task, A_nodes, fold_number = fold_number), 1, prod))

  return(Y*cf_g_trt/g  - Y*cf_g_control/g)

}
#' @export
gradient_generator_late <- function(tmle_task, lik,  node, include_outcome = T, ipw_args = NULL, fold_number){

  task <- tmle_task$get_regression_task(node)
  IC <- ipw_late(tmle_task, lik,  ipw_args, fold_number)
  new_data <- data.table(IC = IC )
  set(new_data, , node, task$Y)
  cols <- task$add_columns(new_data)
  task <- task$clone()
  nodes <- task$nodes
  nodes$outcome <- "IC"
  nodes$covariates <- c(nodes$covariates, node)

  task$initialize(
    task$internal_data,
    nodes = nodes,
    folds = task$folds,
    column_names = cols,
    row_index = task$row_index,
    outcome_type = "continuous"
  )
  return(task)
}

#' @export
generate_npsem_late <- function(baseline_covariates, time_dependent_covariates, time_dependent_treatments,outcome, times ){
  times <- sort(unique(times))
  base_node <- define_node("W", baseline_covariates, time = min(times))
  npsem <- list("W" = base_node)
  parents <- "W"
  for(t in times){
    tcov_name <- c()
    trt_name <-c()
    for(i in seq_along(time_dependent_covariates)){
      name <- paste0("L", t, letters[[i]])
      tcov_name <- c(tcov_name, name)
      cov <- time_dependent_covariates[[i]]
      npsem[[name]] <- define_node(name, cov, parents, time = t)
    }
    for(i in seq_along(time_dependent_treatments)){
      name <- paste0("A", t, letters[[i]])
      trt_name <- c(trt_name, name)
      trt <- time_dependent_treatments[[i]]
      npsem[[name]] <- define_node(name, trt, parents, time = t)
    }
    parents <- c(parents, tcov_name, trt_name)
  }
  npsem[[outcome]] <- define_node("Y", outcome, parents, time = max(times) )
  return(npsem)
}

#' @export
generate_likelihood_late <- function(npsem, trt_learner = make_learner(Lrnr_glm), cov_learner = make_learner(Lrnr_glm), outcome_learner = make_learner(Lrnr_glm)){
  A_nodes <- grep("A", names(npsem), value = T)
  L_nodes <- grep("L", names(npsem), value = T)
  factor_list = list()
  for(node in L_nodes){
    factor_list[[node]] <- LF_fit$new(node, cov_learner)
  }

  for(node in A_nodes){
    factor_list[[node]] <- LF_fit$new(node, trt_learner)
  }
  factor_list[["W"]] <- LF_emp$new("W")
  factor_list[["Y"]] <- LF_fit$new("Y", outcome_learner)
  return(Likelihood$new(factor_list))

}
