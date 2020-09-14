


make_competing_risk_npsem <- function(baseline_covariates,baseline_treatments, censoring_node, risk_nodes, times ) {

  # competing_risks_indicator <- function(data, time, args) {
  #   all_risks <- setdiff(colnames(data), "t")
  #   parents <- args$parents
  #   past_jump <- all(unlist(data[data$t<time, last(.SD), .SDcols = ..all_risks]) == 0)
  #   if(length(parents) > 0){
  #     cur_jump <- all(unlist(data[data$t<=time, last(.SD), .SDcols = ..parents]) == 0)
  #   } else {
  #     cur_jump <- T
  #   }
  #   return(as.numeric(cur_jump & past_jump))
  # }

  competing_risks_indicator <- function(data, time, args, cols) {
    all_risks <- setdiff(cols, "t")
    parents <- args$parents


    past_jump_id <- data[t<time, last(.SD), .SDcols = all_risks, by = id]

    past_jump_id <- past_jump_id$id[rowSums(past_jump_id[, ..all_risks]) == 0]
    if(length(parents) > 0){
      cur_jump_id <- data[id %in% past_jump_id, last(.SD), .SDcols = parents, by = id]
      cur_jump_id <- cur_jump_id$id[rowSums(cur_jump_id[, ..parents]) == 0]
    } else {
      cur_jump_id <- past_jump_id
    }
    set(data, , "keep", as.numeric(data$id %in% cur_jump_id) )
    return(data[, c("id", "keep")])
  }


  npsem <- list()

  npsem[["W"]] <- define_node("W", baseline_covariates, time = 0)
  for(node in baseline_treatments){
    npsem[[node]] <- define_node(node, node, "W", time = 0)
  }
  risk_set_map <- Summary_measure$new(c(censoring_node, risk_nodes, "t"), competing_risks_indicator, args_to_pass = list(parents = risk_nodes), group_by_id = F)
  npsem[[censoring_node]] <-  define_node(censoring_node, censoring_node, c("W", baseline_treatments), time = times, risk_set_map = risk_set_map, missing_row_implies_not_at_risk = F)
  for(i in seq_along(risk_nodes)){
    node <- risk_nodes[[i]]
    nodeName <- node
    if(i != 1){
      nodeName <- paste0(node, "_quasi")
    }
    risk_set_map <- Summary_measure$new(c(censoring_node, risk_nodes, "t"), competing_risks_indicator, args_to_pass = list(parents = risk_nodes[risk_nodes<node]), group_by_id = F)
    npsem[[nodeName]] <- define_node(nodeName, node, c(baseline_covariates, baseline_treatments), time = times, risk_set_map = risk_set_map, missing_row_implies_not_at_risk = F)
  }

  for(i in seq_along(risk_nodes)){
    if(i==1) next
    node <- risk_nodes[[i]]
    nodeName <- node

    risk_set_map <- Summary_measure$new(c(censoring_node, risk_nodes, "t"), competing_risks_indicator, args_to_pass = list(parents = c()), group_by_id = F)
    npsem[[nodeName]] <- define_node(nodeName, node, c(baseline_covariates, baseline_treatments), time = times, risk_set_map = risk_set_map, missing_row_implies_not_at_risk = F)
  }
  return(npsem)
}

make_competing_risk_likelihood <- function(baseline_node, trtment_nodes, censoring_node, competing_risk_nodes, trt_learner = make_learner(Lrnr_glm), competing_risks_learner = make_learner(Lrnr_glm), censoring_learner = make_learner(Lrnr_glm)){
  factor_list <- list()
  factor_list[[baseline_node]] <- LF_emp$new(baseline_node)
  for(node in trtment_nodes){
    factor_list[[node]] <- LF_fit$new(node, trt_learner)
  }
  for(node in competing_risk_nodes){
    factor_list[[node]] <- LF_fit$new(node, competing_risks_learner, is_time_variant = T, type = "mean")
  }
  factor_list[[censoring_node]]  <- LF_fit$new(censoring_node, censoring_learner, is_time_variant = T, type = "mean")
  return(Likelihood$new(factor_list))
}
