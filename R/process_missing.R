#' @export
process_missing <- function(data, node_list){
  # todo - do IPCW instead of dropping these
  drop_vars <- c(node_list$A, node_list$Y)
  to_drop <- data[, apply(is.na(.SD), 1, any), .SDcols=drop_vars]
  filtered <- data[!to_drop]
  
  
  p_missing <- sapply(filtered[,impute_vars, with=FALSE], function(x)mean(is.na(x)))
  to_impute <- names(p_missing[(0 < p_missing) & (p_missing < max_p_missing)])
  no_missing <- names(p_missing[p_missing==0])
  no_missing <- c(no_missing, drop_vars)
  missing_indicators <- filtered[, lapply(.SD, is.na), .SDcols=to_impute]
  missing_names <- sprintf("delta_%s" ,any_missing)
  setnames(missing_indicators, missing_names)
  impute_median <- function(x){
    value <- median(as.numeric(x[!is.na(x)]))
    x[is.na(x)] <- value
    x
  }
  imputed <- filtered[, lapply(.SD, impute_median), .SDcols=to_impute]
  all_observed <- filtered[, no_missing, with=FALSE]
  
  processed <- cbind(all_observed, imputed, missing_indicators)
  updated_nodes <- lapply(node_list, function(node){
    node_no_missing <- intersect(no_missing, node)
    node_imputed <- c(to_impute,missing_names)[to_impute%in%node]
    return(c(node_no_missing, node_imputed))
  })
  
  return(list(data=processed, node_list=updated_nodes))
}