

#' Inherits from \code{\link{LF_base}}; see that page for documentation on likelihood factors in general.
#'
#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @family Likelihood objects
#' @keywords data
#'
#' @return \code{LF_base} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @export
LF_binomial_collapsed <- R6Class(
  classname = "LF_binomial_collapsed",
  portable = TRUE,
  class = TRUE,
  inherit = LF_base,
  public = list(
    initialize = function(name, factor_list, parents, main,  ..., type = "density") {
      type = factor_list[[main]]$type
      private$.is_level_variant <- factor_list[[main]]$is_level_variant
      private$.is_time_variant <- factor_list[[main]]$is_time_variant

      super$initialize(name, ..., type = type)
      private$.parents<- parents
      private$.main<- main
      private$.factor_list<- factor_list

    },
    get_mean = function(tmle_task, fold_number, check_at_risk = T, expand = T, drop_id = F, drop_time = F,...) {
      #TODO ignores degeneracy at the moment
      factors <- self$factor_list
      parent_prob <- setDT(lapply(factors[private$.parents], function(factor) {
        factor$get_density(tmle_task, fold_number, quick_pred = T)
      }))

      parent_prob <- as.vector(apply(parent_prob, 1, function(v) {prod(1-v)}))

      main_prob <- factors[[private$.main]]$get_mean(tmle_task, fold_number, quick_pred = T)

      prob <- parent_prob*main_prob

      learner_task <- tmle_task$get_regression_task(self$name, expand = T, include_bins = F, is_time_variant = self$is_time_variant)
      data <- learner_task$get_data()
      if(check_at_risk & "at_risk" %in% colnames(data) ) {
        # By setting check_at_risk = F then one can obtain the counterfactual predictions
        # conditioned on everyone being at risk.
        names_degen_val <- paste0("degeneracy_value_", learner_task$nodes$outcome)
        # TODO Support multivariate outcome
        assertthat::assert_that(all(names_degen_val %in% colnames(data)), msg = "If at_risk is a column then last_val must be as well.")
        not_at_risk <- which(data$at_risk == 0)
        if(length(not_at_risk)>0){
          degen_val <- data[not_at_risk, names_degen_val, with = F]
          #If not at risk then equal to degenerate value with prob 1
          #TODO  we should probably not compute the likelihood for all the not at risk people?
          #Since we change their value anyway?
          prob[not_at_risk] <- degen_val
        }

      }
      prob <- data.table(t = learner_task$id, id = learner_task$time, prob = prob)
      setnames(prob, c("t", "id", node))
      if(drop_id){
        prob$id <- NULL
      }
      if(drop_time){
        prob$t <- NULL
      }
      if(ncol(prob)==1){
        prob <- unlist(prob, use.names  = F)
      }
      if(!expand & check_at_risk){
        prob <- prob[not_at_risk]
      }

      return(prob)
    },
    get_density = function(tmle_task, fold_number, quick_pred = F, expand = T, drop_id = F, drop_time = F, ...) {

      factors <- self$factor_list
      parent_prob <- setDT(lapply(factors[private$.parents], function(factor) {
        factor$get_density(tmle_task, fold_number, quick_pred = T)
      }))

      parent_prob <- as.vector(apply(parent_prob, 1, function(v) {prod(1-v)}))

      main_prob <- factors[[private$.main]]$get_density(tmle_task, fold_number, quick_pred = T)

      prob <- parent_prob*main_prob
      if(quick_pred){
        return(prob)
      }
      observed <- tmle_task$get_tmle_node(private$.main, format = T)[[1]]
      likelihood <- ifelse(observed == 1, prob, 1 - prob)

      learner_task <- tmle_task$get_regression_task(self$name, expand = T, include_bins = F, is_time_variant = self$is_time_variant)
      data <- learner_task$get_data()
      if(check_at_risk & "at_risk" %in% colnames(data) ) {
        # By setting check_at_risk = F then one can obtain the counterfactual predictions
        # conditioned on everyone being at risk.
        names_degen_val <- paste0("degeneracy_value_", learner_task$nodes$outcome)
        # TODO Support multivariate outcome
        assertthat::assert_that(all(names_degen_val %in% colnames(data)), msg = "If at_risk is a column then last_val must be as well.")
        not_at_risk <- which(data$at_risk == 0)
        if(length(not_at_risk)>0){
          degen_val <- data[not_at_risk, names_degen_val, with = F]
          #If not at risk then equal to degenerate value with prob 1
          #TODO  we should probably not compute the likelihood for all the not at risk people?
          #Since we change their value anyway?
          likelihood[not_at_risk] <- as.numeric(observed[not_at_risk] == degen_val)
        }

      }


    likelihood <- data.table(t = learner_task$id, id = learner_task$time, likelihood = likelihood)
    setnames(prob, c("t", "id", node))
    if(drop_id){
      likelihood$id <- NULL
    }
    if(drop_time){
      likelihood$t <- NULL
    }
    if(ncol(likelihood)==1){
      likelihood <- unlist(likelihood, use.names  = F)
    }
    if(!expand & check_at_risk){
      likelihood <- likelihood[not_at_risk]
    }
      return(likelihood)
    }

  ),
  active = list(
    factor_list = function() {
      return(private$.factor_list)
    },
    is_level_variant = function(){
      return(private$.is_level_variant)
    },
    is_time_variant = function(){
      return(private$.is_time_variant)
    }
  ),
  private = list(
    .name = NULL,
    .factor_list = NULL,
    .main = NULL,
    .parents = NULL,
    .is_level_variant = NULL,
    .is_time_variant = NULL
  )
)
