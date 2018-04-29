#' Defines an update (submodel+loss function)
#'
#' Current Limitations:
#' only updating one node
#' loss function and submodel are hard-coded (need to accept arguments for these)
#' no support for one-step (recursive TMLE)
#' @importFrom R6 R6Class
#'
#' @export
#
tmle3_Update <- R6Class(
  classname = "tmle3_Update",
  portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(){
      
    },
    update_step = function(likelihood, tmle_task) {
      update_node <- self$update_nodes[[1]]
      likelihood_values <- likelihood$get_likelihoods(tmle_task, update_node)
      submodel_data <- self$generate_submodel_data(tmle_task, likelihood_values, update_node)
      new_epsilon <- self$fit_submodel(submodel_data)
      updated_likelihood <- self$apply_submodel(submodel_data, new_epsilon)
      likelihood$update(update_node, updated_likelihood)
    },
    generate_submodel_data = function(tmle_task, likelihood_values, update_node) {
      #todo: support not getting observed for case where we're applying updates instead of fitting them
      clever_covariates <- lapply(self$tmle_params, function(tmle_param) unlist(tmle_param$clever_covariates(tmle_task)))
      dt <- do.call(cbind, clever_covariates)
      # clever_covariates <- tmle_param$clever_covariates(tmle_task)
      submodel_data <- list(
        observed = tmle_task$get_tmle_node(update_node, bound = TRUE),
        H = dt,
        initial = unlist(likelihood_values)
      )
      return(submodel_data)
    },
    fit_submodel = function(submodel_data) {
      # fit submodel
      # submodel function might be predict here, but generally _is_ a function we're trying to fit
      suppressWarnings({
        submodel_fit <- glm(observed~H - 1, submodel_data, offset = qlogis(submodel_data$initial), family = binomial())
      })
      epsilon <- coef(submodel_fit)
      private$.epsilons <- c(private$.epsilons, list(epsilon))
      return(epsilon)
    },
    submodel = function(epsilon, initial, H) {
      plogis(qlogis(initial) + H %*% epsilon)
    },
    # submodel = function(epsilon, initial,H){
    #   initial+H%*%epsilon
    # },
    loss_function = function(estimate, observed) {
      -1 * ifelse(observed == 1, log(estimate), log(1 - estimate))
    },
    apply_submodel = function(submodel_data, epsilon) {
      updated_likelihood <- self$submodel(epsilon, submodel_data$initial, submodel_data$H)
      return(updated_likelihood)
    },
    apply_updates = function(tmle_task, likelihood, initial_likelihood, epsilons) {
      update_node <- self$update_nodes[[1]]
      updated_likelihood <- initial_likelihood
      for (epsilon in epsilons) {
        submodel_data <- self$generate_submodel_data(tmle_task, updated_likelihood, update_node)
        updated_likelihood <- self$apply_submodel(submodel_data, epsilon)
      }

      return(updated_likelihood)
    }
  ),
  active = list(
    epsilons = function() {
      return(private$.epsilons)
    },
    tmle_params = function(new_params = NULL) {
      if(!is.null(new_params)){
        if (inherits(new_params, "Param_base")) {
        new_params <- list(new_params)
      }
      private$.tmle_params <- new_params
      private$.update_nodes <- unique(unlist(lapply(new_params, `[[`, "update_nodes")))

      }
      return(private$.tmle_params)
    },
    update_nodes = function() {
      return(private$.update_nodes)
    },
    step_number = function(){
      return(length(self$epsilons))
    }
  ),
  private = list(
    .epsilons = list(),
    .tmle_params = NULL,
    .update_nodes = NULL
  )
)
