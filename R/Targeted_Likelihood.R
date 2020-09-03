#' Targeted Likelihood
#'
#' Represents a likelihood where one or more likelihood factors has been updated
#' to target a set of parameter(s)
#' @importFrom R6 R6Class
#' @importFrom sl3 Lrnr_base args_to_list
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @family Likelihood objects
#' @keywords data
#'
#' @return \code{Likelihood} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @template Likelihood_extra
#'
#' @export
Targeted_Likelihood <- R6Class(
  classname = "Targeted_Likelihood",
  portable = TRUE,
  class = TRUE,
  inherit = Likelihood,
  public = list(
    initialize = function(initial_likelihood, updater = NULL, submodel_type_by_node = "logistic", check_sync = T, ...) {
      params <- args_to_list()
      params$submodel_type_by_node <- NULL
      params$check_sync <- NULL
      private$.initial_likelihood <- initial_likelihood
      private$.submodel_type_by_node <- submodel_type_by_node
      private$.check_sync <- check_sync
      # handle updater arguments
      if (is.null(updater)) {
        updater <- tmle3_Update$new()
      } else if (inherits(updater, "tmle3_Update")) {
        # do nothing
      } else if (inherits(updater, "list")) {
        # construct updater from list arguments
        updater <- do.call(tmle3_Update$new, updater)
      }
      private$.updater <- updater

      super$initialize(params)
    },
    update = function(new_epsilon, step_number, fold_number = "full", update_node) {

      # todo: rethink which tasks need updates here
      # tasks_at_step <- self$cache$tasks_at_step(step_number)
      tasks_at_step <- self$cache$tasks
      # If task has attr target_nodes then only update these nodes.
      to_update <- sapply(tasks_at_step, function(task) {
        target_nodes <- attr(task, "target_nodes")
        if(is.null(target_nodes)){
          return(T)
        }
        if(update_node %in% c(target_nodes)){
          return(T)
        }
        return(F)
      })
      tasks_at_step <- tasks_at_step[to_update]
      # first, calculate all updates
      task_updates <- lapply(tasks_at_step, self$updater$apply_update, self, fold_number, new_epsilon, update_node)

      # then, store all updates
      for (task_index in seq_along(tasks_at_step)) {
        task <- tasks_at_step[[task_index]]
        updated_values <- task_updates[[task_index]]

        likelihood_factor <- self$factor_list[[update_node]]
        self$cache$set_values(likelihood_factor, task, step_number + 1, fold_number, updated_values, node = update_node)
      }
      # for (task in tasks_at_step) {
      #   all_submodels <- self$updater$generate_submodel_data(self, task, fold_number)
      #   updated_values <- self$updater$apply_submodels(all_submodels, new_epsilons)
      #   for (node in names(updated_values)) {
      #     likelihood_factor <- self$factor_list[[node]]
      #     self$cache$set_values(likelihood_factor, task, step_number + 1, fold_number, updated_values[[node]])
      #   }
      # }
    },
    sync_task = function(tmle_task, fold_number = "full", check = T, max_step = NULL){
      # Takes a task and syncs it with current update status of likelihood
      # Returns task invisibly.

      # Checks if task is up to date
      # Might be useful to set check = F for efficiency ...
      # e.g. if you are montecarlo simulating from a targeted likelihood

      if(check){
        current_step <- self$updater$step_number
        nodes <- self$updater$update_nodes
        target_nodes <- attr(tmle_task, "target_nodes")
        if(!is.null(target_nodes)){
          nodes <- target_nodes
        }
        in_sync <- TRUE
        for(node in nodes){
          likelihood_factor <- self$factor_list[[node]]
          step_number <- self$cache$get_update_step(likelihood_factor, tmle_task, fold_number, node = node)
          if(step_number != current_step) {
            in_sync <- FALSE
          }
        }
        if(in_sync) {
          return(invisible(tmle_task))
        }
      }


      private$.check_sync <- F
      epsilons <- self$updater$epsilons
      step_count <- 0
      if(is.null(max_step)) max_step <- length(epsilons)
      #TODO is this double for loop something to worry about?
      for(eps_step in epsilons){
        if(step_count >= max_step){
          break
        }
        for(node in names(eps_step)){
          target_nodes <- attr(tmle_task, "target_nodes")
          if(!is.null(target_nodes) & !(node %in% c(target_nodes))){
            next
          }
          eps <- eps_step[[node]]
          likelihood_factor <- self$factor_list[[node]]

          step_number <- self$cache$get_update_step(likelihood_factor, tmle_task, fold_number, node = node)
          if(is.null(step_number)){
            step_number <- 0
          }
          if(step_number > step_count){
            next
          }
          updated_likelihoods <-  self$updater$apply_update(tmle_task, self, fold_number, eps, node)
          self$cache$set_values(likelihood_factor, tmle_task, step_number + 1, fold_number, updated_likelihoods, node)

        }
        step_count <- step_count + 1

      }
      return(invisible(tmle_task))
    },
    get_likelihood = function(tmle_task, node, fold_number = "full",  check_sync = NULL, ...) {
      if(is.null(check_sync)){
        check_sync <- private$.check_sync
      }
      if (node %in% self$updater$update_nodes) {
        # self$updater$get_updated_likelihood(self, tmle_task, node)
        likelihood_factor <- self$factor_list[[node]]
        # first check for cached values for this task
        value_step <- self$cache$get_update_step(likelihood_factor, tmle_task, fold_number, node = node)

        if (!is.null(value_step)) {
          # if some are available, grab them
          likelihood_values <- self$cache$get_values(likelihood_factor, tmle_task, fold_number, node = node)
        } else {
          # if not, generate new ones
          likelihood_values <- self$initial_likelihood$get_likelihood(tmle_task, node, fold_number, to_wide = F)
          value_step <- 0
          self$cache$set_values(likelihood_factor, tmle_task, value_step, fold_number, likelihood_values, node = node)
        }

        if (value_step < self$updater$step_number & check_sync) {
          stop(
            "cached likelihood value is out of sync with updates\n",
            "lf_uuid: ", likelihood_factor$uuid, "\n",
            "task_uuid: ", tmle_task$uuid, "\n",
            "node: ", node, " fold_number: ", fold_number, "\n",
            "cached_step: ", value_step, "\n",
            "update_step: ", self$updater$step_number, "\n"
          )
        } else if(value_step < self$updater$step_number){
          warning(
            "cached likelihood value is out of sync with updates\n",
            "lf_uuid: ", likelihood_factor$uuid, "\n",
            "task_uuid: ", tmle_task$uuid, "\n",
            "node: ", node, " fold_number: ", fold_number, "\n",
            "cached_step: ", value_step, "\n",
            "update_step: ", self$updater$step_number, "\n"
          )
        }
        # todo: maybe update here, or error if not already updated
      } else {
        # not a node that updates, so we can just use initial likelihood
        likelihood_values <- self$initial_likelihood$get_likelihood(tmle_task, node, fold_number, ...)
      }

      return(likelihood_values)
    },
    add_factors = function(factor_list) {
      self$initial_likelihood$add_factors(factor_list)
    },
    submodel_type = function(node){
      if(is.character(private$.submodel_type_by_node)){
        return(private$.submodel_type_by_node)
      } else {
        return(private$.submodel_type_by_node[[node]])
      }
    }
  ),
  active = list(
    name = function() {
      node_names <- names(self$intervention_list)
      node_values <- sapply(self$intervention_list, `[[`, "values")
      intervention_name <- paste(sprintf("%s=%s", node_names, as.character(node_values)), collapse = ", ")
      return(intervention_name)
    },
    initial_likelihood = function() {
      return(private$.initial_likelihood)
    },
    updater = function() {
      return(private$.updater)
    },
    factor_list = function() {
      return(self$initial_likelihood$factor_list)
    },
    training_task = function() {
      return(self$initial_likelihood$training_task)
    },
    censoring_nodes = function() {
      return(self$initial_likelihood$censoring_nodes)
    }
  ),
  private = list(
    .initial_likelihood = NULL,
    .updater = NULL,
    .submodel_type_by_node = NULL,
    .check_sync = NULL
  )
)
