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
    initialize = function(initial_likelihood, updater = NULL, ...) {
      params <- args_to_list()

      private$.initial_likelihood <- initial_likelihood

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

      # first, calculate all updates
      task_updates <- lapply(tasks_at_step, self$updater$apply_update, self, fold_number, new_epsilon, update_node)

      # then, store all updates
      for (task_index in seq_along(tasks_at_step)) {
        task <- tasks_at_step[[task_index]]
        updated_values <- task_updates[[task_index]]

        likelihood_factor <- self$factor_list[[update_node]]
        self$cache$set_values(likelihood_factor, task, step_number + 1, fold_number, updated_values)
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
    get_likelihood = function(tmle_task, node, fold_number = "full") {
      if (node %in% self$updater$update_nodes) {
        # self$updater$get_updated_likelihood(self, tmle_task, node)
        likelihood_factor <- self$factor_list[[node]]
        # first check for cached values for this task
        value_step <- self$cache$get_update_step(likelihood_factor, tmle_task, fold_number)

        if (!is.null(value_step)) {
          # if some are available, grab them
          likelihood_values <- self$cache$get_values(likelihood_factor, tmle_task, fold_number)
        } else {
          # if not, generate new ones
          likelihood_values <- self$initial_likelihood$get_likelihood(tmle_task, node, fold_number)
          value_step <- 0
          self$cache$set_values(likelihood_factor, tmle_task, value_step, fold_number, likelihood_values)
        }

        if (value_step < self$updater$step_number) {
          stop(
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
        likelihood_values <- self$initial_likelihood$get_likelihood(tmle_task, node, fold_number)
      }

      return(likelihood_values)
    },
    add_factors = function(factor_list) {
      self$initial_likelihood$add_factors(factor_list)
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
    .updater = NULL
  )
)
