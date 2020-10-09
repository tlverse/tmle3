#' Counterfactual Likelihood
#'
#' Represents a counterfactual likelihood where one or more likelihood factors has been replaced with an intervention as specified by \code{intervention_list}.
#' Inherits from \code{\link{Likelihood}}. Other factors (including their updates) are taken from an underlying \code{observed_likelihood} estimated from observed data.
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
#' @template CF_Likelihood_extra
#'
#' @export
CF_Likelihood <- R6Class(
  classname = "CF_Likelihood",
  portable = TRUE,
  class = TRUE,
  inherit = Likelihood,
  public = list(
    initialize = function(observed_likelihood, intervention_list, ...) {
      private$.observed_likelihood <- observed_likelihood
      private$.training_task <- observed_likelihood$training_task
      if (inherits(intervention_list, "LF_base")) {
        intervention_list <- list(intervention_list)
      }

      intervention_nodes <- sapply(intervention_list, `[[`, "name")
      names(intervention_list) <- intervention_nodes
      private$.intervention_list <- intervention_list
      private$.intervention_nodes <- intervention_nodes
      private$.cf_tasks <- self$enumerate_cf_tasks(observed_likelihood$training_task)
      params <- args_to_list()
      super$initialize(params)
    },
    enumerate_cf_tasks = function(tmle_task) {
      intervention_list <- self$intervention_list

      # hack for no intervention
      if (length(intervention_list) == 0) {
        return(list(tmle_task))
      }

      # todo: extend for stochastic interventions
      # currently assumes each intervention node returns one cf_value vector
      # get factors for nodes
      # todo: need to do this sequentially and build the task up based on time ordering
      # because dynamic rules depend on past rule values
      all_values <- lapply(intervention_list, function(likelihood_factor) {
        likelihood_factor$cf_values(tmle_task)
      })

      cf_data <- as.data.table(all_values)
      cf_task <- tmle_task$generate_counterfactual_task(UUIDgenerate(), cf_data)
      cf_tasks <- list(cf_task)
      return(cf_tasks)
    },

    get_likelihood = function(tmle_task, node, fold_number = "full") {
      # todo: this will not handle the case where the cf_likelihood is based on
      # an updated likelihood factor (e.g. old tx shift)
      if (node %in% self$intervention_nodes) {
        likelihood_values <- super$get_likelihood(tmle_task, node, fold_number)
      } else {
        # dispatch to observed likelihood if not an intervention node
        # that way, we get updates to those nodes
        likelihood_values <- self$observed_likelihood$get_likelihood(tmle_task, node, fold_number)
      }

      return(likelihood_values)
    }
  ),
  active = list(
    name = function() {
      node_names <- names(self$intervention_list)
      node_values <- sapply(self$intervention_list, `[[`, "values")
      intervention_name <- paste(sprintf("%s=%s", node_names, as.character(node_values)), collapse = ", ")
      return(intervention_name)
    },
    observed_likelihood = function() {
      return(private$.observed_likelihood)
    },
    intervention_list = function() {
      return(private$.intervention_list)
    },
    intervention_nodes = function() {
      return(private$.intervention_nodes)
    },
    factor_list = function() {
      fl <- self$observed_likelihood$factor_list
      fl[self$intervention_nodes] <- self$intervention_list[self$intervention_nodes]
      return(fl)
    },
    updater = function() {
      self$observed_likelihood$updater
    },
    cf_tasks = function() {
      return(private$.cf_tasks)
    }
  ),
  private = list(
    .observed_likelihood = NULL,
    .intervention_list = NULL,
    .intervention_nodes = NULL,
    .cf_tasks = NULL
  )
)

#' @param ... Passes all arguments to the constructor. See documentation for the
#'  Constructor below.
#'
#' @rdname CF_Likelihood
#'
#' @export
#
make_CF_Likelihood <- CF_Likelihood$new
