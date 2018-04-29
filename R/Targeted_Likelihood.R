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
#' @template CF_Likelihood_extra
#'
#' @export
Targeted_Likelihood <- R6Class(
  classname = "Targeted_Likelihood",
  portable = TRUE,
  class = TRUE,
  inherit = Likelihood,
  public = list(
    initialize = function(initial_likelihood, updater, ...) {
      params <- args_to_list()
      
      private$.initial_likelihood <- initial_likelihood
      private$.updater <- updater
      super$initialize(params)
    },
    update = function(update_node, updated_likelihood){
      likelihood_factor <- self$factor_list[[update_node]]
      self$cache$set_values(likelihood_factor, self$initial_likelihood$training_task, self$updater$step_number, updated_likelihood)
    },
    get_likelihood = function(tmle_task, node) {
      
      if(node%in%self$updater$update_nodes){
        
        likelihood_factor <- self$factor_list[[node]]
        # first check for cached values for this task
        value_step <- self$cache$get_update_step(likelihood_factor, tmle_task)
        
        if(!is.null(value_step)){
          # if some are available, grab them
          likelihood_values <- self$cache$get_values(likelihood_factor, tmle_task)
        } else {
          # if not, generate new ones 
          likelihood_values <- self$initial_likelihood$get_likelihood(tmle_task, node)
          value_step <- 0
        }
        
        # apply updates if necessary
        # todo: maybe let updater handle this logic
        # think about what happens if we actually need an *older* likelihood value
        if(value_step<self$updater$step_number){
          updates <- seq(from=value_step+1, to=self$updater$step_number)
          epsilons <- self$updater$epsilons[updates]
          likelihood_values <- self$updater$apply_updates(tmle_task, self, likelihood_values, epsilons)
          value_step <- self$updater$step_number
        }
        
        # todo: this sets values even if we haven't changed anything
        self$cache$set_values(likelihood_factor, tmle_task, value_step, likelihood_values)

        
      } else {
        # not a node that updates, so we can just use initials
        likelihood_values <- self$initial_likelihood$get_likelihood(tmle_task, node)
        
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
    initial_likelihood = function() {
      return(private$.initial_likelihood)
    },
    updater = function() {
      return(private$.updater)
    },
    factor_list = function(){
      return(self$initial_likelihood$factor_list)
    },
    training_task = function(){
      return(self$initial_likelihood$training_task)
    }
  ),
  private = list(
    .initial_likelihood = NULL,
    .updater = NULL
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
