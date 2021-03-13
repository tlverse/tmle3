#' Delta Method Parameters
#'
#' These parameters are smooth functionals of one or more other params
#' They are not fit directly with tmle, but are estimated using the delta method
#' todo: better docs
#' They do not return have clever covariates
#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @family Parameters
#' @keywords data
#'
#' @export
Param_delta <- R6Class(
  classname = "Param_delta",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, delta_param, parent_parameters,
                          ..., outcome_node = NA) {
      super$initialize(observed_likelihood, ..., outcome_node = outcome_node)
      private$.delta_param <- delta_param
      private$.parent_parameters <- parent_parameters
    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full") {
      return(list())
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      estimates <- lapply(
        self$parent_parameters,
        function(tmle_param) {
          tmle_param$estimates(tmle_task, fold_number)
        }
      )

      psis <- lapply(estimates, `[[`, "psi")
      ICs <- lapply(estimates, `[[`, "IC")
      psi <- self$delta_param$f(x = psis, dx = ICs)
      IC <- self$delta_param$df(x = psis, dx = ICs)

      list(
        psi = psi, IC = IC, name = self$name,
        transform = self$delta_param$transform
      )
    }
  ),
  active = list(
    name = function() {
      param_names <- sapply(
        self$parent_parameters,
        function(tmle_param) tmle_param$name
      )
      name <- self$delta_param$name(param_names)
      return(name)
    },
    type = function() {
      return(self$delta_param$type)
    },
    delta_param = function() {
      return(private$.delta_param)
    },
    parent_parameters = function() {
      return(private$.parent_parameters)
    },
    update_nodes = function() {
      return(NULL)
    }
  ),
  private = list(
    .delta_param = NULL,
    .parent_parameters = list()
  )
)
