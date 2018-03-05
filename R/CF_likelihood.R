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

      if (inherits(intervention_list, "LF_base")) {
        intervention_list <- list(intervention_list)
      }

      intervention_nodes <- sapply(intervention_list, `[[`, "name")
      names(intervention_list) <- intervention_nodes

      private$.intervention_list <- intervention_list
      params <- args_to_list()
      super$initialize(params)
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
    factor_list = function() {
      fl <- self$observed_likelihood$factor_list
      fl[names(self$intervention_list)] <- self$intervention_list[names(self$intervention_list)]
      return(fl)
    },
    update_list = function() {
      self$observed_likelihood$update_list
    }
  ),
  private = list(
    .observed_likelihood = NULL,
    .intervention_list = NULL
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
