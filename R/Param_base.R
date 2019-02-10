#' Base Class for Defining Parameters
#'
#' A parameter is a function of the likelihood. Once given a \code{\link{Likelihood}} object, a parameter will a value.
#' These objects also contain information about the efficient influence function (EIF) of a parameter, as well as its clever covariate(s).
#'
#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @family Parameters
#' @keywords data
#'
#' @return \code{Param_base} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @template Param_base_extra
#' @export
Param_base <- R6Class(
  classname = "Param_base",
  portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(observed_likelihood, ..., outcome_node = "Y") {
      private$.observed_likelihood <- observed_likelihood
      private$.outcome_node <- outcome_node
    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full") {
      stop("Param_base is a base class")
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      stop("Param_base is a base class")
    }
  ),
  active = list(
    type = function() {
      return(private$.type)
    },
    observed_likelihood = function() {
      return(private$.observed_likelihood)
    },
    outcome_node = function() {
      return(private$.outcome_node)
    }
  ),
  private = list(
    .type = "undefined",
    .observed_likelihood = NULL,
    .outcome_node = NULL
  )
)


#' Define a Parameter
#'
#' @param Param_class the class of the Parameter. Should inherit from \code{\link{Param_base}}
#' @param ... arguments that define the parameter See the constructor for the specified \code{Parameter}.
#' @family Parameters
#' @export
#
define_param <- function(Param_class, ...) {
  return(Param_class$new(...))
}
