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

      if (!is.null(observed_likelihood$censoring_nodes[[outcome_node]])) {
        if (!self$supports_outcome_censoring) {
          stop(sprintf(
            "%s has censoring mechanism, but %s does not yet support outcome node censoring",
            outcome_node, class(self)[1]
          ))
        }
      }

      if (inherits(observed_likelihood, "Targeted_Likelihood")) {
        # register parameter with updater
        observed_likelihood$updater$register_param(self)
      } else if (inherits(observed_likelihood, "Likelihood")) {
        warning("Parameter was passed a non-Targeted Likelihood object so estimates cannot be updated from initial")
      } else {
        stop("Invalid Likelihood class: ", class(observed_likelihood))
      }
    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full", submodel_type = "logistic") {
      # TODO Returns clever covariates and component-specific ICs ?
      stop("Param_base is a base class")
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      #Returns full IC and estimate
      stop("Param_base is a base class")
    },

    print = function() {
      cat(sprintf("%s: %s\n", class(self)[1], self$name))
    },
    supports_submodel_type = function(submodel_type, node){

      if(is.list(self$submodel_type_supported)) {
        if( !is.null(self$submodel_type_supported[[node]]) & self$submodel_type_supported[[node]]!=submodel_type) {
          stop(sprintf("This Param does not support the optimization strategy: %s", submodel_type))

        }
      }
      else if(!(submodel_type %in% self$submodel_type_supported)){
        stop(sprintf("This Param does not support the optimization strategy: %s", submodel_type))
      }

    }
  ),
  active = list(
    name = function() {
      return(private$.type)
    },
    type = function() {
      return(private$.type)
    },

    observed_likelihood = function() {
      return(private$.observed_likelihood)
    },
    outcome_node = function() {
      return(private$.outcome_node)
    },
    supports_outcome_censoring = function() {
      return(private$.supports_outcome_censoring)
    },
    supports_weights = function() {
      return(TRUE)
    },
    targeted = function() {
      return(private$.targeted)
    },
  submodel_type_supported = function() {
    return(private$.submodel_type_supported)
  }
  ),
  private = list(
    .type = "undefined",
    .observed_likelihood = NULL,
    .outcome_node = NULL,
    .targeted = TRUE,
    .supports_outcome_censoring = FALSE,
    .submodel_type_supported = c("logistic"),
    .supports_weights = T
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
