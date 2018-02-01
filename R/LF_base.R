#' Base Class for Defining Likelihood Factors
#'
#' A Likelihood factor models a conditional density (continuous variable)
#' or a conditional probability (discrete variable) function for a particular tmle node.
#' The conditioning set is defined as all parent nodes (defined in tmle_task). In the case of a continuous
#' outcome variable, where a full density isn't needed, this can also model a conditional mean.
#'
#' @importFrom R6 R6Class
#'
#' @export
#
LF_base <- R6Class(
  classname = "LF_base",
  portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(name, type = "density", ...) {
      private$.name <- name
      private$.type <- type
    },
    delayed_train = function(tmle_task){
      return(list())
    },
    train = function(tmle_task, ...) {
      # get possible values from task if discrete
      tmle_node <- tmle_task$tmle_nodes[[self$name]]
      private$.variable_type <- tmle_node$variable_type

      # subclasses may do more, like fit sl3 models
    },
    get_density = function(tmle_task, only_observed = FALSE) {
      stop("this is a base class")
    },
    get_mean = function(task) {
      stop("this is a base class")
    },
    print = function() {
      cat(sprintf("%s: %s\n", self$name, class(self)[1]))
    }
  ),
  active = list(
    name = function() {
      return(private$.name)
    },
    variable_type = function() {
      return(private$.variable_type)
    },
    type = function() {
      return(private$.type)
    },
    values = function() {
      variable_type <- self$variable_type
      if (!is.null(variable_type)) {
        return(variable_type$levels)
      } else {
        return(NULL)
      }
    }
  ),
  private = list(
    .name = NULL,
    .variable_type = c(),
    .type = NULL
  )
)

#' Defining Likelihood Functionals
#'
#' @param LF_class the class of likelihood
#' @param ... arguments that define the likelihood functional
#'
#' @export
#
define_lf <- function(LF_class, ...) {
  return(LF_class$new(...))
}
