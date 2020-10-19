#' Base Class for Defining Likelihood Factors
#'
#' A Likelihood factor models a conditional density function.
#' The conditioning set is defined as all parent nodes (defined in \code{\link{tmle3_Task}}). In the case of a continuous
#' outcome variable, where a full density isn't needed, this can also model a conditional mean. This is the base class, which
#' is intended to be abstract. See below for a list of possible likelihood factor classes.
#'
#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @family Likelihood objects
#' @keywords data
#'
#' @return \code{LF_base} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @template LF_base_extra
#'
#' @export
LF_base <- R6Class(
  classname = "LF_base",
  portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(name, bound = NULL, ..., type = "density", cache = TRUE) {
      private$.name <- name
      private$.type <- type
      private$.bound <- bound
      private$.uuid <- UUIDgenerate(use.time = TRUE)
      private$.cache <- cache
    },
    delayed_train = function(tmle_task) {
      return(list())
    },
    train = function(tmle_task, ...) {
      # get possible values from task if discrete
      tmle_node <- tmle_task$npsem[[self$name]]
      private$.variable_type <- tmle_node$variable_type
      private$.training_task <- tmle_task
      # subclasses may do more, like fit sl3 models
    },
    get_density = function(tmle_task, fold_number) {
      stop("density not supported")
    },
    get_mean = function(tmle_task, fold_number) {
      stop("mean not supported")
    },
    get_likelihood = function(tmle_task, fold_number = "full") {
      if (self$type == "mean") {
        values <- self$get_mean(tmle_task, fold_number)
      } else {
        values <- self$get_density(tmle_task, fold_number)
      }
      if (!is.null(self$bound)) {
        values <- bound(values, self$bound)
      }

      return(values)
    },
    sample = function(tmle_task, n_samples = NULL, fold_number = "full") {
      stop("sampling not supported")
    },
    cf_values = function(tmle_task) {
      stop(sprintf("%s is not a valid intervention type", class(self)[1]))
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
    },
    uuid = function() {
      return(private$.uuid)
    },
    bound = function() {
      return(private$.bound)
    },
    cache = function() {
      return(private$.cache)
    },
    training_task = function() {
      return(private$.training_task)
    }
  ),
  private = list(
    .name = NULL,
    .variable_type = c(),
    .memoized_values = list(),
    .type = NULL,
    .uuid = NULL,
    .bound = NULL,
    .cache = TRUE,
    .training_task = NULL
  )
)

#' Define a Likelihood Factor
#'
#' @param LF_class the class of likelihood factor. Should inherit from \code{\link{LF_base}}
#' @param ... arguments that define the likelihood factor. See the constructor for the specified \code{LF_class}.
#' @family Likelihood objects
#' @export
#
define_lf <- function(LF_class, ...) {
  return(LF_class$new(...))
}
