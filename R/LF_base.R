#' Base Class for Likelihood Functionals
#'
#' @importFrom R6 R6Class
#'
#' @export
#
LF_base <- R6Class(classname = "LF_base",
                     portable = TRUE,
                     class = TRUE,
  public = list(
    initialize <- function(name, ...) {
      private$.name <- name
    },
    get_likelihood <- function(task, only_observed = FALSE) {
      stop("this is a base class")
    },
    print <- function() {
      cat(sprintf("%s: %s\n", self$name, class(self)[1]))
    }
  ),
  active = list(
    name <- function() {
      return(private$.name)
    },
    is_degenerate <- function() {
      return(private$.is_degenerate)
    }
  ),
  private = list(
    .name = NULL,
    .is_degenerate = FALSE
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

