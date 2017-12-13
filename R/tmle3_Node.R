#' Class for Nonparametric Structural Equation Models
#'
#' @importFrom R6 R6Class
#'
#' @export
#
tmle3_Node <- R6Class(
  classname = "tmle3_Node",
  portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(name, variables, parents = c(),
                          variable_type = NULL) {
      private$.name <- name
      private$.variables <- variables
      private$.parents <- parents
      private$.variable_type <- variable_type
    },
    print = function() {
      node_class <- class(self)[1]
      cat(sprintf("%s: %s\n", node_class, self$name))
      cat(sprintf("\tVariables: %s\n", paste(self$variables, collapse = ", ")))
      cat(sprintf("\tParents: %s\n", paste(self$parents, collapse = ", ")))
    },
    guess_variable_type = function(variable_data) {
      private$.variable_type <- variable_type(x = variable_data)
    }
  ),
  active = list(
    name = function() {
      return(private$.name)
    },
    variables = function() {
      return(private$.variables)
    },
    parents = function() {
      return(private$.parents)
    },
    variable_type = function() {
      return(private$.variable_type)
    }
  ),
  private = list(
    .name = NULL,
    .variables = NULL,
    .parents = NULL,
    .variable_type = NULL
  )
)

#' Define a Node (set of variables) in an NPSEM
#'
#' @param name A character, the name of node
#' @param variables A character vector of the variables that comprise the node
#' @param parents A character vector of the names of the parent nodes. If missing, node is assumed to have no parents.
#' @param variable_type A \code{sl3} \code{variable_type} object specifying the data type of this variable.
#' If missing, variable_type is guessed from the data
#' @export
define_node <- function(name, variables, parents=c(), variable_type = NULL) {
  tmle3_Node$new(name, variables, parents, variable_type)
}
