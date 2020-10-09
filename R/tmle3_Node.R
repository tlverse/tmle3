#' A Node (set of variables) in an NPSEM
#'
#' This class defines a node in an NPSEM
#'
#' @importFrom R6 R6Class
#' @export
#'
#' @keywords data
#'
#' @return \code{tmle3_Node} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @section Constructor:
#'   \code{make_tmle3_task(name, variables, parents = c(), variable_type = NULL)}
#'
#'   \describe{
#'     \item{\code{name}}{character, the name of node
#'     }
#'     \item{\code{variables}}{character vector, the names of the variables that comprise the node
#'     }
#'     \item{\code{parents}}{character vector, the names of the parent nodes. If censoring, node is assumed to have no parents.
#'     }
#'     \item{\code{variable_type}}{\code{\link[sl3]{variable_type}} object, specifying the data type of this variable.
#'     If censoring, variable_type will be guessed later from the data.
#'     }
#'     }
#'
#' @section Methods:
#'
#' \describe{
#' \item{\code{guess_variable_type(variable_data)}}{
#'   Guesses the \code{\link[sl3]{variable_type}} from the provided data.
#'   This will be called by the \code{\link{tmle3_Task}} constructor if no variable_type was provided.
#'
#'   \itemize{
#'     \item{\code{variable_data}: the observed variable data.
#'     }
#'   }
#'   }
#' }
#'
#' @section Fields:
#' \describe{
#'     \item{\code{name}}{character, the name of node
#'     }
#'     \item{\code{variables}}{character vector, the names of the variables that comprise the node
#'     }
#'     \item{\code{parents}}{character vector, the names of the parent nodes. If censoring, node is assumed to have no parents.
#'     }
#'     \item{\code{variable_type}}{\code{\link[sl3]{variable_type}} object, specifying the data type of this variable.}
#' }
tmle3_Node <- R6Class(
  classname = "tmle3_Node",
  portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(name, variables, parents = c(),
                          variable_type = NULL, censoring_node = NULL, scale = FALSE) {
      private$.name <- name
      private$.variables <- variables
      private$.parents <- parents
      private$.variable_type <- variable_type
      private$.scale <- scale
      private$.censoring_node <- censoring_node
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
    censoring_node = function(new_censoring_node = NULL) {
      if (!is.null(new_censoring_node)) {
        private$.censoring_node <- new_censoring_node
      }
      return(private$.censoring_node)
    },
    parents = function() {
      return(private$.parents)
    },
    scale = function() {
      return(private$.scale)
    },
    variable_type = function(new_variable_type = NULL) {
      if (!is.null(new_variable_type)) {
        private$.variable_type <- new_variable_type
      }
      return(private$.variable_type)
    }
  ),
  private = list(
    .name = NULL,
    .variables = NULL,
    .censoring_node = NULL,
    .parents = NULL,
    .variable_type = NULL,
    .scale = NULL
  )
)

#' @param ... Passes all arguments to the constructor. See documentation for the
#'  Constructor below.
#' @export
#' @rdname tmle3_Node
define_node <- tmle3_Node$new

#' Helper functions for the NPSEM
#'
#' \code{all_ancestors} returns a list of all_ancestors of the specified node.
#' \code{time_ordering} attempts to find a time_ordering for the variables.
#'
#' @param node_name the node to search for ancestors of
#' @param npsem the NPSEM, defined by a list of \code{\link{tmle3_Node}} objects.
#' @rdname npsem_helpers
#' @export
all_ancestors <- function(node_name, npsem) {
  this_node <- npsem[[node_name]]
  if (length(this_node$parents) > 0) {
    ancestors <- unlist(lapply(this_node$parents, all_ancestors, npsem))
    ancestors <- unique(c(this_node$parents, ancestors))
    return(ancestors)
  } else {
    return(NULL)
  }
}

#' Get time ordering of nodes
#'
#' @rdname npsem_helpers
#' @export
time_ordering <- function(npsem) {
  node_names <- lapply(npsem, `[[`, "name")
  node_parents <- lapply(npsem, `[[`, "parents")
  n_nodes <- length(npsem)
  ordered_nodes <- c()
  while (length(ordered_nodes) < n_nodes) {
    parents_ordered <- sapply(node_parents, function(parents) all(parents %in% ordered_nodes))
    new_ordered <- setdiff(unlist(node_names[parents_ordered]), ordered_nodes)
    ordered_nodes <- c(ordered_nodes, new_ordered)
    if (length(new_ordered) == 0) {
      stop("failed to order nodes")
    }
  }

  return(ordered_nodes)
}
