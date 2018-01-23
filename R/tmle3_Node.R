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

#' Find all ancestors of given node
#' 
#' @param node_name the node to search for ancestors of
#' @param tmle_nodes the list of nodes and parents (i.e. the graph)
all_ancestors <- function(node_name,tmle_nodes){
  this_node <- tmle_nodes[[node_name]]
  if(length(this_node$parents)>0){
    ancestors <- unlist(lapply(this_node$parents,all_ancestors, tmle_nodes))
    ancestors <- unique(c(this_node$parents, ancestors))
    return(ancestors)
  } else{
    return(NULL)
  }
}

#' Get time ordering of nodes
#' 
#' @param tmle_nodes the list of nodes and parents (i.e. the graph)
time_ordering <- function(tmle_nodes){
  node_names <- lapply(tmle_nodes, `[[`, "name")
  node_parents <- lapply(tmle_nodes, `[[`, "parents")
  n_nodes <- length(tmle_nodes)
  ordered_nodes <- c()
  while(length(ordered_nodes) < n_nodes){
    parents_ordered <- sapply(node_parents, function(parents)all(parents%in%ordered_nodes))
    new_ordered <- setdiff(unlist(node_names[parents_ordered]), ordered_nodes)
    ordered_nodes <- c(ordered_nodes, new_ordered)
    if(length(new_ordered)==0){
      stop("failed to order nodes")
    }
  }
  
  return(ordered_nodes)
}

