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
    initialize = function(name, variables, parents = c(), time = NULL, summary_functions = NULL, risk_set_map = NULL, degeneracy_type = "last", missing_row_implies_not_at_risk = T,
                              variable_type = NULL, censoring_node = NULL, scale = FALSE) {


      if(is.null(time)){
        time <- 0
      }
      if(!is.null(summary_functions) & !is.list(summary_functions)){
        summary_functions <- list(summary_functions)
      }
      if (!is.null(risk_set_map)){
        if(!is.character(risk_set_map)){
          risk_set_map$set_name(paste(paste(variables, collapse = "_"), "at_risk", sep = "_") )
        }
      }
      private$.ltmle_params <- list(degeneracy_type = degeneracy_type, missing_row_implies_not_at_risk = missing_row_implies_not_at_risk, risk_set_map = risk_set_map, time = time, summary_functions = summary_functions)

      private$.name <- name
      private$.variables <- variables
      private$.parents <- parents
      private$.variable_type <- variable_type
      private$.scale <- scale
      private$.censoring_node <- censoring_node
    },
    print = function() {
      if(is.character(self$risk_set_map)){
        risk_name <- self$risk_set_map
      } else {
        risk_name <- self$risk_set_map$name
      }
      node_class <- class(self)[1]
      cat(sprintf("%s: %s\n", node_class, self$name))
      cat(sprintf("\tVariables: %s\n", paste(self$variables, collapse = ", ")))
      cat(sprintf("\tParents: %s\n", paste(self$parents, collapse = ", ")))
      cat(sprintf("\tTime: %s\n", paste(self$time, collapse = ", ")))
      cat(sprintf("\tSummary Measures: %s\n", paste(unlist(sapply(self$summary_functions, function(f){f$name})), collapse = ", ")))
      cat(sprintf("\tRisk-set Map: %s\n",risk_name))

    },
    guess_variable_type = function(variable_data) {
      private$.variable_type <- sl3::variable_type(x = variable_data)
    },
    risk_set = function(data, time, subset_time = T){

      if(subset_time) data <- data[t <= time]
      #Assumes data == data[t <= time,] and time is single number
      #Computes, for this node, the id's of those in data at risk of changing their value at this time
      at_risk_map <- self$risk_set_map
      missing_not_at_risk <- private$.ltmle_params$missing_row_implies_not_at_risk
      if(missing_not_at_risk){
        keep_id <- unique(data[t == time, id])
        data <- data[id %in% keep_id,]
      }
      if(is.null(at_risk_map)) {
        risk_set <- unique(data$id)
        return(risk_set)
        }
      if(is.character(at_risk_map)) {
        if(missing_not_at_risk){
          #If those missing rows are not at risk
          #then only check for those with rows at this time
          risk_set <- data[t == time & data[,at_risk_map,with = F, drop = T]==1, "id", with = F][["id"]]
        } else{
          #Otherwise find the last value of risk indicator
          data <- data[, last(.SD), by = id, .SDcols = at_risk_map]
          risk_set <- data$id[data[[at_risk_map]] ==1]
        }
      } else if (inherits(at_risk_map, "Summary_measure")){

        risk_set <-at_risk_map$summarize(data,time)# [,at_risk_map$name, with = F, drop = T]

        risk_set <- risk_set$id[risk_set[[at_risk_map$name]]==1]
      } else {
        risk_set <- risk_set(data, time)
      }
      return(unlist(risk_set, use.names = F))

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
    summary_functions = function(){
      return(private$.ltmle_params$summary_functions)
    },
    missing_not_at_risk = function(){
      private$.ltmle_params$missing_row_implies_not_at_risk
    },
    time = function(time = NULL){
      if(!is.null(time)){
        private$.ltmle_params$time <- time
      }
      private$.ltmle_params$time
    },
    risk_set_map = function(){
      private$.ltmle_params$risk_set_map
    },
    degeneracy_type = function(){
      private$.ltmle_params$degeneracy_type
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
    .scale = NULL,
    .ltmle_params = NULL
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
