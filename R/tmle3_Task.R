#' Class for Storing Data and NPSEM for TMLE
#'
#' This class inherits from \code{\link[sl3]{sl3_Task}}. In addition to all the methods supported by \code{\link[sl3]{sl3_Task}}, it supports the following.
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#' @importFrom sl3 sl3_Task
#' @import data.table
#'
#' @export
#'
#' @keywords data
#'
#' @return \code{tmle3_Task} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @template tmle3_Task_extra
#'
#' @export

tmle3_Task <- R6Class(
  classname = "tmle3_Task",
  portable = TRUE,
  class = TRUE,
  inherit = sl3_Task,
  public = list(
    initialize = function(data, tmle_nodes, ...) {
      super$initialize(data, covariates = c(), outcome = NULL, ...)

      node_names <- sapply(tmle_nodes, `[[`, "name")
      names(tmle_nodes) <- node_names
      for (node_name in node_names) {
        variables <- tmle_nodes[[node_name]]$variables
        variable_data <- super$get_data(, variables)
        if (ncol(variable_data) == 1) {
          variable_data <- unlist(variable_data, use.names = FALSE)
        }
        if (is.null(tmle_nodes[[node_name]]$variable_type)) {
          tmle_nodes[[node_name]]$guess_variable_type(variable_data)
        }
      }
      private$.tmle_nodes <- tmle_nodes
    },
    get_tmle_node = function(node_name, bound = FALSE) {
      tmle_node <- self$tmle_nodes[[node_name]]
      node_var <- tmle_node$variables

      data <- self$get_data(, node_var)

      if (bound) {
        bounds <- tmle_node$variable_type$bounds
        if (!is.null(bounds)) {
          scale <- bounds[2] - bounds[1]
          shift <- bounds[1]
          data <- data[, lapply(.SD, function(vals) (vals - shift) / scale)]
        }
      }


      if (ncol(data) == 1) {
        return(unlist(data, use.names = FALSE))
      } else {
        return(data)
      }
    },
    get_regression_task = function(target_node) {
      tmle_nodes <- self$tmle_nodes
      target_node <- tmle_nodes[[target_node]]
      parent_names <- target_node$parents
      parent_nodes <- tmle_nodes[parent_names]
      outcome <- target_node$variables
      covariates <- unlist(lapply(parent_nodes, `[[`, "variables"))

      # todo: consider if self$data isn't a better option here
      data <- self$raw_data

      # bound continuous outcome if bounds are specified to variable_type
      variable_type <- target_node$variable_type
      column_names <- self$column_names
      if ((variable_type$type == "continuous") && (!is.na(variable_type$bounds))) {
        # todo: make quasibinomial, make more learners play nice with quasibinomial outcomes
        bounded_vals <- self$get_tmle_node(target_node$name, bound = TRUE)
        col_name <- sprintf("__%s_bounded", target_node$name)
        set(data, , col_name, bounded_vals)
        column_names[col_name] <- col_name
        outcome <- col_name
      }

      nodes <- self$nodes
      nodes$outcome <- outcome
      nodes$covariates <- covariates

      # todo: make sure folds transfer
      return(sl3_Task$new(
        data, nodes = nodes,
        outcome_type = variable_type,
        column_names = column_names
      ))
    },
    generate_counterfactual_task = function(uuid, new_data) {
      # for current_factor, generate counterfactual values
      node_names <- names(new_data)
      node_variables <- sapply(node_names, function(node_name) self$tmle_nodes[[node_name]]$variables)
      setnames(new_data, node_names, node_variables)

      new_task <- self$clone()
      new_column_names <- new_task$add_columns(uuid, new_data)
      new_task$initialize(
        self$raw_data, self$tmle_nodes,
        column_names = new_column_names
      )
      return(new_task)
    },
    next_in_chain = function(...) {
      return(super$next_in_chain(tmle_nodes = self$tmle_nodes, ...))
    },
    print = function() {
      cat(sprintf("A sl3 Task with %d obs and these nodes:\n", self$nrow))
      print(self$tmle_nodes)
    }
  ),
  active = list(
    tmle_nodes = function() {
      return(private$.tmle_nodes)
    },
    data = function() {
      all_variables <- unlist(lapply(self$tmle_nodes, `[[`, "variables"))
      self$get_data(columns = all_variables)
    }
  ),
  private = list(
    .tmle_nodes = NULL
  )
)

#' @param ... Passes all arguments to the constructor. See documentation for the
#'  Constructor below.
#'
#' @rdname tmle3_Task
#'
#' @export
#
make_tmle3_Task <- tmle3_Task$new
