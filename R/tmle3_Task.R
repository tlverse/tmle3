#' Class for Storing Data and NPSEM for TMLE
#'
#' This class inherits from \code{\link[sl3]{sl3_Task}}. In addition to all the
#'  methods supported by \code{\link[sl3]{sl3_Task}}, it supports the following.
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
#
tmle3_Task <- R6Class(
  classname = "tmle3_Task",
  portable = TRUE,
  class = TRUE,
  inherit = sl3_Task,
  public = list(
    initialize = function(data, npsem, ...) {
      super$initialize(data, covariates = c(), outcome = NULL, ...)

      node_names <- sapply(npsem, `[[`, "name")
      names(npsem) <- node_names

      # process nodes
      for (node_name in node_names) {
        current_node <- npsem[[node_name]]

        # get variable data and censoring indicator
        variables <- current_node$variables

        if (length(variables) == 0) {
          next
        }
        variable_data <- super$get_data(, variables)
        censoring <- apply(is.na(variable_data), 1, any)

        if (ncol(variable_data) == 1) {
          variable_data <- unlist(variable_data, use.names = FALSE)
        }

        # determine variable type
        if (is.null(current_node$variable_type)) {
          current_node$guess_variable_type(variable_data)
        }

        current_type <- current_node$variable_type

        # setup bounds for scaling of bounded continuous outcome if necessary


        if ((current_node$scale) &&
          (current_type$type == "continuous") &&
          (is.null(current_type$bounds))) {
          min_x <- min(variable_data)
          max_x <- max(variable_data)
          range <- max_x - min_x
          lower <- min_x #- 0.1 * range
          upper <- max_x #+ 0.1 * range
          bounded_variable_type <- variable_type(
            type = "continuous",
            bounds = c(lower, upper)
          )
          current_node$variable_type <- bounded_variable_type
        }

        # create or identify censoring node
        if (any(censoring)) {

          # first, look for explicitly denoted censoring node
          censoring_node <- current_node$censoring_node

          # next look in the npsem with the naming convention delta_X

          if (is.null(censoring_node)) {
            censoring_node_name <- sprintf("delta_%s", current_node$name)
            censoring_node <- npsem[[censoring_node_name]]
          } else {
            censoring_node_name <- censoring_node$name
          }

          # if we can't find a node, create one automatically

          if (is.null(censoring_node)) {

            # add censoring indicator to data
            censoring_dt <- data.table(!censoring)
            names(censoring_dt) <- censoring_node_name
            new_column_names <- super$add_columns(censoring_dt, uuid::UUIDgenerate())
            private$.column_names <- new_column_names

            censoring_node <- tmle3_Node$new(
              name = censoring_node_name,
              variables = censoring_node_name,
              parents = current_node$parents,
              variable_type = variable_type("binomial"),
              censoring_node = NULL,
              scale = FALSE
            )
          }

          # add censoring node to npsem and to current node
          current_node$censoring_node <- censoring_node
          npsem[[censoring_node_name]] <- censoring_node
        } else {
          # do we want to delete missingness node here?
        }

        # update npsem
        npsem[[node_name]] <- current_node
      }

      private$.npsem <- npsem
      private$.node_cache <- new.env()
    },
    get_tmle_node = function(node_name, format = FALSE, impute_censoring = FALSE) {
      cache_key <- sprintf("%s_%s_%s", node_name, format, impute_censoring)

      cached_data <- get0(cache_key, private$.node_cache, inherits = FALSE)
      if (!is.null(cached_data)) {
        return(cached_data)
      }
      tmle_node <- self$npsem[[node_name]]
      node_var <- tmle_node$variables
      if (is.null(node_var)) {
        return(data.table(NULL))
      }

      node_type <- tmle_node$node_type
      data <- self$get_data(self$row_index, node_var)

      if ((ncol(data) == 1)) {
        data <- unlist(data, use.names = FALSE)
      }

      if (format == TRUE) {
        var_type <- tmle_node$variable_type
        data <- var_type$format(data)
        data <- self$scale(data, node_name)
        data <- data.table(data)
        setnames(data, node_var)
      }

      censoring_node <- tmle_node$censoring_node

      if (is(censoring_node, "tmle3_Node") && impute_censoring) {
        observed <- self$get_tmle_node(censoring_node$name)
        censoring <- !observed

        # impute arbitrary value for node Need to keep the data shape the same,
        # but value should not matter here as this will only be used for prediction
        # and for generating values for ICs (which will then be cancelled by 0)
        impute_value <- data[which(!censoring)[1]]
        if (is.data.table(data)) {
          set(data, which(censoring), names(data), as.list(impute_value))
        } else {
          data[censoring] <- impute_value
        }
      }



      assign(cache_key, data, private$.node_cache)

      return(data)
    },
    get_regression_task = function(target_node, scale = FALSE, drop_censored = FALSE) {
      npsem <- self$npsem
      target_node_object <- npsem[[target_node]]
      parent_names <- target_node_object$parents
      parent_nodes <- npsem[parent_names]

      outcome_data <- self$get_tmle_node(target_node, format = TRUE)
      all_covariate_data <- lapply(parent_names, self$get_tmle_node, format = TRUE)

      outcome <- target_node_object$variables
      covariates <- unlist(lapply(parent_nodes, `[[`, "variables"))



      nodes <- self$nodes
      node_data <- self$get_data(, unlist(nodes))
      nodes$outcome <- outcome
      nodes$covariates <- covariates


      regression_data <- do.call(cbind, c(all_covariate_data, outcome_data, node_data))

      censoring_node <- target_node_object$censoring_node

      if (is(censoring_node, "tmle3_Node")) {
        observed <- self$get_tmle_node(censoring_node$name)
        censoring <- !observed
      } else {
        censoring <- rep(FALSE, nrow(regression_data))
      }

      folds <- self$folds
      if (drop_censored) {
        regression_data <- regression_data[!censoring, ]
        folds <- sl3::subset_folds(self$folds, which(!censoring))
      } else {
        # impute arbitrary value for node Need to keep the data shape the same,
        # but value should not matter here as this will only be used for prediction
        # and for generating values for ICs (which will then be cancelled by 0)
        impute_value <- regression_data[which(!censoring)[1], outcome, with = FALSE]
        set(regression_data, which(censoring), outcome, impute_value)
      }

      suppressWarnings({
        regression_task <- sl3_Task$new(
          regression_data,
          nodes = nodes,
          outcome_type = target_node_object$variable_type,
          folds = folds
        )
      })
      return(regression_task)
    },
    generate_counterfactual_task = function(uuid, new_data) {
      # for current_factor, generate counterfactual values
      node_names <- names(new_data)
      node_variables <- sapply(
        node_names,
        function(node_name) {
          self$npsem[[node_name]]$variables
        }
      )
      setnames(new_data, node_names, node_variables)

      new_task <- self$clone()
      new_column_names <- new_task$add_columns(new_data, uuid)
      new_task$initialize(
        self$internal_data, self$npsem,
        column_names = new_column_names,
        folds = self$folds,
        row_index = self$row_index
      )
      return(new_task)
    },
    next_in_chain = function(...) {
      return(super$next_in_chain(npsem = self$npsem, ...))
    },
    print = function() {
      cat(sprintf("A sl3 Task with %d obs and these nodes:\n", self$nrow))
      print(self$npsem)
    },
    get_node_bounds = function(node) {
      npsem <- self$npsem
      node_object <- npsem[[node]]
      variable_type <- node_object$variable_type
      return(variable_type$bounds)
    },
    scale = function(x, node) {
      bounds <- self$get_node_bounds(node)

      # nothing to do if no bounds, so return untransformed
      if (is.null(bounds)) {
        return(x)
      }

      scale <- bounds[2] - bounds[1]
      shift <- bounds[1]
      x_scaled <- (x - shift) / scale

      return(x_scaled)
    },
    unscale = function(x_scaled, node) {
      bounds <- self$get_node_bounds(node)

      # nothing to do if no bounds, so return untransformed
      if (is.null(bounds)) {
        return(x_scaled)
      }

      scale <- bounds[2] - bounds[1]
      shift <- bounds[1]
      x <- (x_scaled * scale) + shift

      return(x)
    },
    subset_task = function(row_index, drop_folds = FALSE) {
      if (is.logical(row_index)) {
        row_index <- which(row_index)
      }
      old_row_index <- private$.row_index
      if (!is.null(old_row_index)) {
        # index into the logical rows of this task
        row_index <- old_row_index[row_index]
      }
      new_task <- self$clone()
      if (drop_folds) {
        new_folds <- NULL
      } else {
        new_folds <- sl3::subset_folds(self$folds, row_index)
      }

      new_task$initialize(
        self$internal_data, self$npsem,
        column_names = self$column_names,
        folds = new_folds,
        row_index = row_index
      )
      return(new_task)
    }
  ),
  active = list(
    npsem = function() {
      return(private$.npsem)
    },
    data = function() {
      all_variables <- unlist(lapply(self$npsem, `[[`, "variables"))
      self$get_data(columns = all_variables)
    }
  ),
  private = list(
    .npsem = NULL,
    .node_cache = NULL
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
