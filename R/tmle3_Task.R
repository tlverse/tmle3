#' Class for Storing Data and NPSEM for TMLE
#'
#' This class inherits from \code{\link[sl3]{sl3_Task}}. In addition to all the
#'  methods supported by \code{\link[sl3]{sl3_Task}}, it supports the following.
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#' @importFrom sl3 sl3_Task
#' @importFrom digest digest
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
    initialize = function(data, npsem,  intervene_all_at_risk = F, summary_measure_columns = NULL, ...) {

      if(!inherits(data, "Shared_Data")){
        # For ease of coding and cleanness of code (and working with data.tables)
        # I assume that the id and time columns are "id" and "t" respectively.

        if(!is.null(nodes)){
          id <- nodes$id
          time <- nodes$time
        }

        if(is.null(id)){
          set(data, , "id", 1:nrow(data))
        }
        if(is.null(time)){
          set(data, , "t", rep(0,nrow(data)))
        }

        if(id!="id"){
          data[,"id", with = F] <- data[,id, with = F]
          id = "id"
        }
        if(time!="t"){
          data[, "t", with = F] <- data[,time, with = F]
          time = "t"
        }
        data <- setkey(data, id, t)
        shared_data <- data
      } else{
        # This assumes preprocessing has been done (e.g. sorting by id and t)
        shared_data <- data
        if(key(shared_data$data) != c("id", "t")){
          stop("Shared_Data object passed does not have a (id, t) key set.")
        }
      }

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
        # TODO Should rethink this and how it will work with the risk_set_map
        # In principle,
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
      private$.uuid <- digest(self$data)
    },
    get_tmle_node = function(node_name, format = FALSE, impute_censoring = FALSE, include_time = F, include_id = F, force_time_value = NULL, expand = F, compute_risk_set = T) {
      if(private$.force_at_risk) expand <- T

      if(is.null(force_time_value)) force_time_value <- F
      cache_key <- sprintf("%s_%s_%s_%s", node_name, format, force_time_value, expand)

      cached_data <- get0(cache_key, private$.node_cache, inherits = FALSE)
      if (!is.null(cached_data)) {
        if(!include_time){
          cached_data$t <- NULL
        }
        if(!include_id){
          cached_data$id <- NULL
        }
        return(cached_data)
      }
      tmle_node <- self$npsem[[node_name]]
      node_var <- tmle_node$variables

      if (is.null(node_var)) {
        return(data.table(NULL))
      }



      if(is.numeric(force_time_value)){
        time <- force_time_value
      } else {
        time <- tmle_node$time
      }
      if(is.null(time)) time <- 0

      if(length(time) > 1){
        at_risk_map <- tmle_node$at_risk_map
        if(expand  | !is.null(tmle_node$at_risk_map) | !tmle_node$missing_not_at_risk){
          #TODO, when to get value at all times by repeeatedly calling get_tmle_node with force_time_value argument??
          # The main issue is that computing the at_risk indicator requires applying a function to data[t <= time]
          # so there isn't any general shortcut exploiting the long format of the data
          data <- lapply(time, self$get_tmle_node, node_name= node_name, format = format, include_time = T, include_id = T, expand = expand)
          #setkey(data, id , t)
          return(data)
        }
        else {
          data <-  self$get_data(self$row_index,  c("id", "t", node_var))
          data <- data[t %in% time]
        }

      } else {
        #The at_risk summary measure might need other columns so grab all
        data <-  self$get_data(self$row_index,)
        data <- data[t <= time]
        if(compute_risk_set & !private$.force_at_risk){
          risk_set <- tmle_node$risk_set(data, time)
        }
        data <- data[, c("id", "t", node_var), with = F]

        if(expand){
          # Get most recent value for all
          data <- data[, last(.SD), by = id]
          if(compute_risk_set){
            if(private$.force_at_risk) {
              data$at_risk <- 1
            } else {
              data$at_risk <- as.numeric(data$id %in% risk_set)
            }
            if(time > 0) {
              last_vals <- self$data[t <= (time - 1), last(.SD), by = id, .SDcols = c(node_var)][,c(node_var),with = F]
            } else {
              last_vals <- NA
            }
            set(data, , paste0("last_val_",node_var) , last_vals)
          }

        } else {
          # Get most recent value for all those at risk
          if(compute_risk_set){
            data <- data[id %in% risk_set, last(.SD), by = id]
          } else {
            data <- data[, last(.SD), by = id]
          }

        }

      }
      set(data,, "t", time)
      if (format == TRUE) {
        data_node <- data[, node_var, with = F]
        if ((ncol(data_node) == 1)) {
          data_node <- unlist(data_node, use.names = FALSE)
        }
        var_type <- tmle_node$variable_type
        data_node <- var_type$format(data_node)
        data_node <- self$scale(data_node, node_name)
        set(data,, node_var, data_node)
      }

      assign(cache_key, data, private$.node_cache)

      if(!include_time){
        data$t <- NULL
      }
      if(!include_id){
        data$id <- NULL
      }

      censoring_node <- tmle_node$censoring_node

      # TODO So I think we should treat censoring and risk_set's differently
      # We say someone is no longer at_risk if their value will not change in time
      # We say someone is censored if their value is truly missing/unobserved
      # For the case of Param_survival, if we define our nodes
      # as the observed censoring and failure counting processes:
      #N(t) = 1(Ttilde <=t, Delta = 1), A(t) = 1(Ttilde <=t, Delta = 0)
      #Then we do not actually have censoring/outcome missingness in the sense that
      # the observed data we need for estimation isn't actually missing.
      # In this case, the risk_set_map is what we need to describe the nodes
      # (i.e. when any of the counting processes jumps, then they both remain the same value with prob 1)
      # On the other hand, if we have (W, A Y) and Y is missing then this is true censoring
      # We are really missing the observed data Y.
      # And it does not make sense to say this individual is no longer at risk/their value of Y stays the same
      #So I think we still need the notion of a censoring node but just need to be careful how we use it.
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
    # TODO: add time_variance
    get_regression_task = function(target_node, scale = FALSE, drop_censored = FALSE, is_time_variant = FALSE) {
      npsem <- self$npsem
      target_node_object <- npsem[[target_node]]
      parent_names <- target_node_object$parents
      parent_nodes <- npsem[parent_names]

      outcome_data <- self$get_tmle_node(target_node, format = TRUE)
      all_covariate_data <- lapply(parent_names, self$get_tmle_node, format = TRUE)

      outcome <- target_node_object$variables
      # TODO: check
      cov_nodes <- parent_nodes
      covariates <- unlist(lapply(cov_nodes, `[[`, "variables"))



      nodes <- self$nodes
      node_data <- self$get_data(, unlist(nodes))
      nodes$outcome <- outcome
      nodes$covariates <- covariates


      regression_data <- do.call(cbind, c(all_covariate_data, outcome_data, node_data))

      if ((is_time_variant) && (!is.null(self$nodes$time))) {
        regression_data$time <- self$time
        nodes$covariates <- c(nodes$covariates, "time")
      }

      censoring_node <- target_node_object$censoring_node

      indices <- seq_len(self$nrow)
      if (is(censoring_node, "tmle3_Node")) {
        observed <- self$get_tmle_node(censoring_node$name)
        censoring <- !observed
      } else {
        censoring <- rep(FALSE, nrow(regression_data))
      }

      if (drop_censored) {
        indices <- intersect(indices, which(!censoring))
      } else {
        # impute arbitrary value for node Need to keep the data shape the same,
        # but value should not matter here as this will only be used for prediction
        # and for generating values for ICs (which will then be cancelled by 0)
        impute_value <- regression_data[which(!censoring)[1], outcome, with = FALSE]
        set(regression_data, which(censoring), outcome, impute_value)
      }



      if ((!is_time_variant) && (!is.null(self$nodes$time))) {
        time_data <- self$time
        indices <- which(time_data == 1)
        indices <- intersect(indices, which(time_data == 1))
      }

      folds <- self$folds
      if (length(indices) < self$nrow) {
        regression_data <- regression_data[indices, ]
        folds <- sl3::subset_folds(folds, indices)
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
        nodes = self$nodes,
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
