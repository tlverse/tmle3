#' Class for Storing Data and NPSEM for TMLE
#'
#' This class inherits from \code{\link[sl3]{sl3_Task}}. In addition to all the
#'  methods supported by \code{\link[sl3]{sl3_Task}}, it supports the following.
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#' @importFrom sl3 sl3_Task
#' @importFrom sl3 Shared_Data
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
    initialize = function(data, npsem, summary_measure_columns = NULL, id = NULL, time = NULL, force_at_risk = F, long_format = NULL, folds_for_ids = NULL, ...) {


      dot_args <- list(...)
      if(is.null(id)){
        id <- dot_args$nodes$id
      }
      if(is.null(time)){
        time <- dot_args$nodes$time
      }


      if(!inherits(data, "Shared_Data")){
        # For ease of coding and cleanness of code (and working with data.tables)
        # I assume that the id and time columns are "id" and "t" respectively.

        if(!is.null(long_format)){
          long_format <- long_format
        }
        else if(!is.null(time) | "t" %in% colnames(data)) {
          long_format <- T
        } else {
          long_format <- F
        }


        if(!is.data.table(data)) data <- as.data.table(data)
        #TODO if passed through nodes arg

        if(is.null(id)){
          if("id" %in% colnames(data)){
            warning("id argument not specified but id column found in data. Using this as id.")
            id <- "id"
          } else {
            set(data, , "id", 1:nrow(data))
            id <- "id"
          }

        }
        if(is.null(time)){
          if("t" %in% colnames(data)){
            warning("time argument not specified but column labeled t found in data. Using this as time")
            time <- "t"
          } else {
            set(data, , "t", rep(0,nrow(data)))
            time <- "t"
          }

        }

        if(id!="id"){
          data$id <- data[,id, with = F]
          id <- "id"
        }
        if(time!="t"){
          data$t <- data[,time, with = F]
          time <- "t"
        }
        # Ensure sorting is preserved if id is cast to factor later down line
        #TODO think about this
        data[, id := as.factor(id)]
        data <- setkey(data, id, t)
        shared_data <- data


      } else{
        # This assumes preprocessing has been done (e.g. sorting by id and t)
        shared_data <- data
        id <- "id"
        time <- "t"
        if(length(unique(shared_data$raw_data$t)) > 1){
          long_format <- T
        } else {
          long_format <- F
        }
        if(!all(key(shared_data$raw_data) == c("id", "t"))){
          setkey(shared_data$raw_data, "id", "t")
          stop("Shared_Data object passed does not have a (id, t) key set.")
        }
      }

      super$initialize(shared_data, covariates = c(), outcome = NULL, id = id, time  = time,  ...)
      if(!is.null(folds_for_ids)){
        private$.folds <- folds_for_ids
      }
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
              time = current_node$time,
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
      private$.force_at_risk <- force_at_risk
      private$.summary_measure_columns <- summary_measure_columns
      private$.long_format <- long_format
      if(private$.force_at_risk){
        add <- "T"
      } else{
        add <- ""
      }

      private$.uuid <- paste0(add, digest(self$data))

    },
    get_tmle_node = function(node_name, format = FALSE, impute_censoring = FALSE, include_time = F, include_id = F, force_time_value = NULL, expand = T, compute_risk_set = F) {
      force_at_risk <- private$.force_at_risk


      if(is.null(force_time_value)) force_time_value <- F
      cache_key <- sprintf("%s_%s_%s_%s_%s_%s_%s", node_name, format, impute_censoring, force_time_value, expand, compute_risk_set, force_at_risk)

      cached_data <- get0(cache_key, private$.node_cache, inherits = FALSE)
      if (!is.null(cached_data)) {
        if(!include_time){
          cached_data <- cached_data[, setdiff( names(cached_data), c("t")), with = F]
        }
        if(!include_id){
          cached_data <- cached_data[, setdiff( names(cached_data), c("id")), with = F]
        }
        if (!format & (ncol(cached_data) == 1)) {
          cached_data <- unlist(cached_data, use.names = FALSE)
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
          data <- lapply(time, function(t) self$get_tmle_node( force_time_value = t,node_name= node_name, format = format, include_time = T, include_id = T, expand = expand, compute_risk_set = compute_risk_set))
          data <- rbindlist(data)
          #Needed to ensure order matches regrression task outcome
          setkey(data, id , t)
          if(!include_time){
            data <- data[, setdiff( names(data), c("t")), with = F]
          }
          if(!include_id){
            data <- data[, setdiff( names(data), c("id")), with = F]
          }
          if(ncol(data) == 1 & !format) data <- unlist(data, use.names = F)
          return(data)
        }
        else {
          data <-  self$get_data(self$row_index,  c("id", "t", node_var))
          data <- data[t %in% time]

        }

      } else {
        #The at_risk summary measure might need other columns so grab all
        data <-  self$data
        data <- data[t <= time, ]
        if(compute_risk_set & !private$.force_at_risk){

          risk_set <- tmle_node$risk_set(data, time, subset_time = F)

        }
        data <- data[, c("id", "t", node_var), with = F]

        if(expand){
          # Get most recent value for all
          data <- data[, last(.SD), by = id]
          if(compute_risk_set){

            if(private$.force_at_risk) {
              set(data, , "at_risk", 1)
            } else {
              set(data, , "at_risk", as.numeric(data$id %in% risk_set))
            }

            if(tmle_node$degeneracy_type == "last" & time > 0){
              degeneracy_value <- self$data[t < time , last(.SD), by = id, .SDcols = c(node_var)][,c(node_var),with = F]
            } else if(is.numeric(tmle_node$degeneracy_type) & length(tmle_node$degeneracy_type)==1){
              degeneracy_value <- tmle_node$degeneracy_type
            } else {
              if(time > 0) {
                degeneracy_value <- self$data[t < time , last(.SD), by = id, .SDcols = c(node_var)][,c(node_var),with = F]
              } else {
                degeneracy_value <- NA
                }
            }

            set(data, , paste0("degeneracy_value_",node_var) , degeneracy_value)

          }

        } else {
          # Get most recent value for all those at risk
          if(compute_risk_set){
            data <- data[id %in% risk_set, last(.SD), by = id]
          } else {
            data <- data[, last(.SD), by = id]
          }

        }
        set(data,, "t", time)
      }

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



      censoring_node <- tmle_node$censoring_node


      if (is(censoring_node, "tmle3_Node") && impute_censoring) {
        observed <- self$get_tmle_node(censoring_node$name)
        censoring <- !observed

        # impute arbitrary value for node Need to keep the data shape the same,
        # but value should not matter here as this will only be used for prediction
        # and for generating values for ICs (which will then be cancelled by 0)
        impute_value <- 0
        if (is.data.table(data)) {
          set(data, which(censoring), node_var, impute_value)
        } else {
          data[censoring] <- impute_value
        }
      }

      assign(cache_key, data, private$.node_cache)

      if(!include_time){
        data <- data[, setdiff( names(data), c("t")), with = F]
      }
      if(!include_id){
        data <- data[, setdiff( names(data), c("id")), with = F]
      }


      if (!format & (ncol(data) == 1)) {
        data <- unlist(data, use.names = FALSE)
      }

      return(data)
    },
    # TODO: add time_variance
    get_regression_task = function(target_node, scale = FALSE, drop_censored = FALSE, is_time_variant = FALSE,  force_time_value = NULL, expand = T, cache_task = T, include_bins = F, bin_num = NULL) {
      long_format <- private$.long_format
      if(!is.numeric(force_time_value) & cache_task){
        cache_key <- sprintf("%s_%s_%s_%s_%s_%s", paste0(target_node, collapse = "%"), scale, drop_censored, is_time_variant, expand, self$force_at_risk)

        cached_data <- get0(cache_key, private$.node_cache, inherits = FALSE)
        if (!is.null(cached_data)) {
          return(cached_data)
        }
      }
      if(length(target_node)>1){
        all_tasks <- lapply(seq_along(target_node), function(i){
          self$get_regression_task(target_node[[i]] , scale, drop_censored , is_time_variant, expand = expand, include_bins = include_bins, bin_num = i)
        })
        all_nodes <- lapply(all_tasks, function(task) task$nodes)
        # Checks that names are compatible
        regression_X <- rbindlist(lapply(seq_along(all_tasks), function(i) {

          task <- all_tasks[[i]]
          #TODO all this sorting to ensure matching regression tasks for pooling seems inconvenient.
          vars <- sort(unique(c(task$nodes$covariates, task$nodes$id, task$nodes$t, task$nodes$weights)))

          data <- task$get_data(,vars)

          return(data)
          }), use.names = T)

        outcomes <- list()
        regression_Y <- rbindlist(lapply(seq_along(all_tasks), function(i) {
          task <- all_tasks[[i]]
          Y <- data.table(task$Y)
          outcomes[[i]] <<- task$nodes$outcome
          return(Y)
        }))
        outcomes <- transpose(outcomes)

        new_names <- unlist(sapply(outcomes, function(outcome) {paste(unique(outcome), collapse = "%")}))
        setnames(regression_Y, new_names)

        regression_data <- cbind(regression_X, regression_Y)
        setkey(regression_data, id, t)
        nodes <- all_nodes[[1]]

        nodes$outcome <-  colnames(regression_Y)
        # Make sure time is included as covariate
        if(is_time_variant){
          nodes$covariates <- union("t", nodes$covariates)
        }
        if(include_bins) {
          #Otherwise add bin number
          nodes$covariates <- union(nodes$covariates, "bin_num")
        }

        if(long_format){
          nodes$covariates <- sort(nodes$covariates)
          setcolorder(regression_data, order(colnames(regression_data)))

        }


        folds <- self$folds
        folds <- origami::id_folds_to_folds(folds, match(regression_data$id, unique(self$id)))
        if (FALSE & nrow(regression_data) < self$nrow) {
          data_id_t <- self$data[, c("id", "t"), with = F]
          indices <- data_id_t[regression_data[, c("id", "t"), with = F],  which =  T]
          folds <- sl3::subset_folds(folds, indices)
        }

        pooled_regression_task <- sl3_Task$new(
          regression_data,
          nodes = nodes,
          outcome_type = self$npsem[[target_node[1]]]$variable_type,
          folds = folds
        )

        if(!is.numeric(force_time_value)){
          assign(cache_key, pooled_regression_task, private$.node_cache)
        }
        return(pooled_regression_task)

      }


      npsem <- self$npsem
      target_node_object <- npsem[[target_node]]
      target_node <- target_node_object$name
      outcome <- target_node_object$variables
      if(is.numeric(force_time_value)){
        time <- force_time_value
      } else {
        time <- target_node_object$time
      }

      past_data <- self$data

      if(length(time) > 1){
        # TODO summary measures ae expensive to compute. The task cache helps.

        # If node is pooled across time then get pooled regression task

        all_tasks <- lapply(time, function(t) self$get_regression_task(force_time_value = t, target_node = target_node, scale = scale, drop_censored = drop_censored, is_time_variant = is_time_variant, expand = expand, include_bins = include_bins, bin_num = bin_num ))
        all_nodes <- lapply(all_tasks, function(task) task$nodes)
        regression_data <- rbindlist(lapply(all_tasks, function(task) task$get_data()), use.names = T)
        nodes <- all_nodes[[1]]
        if(is_time_variant) {
          nodes$covariates <- union("t", nodes$covariates)
        }
        if(long_format){
          setcolorder(regression_data, order(colnames(regression_data)))
          nodes$covariates <- sort(nodes$covariates)
        }

        setkey(regression_data, id, t)

        folds <- self$folds
        folds <- origami::id_folds_to_folds(folds, match(regression_data$id, unique(self$id)))
        if (FALSE & nrow(regression_data) < self$nrow) {
          data_id_t <- self$data[, c("id", "t"), with = F]
          indices <- data_id_t[regression_data[, c("id", "t"), with = F],  which =  T]
          folds <- sl3::subset_folds(folds, indices)
        }

        pooled_regression_task <- sl3_Task$new(
          regression_data,
          nodes = nodes,
          outcome_type = self$npsem[[target_node[1]]]$variable_type,
          folds = folds
        )

        if(!is.numeric(force_time_value)){
          #Store tasks
          assign(cache_key, pooled_regression_task, private$.node_cache)
        }
        return(pooled_regression_task)
      }

      parent_names <- target_node_object$parents
      #Ensure that output is independent of order of parent_names
      parent_names <- sort(parent_names)
      parent_nodes <- npsem[parent_names]

      if(is.null(unlist(target_node_object$summary_functions))){
        # No summary functions so simply stack node values of parents
        outcome_data <- self$get_tmle_node(target_node, format = TRUE, include_id = T, include_time = T, force_time_value = force_time_value, expand = T, compute_risk_set = T)
        if(length(outcome_data)==0){
          suppressWarnings({
            regression_task <- sl3_Task$new(
              outcome_data,
              covariates = c(),
              outcome_type = target_node_object$variable_type

            )
          })
          return(regression_task)
        }
        if(length(parent_names) > 0){
          parent_times <- lapply(parent_nodes, `[[`, "time")
          parent_data <-   lapply(parent_names, self$get_tmle_node, include_id = F, include_time = F, format = T, expand = T, compute_risk_set = F) #%>% purrr::reduce(merge, "id")
          parent_data <- setDT(unlist(parent_data, recursive = F))


          # This should ensure column names are unique
          if(include_bins){
            safe_names <- paste0(colnames(parent_data), "%", seq_along(colnames(parent_data)))
            setnames(parent_data, safe_names)#paste0(colnames(parent_data), "%", seq_along(colnames(parent_data))))
            # TODO this should no longer be needed
            setnames(parent_data, make.unique(colnames(parent_data)))
          }

          if(any(duplicated(colnames(parent_data)))) {

            warning("Duplicate names found in parent data. Making unique.")
            safe_names <- paste0(colnames(parent_data), "f", seq_along(colnames(parent_data)))
            setnames(parent_data, safe_names)#paste0(colnames(parent_data), "%", seq_along(colnames(parent_data))))
            # TODO this should no longer be needed
            setnames(parent_data, make.unique(colnames(parent_data)))
          }
          #By adding times to column names they become unique

        } else {
          parent_data <- data.table(NULL)
        }

        covariates <- colnames(parent_data) #unlist(lapply(parent_nodes, `[[`, "variables"))
        outcome = setdiff(colnames(outcome_data), c("id", "t", grep("degeneracy_value", colnames(outcome_data), value = T), "at_risk"))

        outcome_index <-  1:length(outcome)
        if(length(covariates)>0){

          cov_index <- length(outcome) + 1:length(covariates)
        } else {
          cov_index <- c()
        }
        #Due to time indexing, we do not have unique column names.
        #In order to support pooling across time, we shouldn't use node names as column names
        #important that outcome variable name doesnt change

        uniq_names <- make.unique(c(outcome,covariates))
        covariates <- uniq_names[cov_index]
        outcome <- uniq_names[outcome_index]
        setnames(parent_data, covariates)

        if((length(time) >1)){
          covariates <- c(covariates, "t")
        }
        all_covariate_data <- parent_data


      } else {


        all_vars <- unique(unlist(lapply(npsem, `[[`, "variables")))

        times <- as.vector(sapply(parent_nodes, function(node) node$time))
        parent_covariates <- as.vector(sapply(parent_nodes, function(node) node$variables))

        # Note that those with missing rows will be included in outcome_data.
        # There value will be set to last measured value.
        outcome_data <- self$get_tmle_node(target_node, format = TRUE, include_id = T, include_time = (time == "pooled"), force_time_value = force_time_value, expand = T, compute_risk_set = T)
        if(length(outcome_data)==0){
          suppressWarnings({
            regression_task <- sl3_Task$new(
              outcome_data,
              covariates = c(),
              outcome_type = target_node_object$variable_type

            )
          })
          return(regression_task)
        }
        past_data <- past_data[t <= time & id %in% outcome_data$id,]


        if(length(parent_covariates) != 0 | !is.null(unlist(target_node_object$summary_functions))){
          summary_measures <- target_node_object$summary_functions
          all_covariate_data <- lapply(summary_measures, function(fun){
            return(fun$summarize(past_data, time, add_id = F))}
          )
          #all_covariate_data <- all_covariate_data %>% purrr::reduce(merge, by = "id")
          #all_covariate_data <- do.call(cbind, all_covariate_data)
          all_covariate_data <- setDT(unlist(all_covariate_data, recursive = F))

          #set(all_covariate_data, , (which(duplicated(names(all_covariate_data)))), NULL)
          covariates <- setdiff(colnames(all_covariate_data), "id")
          if("t" %in% names(all_covariate_data))  set(all_covariate_data, "t", NULL)

        } else {
          all_covariate_data <-  data.table(NULL)

          covariates <- c()
        }
      }

      nodes <- self$nodes

      if(length(setdiff(nodes, c("id", "time"))) > 0){
        node_data <- self$get_data(self$row_index, unlist(nodes))
        #TODO since number of time rows vary per person, only time-indepdent nodes make sense
        # Keep only node_data for each individual at the time of this tmle node
        node_data <- node_data[node_data$id %in% outcome_data$id & node_data$t <= time, last(.SD), by = id]
        setorder(node_data, id)
        node_data$t <- NULL
        node_data$id <- NULL
      } else {
        node_data <- data.table(NULL)
      }

      nodes$outcome <- outcome
      nodes$covariates <- covariates
      if(ncol(all_covariate_data) == 0){
        regression_data <-  list(outcome_data, node_data)# %>% purrr::reduce(merge, "id")
        regression_data <- setDT(unlist(regression_data, recursive = F))
      }
      else {
        regression_data <-   list(all_covariate_data, outcome_data, node_data)
        regression_data <- setDT(unlist(regression_data, recursive = F))
           }
      if(!expand){
        regression_data <- regression_data[regression_data$at_risk ==1]
      }
      set(regression_data, , "t" , time)

      setkey(regression_data, id, t)
      # Necessary for pooled regression tasks and unpooled to be compatible.
      if(long_format){
        setcolorder(regression_data, order(colnames(regression_data)))
      }

      censoring_node <- target_node_object$censoring_node


      if (is(censoring_node, "tmle3_Node")) {
        #This node should share the same time/ riskset
        observed <- self$get_tmle_node(censoring_node$name, expand = T, include_id = T, include_time = T, force_time_value = force_time_value, compute_risk_set = F)
        censoring_ids <- observed[as.numeric(observed[[censoring_node$variables]]) == 0, c("id", "t"), with = F]
        #Subset to (id, t) key pairs that are not censored.

        if(drop_censored) {
          regression_data <- regression_data[!.(censoring_ids$id, censoring_ids$t) ]
        } else {
          #Impute to 0
          regression_data[.(censoring_ids$id, censoring_ids$t), (outcome) := 0]

        }
      }




      folds <- self$folds
      #convert folds for ids to rows
      folds <- origami::id_folds_to_folds(folds, match(regression_data$id, unique(self$id)))
      if (FALSE & nrow(regression_data) < self$nrow) {
        #regression_data <- regression_data[indices, ]
        data_id_t <- self$data[, c("id", "t"), with = F]
        #This should
        indices <- data_id_t[regression_data[, c("id", "t"), with = F],  which =  T]
        folds <- sl3::subset_folds(folds, indices)
      }


      #regression_data <- Shared_Data$new(regression_data, force_copy = F)

      if(is_time_variant){
        nodes$covariates <- union(nodes$covariates, "t")
      } else {
        nodes$covariates <- setdiff(nodes$covariates, "t")
      }

      if(is.numeric(bin_num) & include_bins){

        #TODO collision variables?
        regression_data$bin_num <- bin_num
        nodes$covariates <- union(nodes$covariates, "bin_num")
      }
      # Not needed if rbindlist uses names
      if(long_format) {
        setcolorder(regression_data, order(colnames(regression_data)))
        nodes$covariates <- sort(nodes$covariates)
      }

      suppressWarnings({
        regression_task <- sl3_Task$new(
          regression_data,
          nodes = nodes,
          outcome_type = target_node_object$variable_type,
          folds = folds
        )
      })
      if(!is.numeric(force_time_value) & cache_task){
        assign(cache_key, regression_task, private$.node_cache)
      }

      return(regression_task)
    },
    generate_counterfactual_task = function(uuid, new_data,  force_at_risk = NULL, through_data =  F , remove_rows = F, long_format = self$long_format) {
      # for current_factor, generate counterfactual values

      if(is.null(force_at_risk)){
        force_at_risk <- private$.force_at_risk
      }
      if(nrow(new_data)==1){
        node <- colnames(new_data)
        node_var <- sapply(
          node,
          function(node_name) {
            self$npsem[[node_name]]$variables
          }
        )
        nrow <- nrow(self$data)
        new_data <- new_data[rep(1,nrow)]
        setnames(new_data, node, node_var)
        new_task <- self$clone()
        new_column_names <- new_task$add_columns(new_data, uuid)
        new_task$initialize(
          self$internal_data, self$npsem,
          nodes = self$nodes,
          column_names = new_column_names,
          folds = self$folds,
          row_index = self$row_index,
          force_at_risk = force_at_risk,
          summary_measure_columns = private$.summary_measure_columns,
          long_format = long_format
        )
        return(new_task)

      }



      if(!("t" %in% colnames(new_data)) | !("id" %in% colnames(new_data))){

        if(nrow(new_data) == self$nrow){

          node_names <- setdiff(colnames(new_data), c("id", "t"))
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
            row_index = self$row_index,
            force_at_risk = force_at_risk,
            summary_measure_columns = private$.summary_measure_columns,
            long_format = self$long_format
          )
          return(new_task)
        } else {
          through_data = T
        }
      }

      if(!through_data){

        if(!("t" %in% colnames(new_data)) | !("id" %in% colnames(new_data))){
          stop("t and id column not found")
        }


        data <- data.table::copy(self$get_data(self$row_index,))
        node <-  setdiff(colnames(new_data), c("id", "t"))
        vars <- unlist(lapply(node, function(node){
          self$npsem[[node]]$variables
        }))

        if(remove_rows){
          id_t_ex <- fsetdiff(data[t %in% unique(new_data$t), c("id", "t"), with = F], new_data[, c("id", "t"), with = F])
          data <- data[!.(id_t_ex$id, id_t_ex$t), node, with = F]
        } else {
          id_t_ex <- fsetdiff(data[t %in% unique(new_data$t), c("id", "t"), with = F], new_data[, c("id", "t"), with = F])

          data <- data[.(id_t_ex$id, id_t_ex$t), (vars) := NA]
        }


        has_row <- which(unlist(data[.(new_data$id, new_data$t), !is.na(.SD), .SDcols = vars], use.names = F))

        if(length(has_row) == 0){
          append_row_data <- new_data
          alter_row_data <- data.table(NULL)
        } else {
          append_row_data <- new_data[-has_row]
          alter_row_data <- new_data[has_row]
          data[.(alter_row_data$id, alter_row_data$t), (vars) :=  unlist(alter_row_data[, node, with = F])]
        }



        if(nrow(append_row_data) > 0){
          setnames(append_row_data, node, vars)
          data <- rbind(data, append_row_data, fill = T)
          setkey(data, id, t)
        }
        if(any(is.na(data))) {
          numeric_cols <- unlist(data[, lapply(.SD, is.numeric)], use.names = F)
          setnafill(data, "locf", cols = which(numeric_cols))
        }

        new_task <- self$clone()

        #TODO regenerate folds?? But preserve id division? We are adding rows of time

        new_task$initialize(
          data, self$npsem,
          folds = self$folds,
          #row_index = self$row_index,
          t = "t",
          id = "id",
          nodes = self$nodes,
          force_at_risk = force_at_risk,
          summary_measure_columns = private$.summary_measure_columns,
          long_format = self$long_format

        )
        return(new_task)
      }



      # for current_factor, generate counterfactual values
      node_names <- names(new_data)

      node_variables <- sapply(
        node_names,
        function(node_name) {
          self$npsem[[node_name]]$variables
        }
      )

      node_times <- sapply(
        node_names,
        function(node_name) {
          time <- self$npsem[[node_name]]$time
          }
      )

      node_index <- lapply(
        node_times,
        function(time) {
          if(is.null(time)) return(1:nrow(new_data))
          sort(which(self$data$t %in% time))
        }
      )

      old_data <- data.table::copy(self$data[, unique(node_variables), with = F])

      lapply(seq_along(node_index), function(i){
        index <- node_index[[i]]
        var <- node_variables[[i]]

        set(old_data, index, var, new_data[,node_names[[i]],with=F])
      })

      new_data <- old_data

      #setnames(new_data, node_names, node_variables)

      new_task <- self$clone()

      new_column_names <- new_task$add_columns(new_data, uuid)


      new_task$initialize(
        self$internal_data, self$npsem,
        column_names = new_column_names,
        folds = self$folds,
        row_index = self$row_index,
        force_at_risk = force_at_risk,
        summary_measure_columns = private$.summary_measure_columns,
        long_format = self$long_format
      )

      return(new_task)

    },
    next_in_chain = function(...) {
      return(super$next_in_chain(npsem = self$npsem, force_at_risk = self$force_at_risk,
                                 summary_measure_columns = self$summary_measure_columns,
                                 long_format = self$long_format, ...))
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
        #Only worry about removing unque ids, not rows
        new_ids <- self$id[row_index]
        new_fold_index <- match(unique(new_ids), unique(self$id))
        new_folds <- sl3::subset_folds(self$folds, new_fold_index)
      }

      new_task$initialize(
        self$internal_data, self$npsem,
        column_names = self$column_names,
        folds = new_folds,
        row_index = row_index,
        id = "id",
        time = "t",
        force_at_risk = private$.force_at_risk,
        summary_measure_columns = private$.summary_measure_columns,
        long_format = self$long_format
      )
      return(new_task)
    }
  ),
  active = list(
    folds = function(new_folds) {
      if (!missing(new_folds)) {
        private$.folds <- new_folds
      } else if (is.numeric(private$.folds)) {

        # if an integer, create new_folds object but pass integer to V argument
        if(TRUE) {
          new_folds <- origami::make_folds(length(unique(self$id)), V = private$.folds)
        }
        else if (self$has_node("id")) {
          new_folds <- origami::make_folds(
            cluster_ids = self$id,
            V = private$.folds
          )
        } else {
          new_folds <- origami::make_folds(n = self$nrow, V = private$.folds)
        }
        private$.folds <- new_folds
      } else if (is.null(private$.folds)) {
        # generate folds now if never specified
        if(TRUE) {
          new_folds <- origami::make_folds(length(unique(self$id)))
        }
        else if (self$has_node("id")) {
          new_folds <- origami::make_folds(cluster_ids = self$id)
        } else {
          new_folds <- origami::make_folds(n = self$nrow)
        }
        private$.folds <- new_folds
      }
      return(private$.folds)
    },
    long_format = function(){
      private$.long_format
    },
    npsem = function() {
      return(private$.npsem)
    },
    data = function() {
      all_variables <- unlist(lapply(self$npsem, `[[`, "variables"))
      # I need self$data to give me t and id, so lets include nodes
      all_variables <- union(all_variables, c(unlist(self$nodes), private$.summary_measure_columns))
      self$get_data(columns = all_variables)
    },
    summary_measure_columns = function(){
      private$.summary_measure_columns
    },
    force_at_risk = function(force_at_risk = NULL){
      # By changing this value (after cloning original task)
      # The user can quickly construct a counterfactual version of task where everyone is at risk always
      # Useful when you want entire hazard function in one go (ignoring degeneracy).
      if(!is.null(force_at_risk)){
        private$.force_at_risk <- force_at_risk
        if(private$.force_at_risk){
          add <- "T"
        } else{
          add <- ""
        }

        private$.uuid <- paste0(add, digest(self$data))

      }
      private$.force_at_risk
    }
  ),
  private = list(
    .npsem = NULL,
    .node_cache = NULL,
    .force_at_risk = F,
    .summary_measure_columns = NULL,
    .long_format = NULL
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
