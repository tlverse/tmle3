
Lrnr_thresh <- R6::R6Class(
  classname = "Lrnr_thresh", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(lrnr, strata_variable, cutoffs,
                          ...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("continuous", "binomial"),

    .train = function(task) {
      args <- self$params
      strata_variable <- args$strata_variable
      lrnr <- args$lrnr
      outcome_type <- self$get_outcome_type(task)
      data <- task$data
      cutoffs <- args$cutoffs
      data_list <- list()

      for(cutoff in cutoffs) {
        Xcopy <- copy(data)
        Xcopy$bin <- cutoff
        Xcopy$Ind <- as.numeric(Xcopy[[strata_variable]] >= cutoff)
        Xcopy[[strata_variable]] <- NULL
        data_list[[as.character(cutoff)]] <- Xcopy
      }
      data <- rbindlist(data_list)

      nodes <- task$nodes
      nodes$covariates <- union(setdiff(task$nodes$covariates, strata_variable), c("Ind", "bin"))
      task <- sl3_Task$new(data, nodes = nodes)
      lrnr <- lrnr$train(task)

      return(list(lrnr = lrnr, task = task))
    },
    .predict = function(task = NULL) {

      args <- self$params
      cutoffs <- args$cutoffs

      strata_variable <- args$strata_variable

      outcome_type <- self$get_outcome_type(task)
      data <- task$data

      data_list <- list()

      for(cutoff in cutoffs) {
        Xcopy <- copy(data)
        Xcopy$bin <- cutoff
        Xcopy$Ind <- as.numeric(Xcopy[[strata_variable]] >= cutoff)
        Xcopy[[strata_variable]] <- NULL
        data_list[[as.character(cutoff)]] <- Xcopy
      }
      data <- rbindlist(data_list)
      nodes <- task$nodes
      nodes$covariates <- union(setdiff(task$nodes$covariates, strata_variable), c("Ind", "bin"))
      task <- sl3_Task$new(data, nodes = nodes)
      predictions <- self$fit_object$lrnr$predict(task)
      return(as.vector(predictions))
    }
  )
)

Lrnr_CDF <- R6::R6Class(
  classname = "Lrnr_CDF", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(lrnr, num_bins, threshs, type = "left-continuous",
                          ...) {

      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("continuous", "binomial"),

    .train = function(task) {
      args <- self$params

      lrnr <- args$lrnr
      num_bins <- args$num_bins

      outcome_type <- self$get_outcome_type(task)
      data <- task$data
      Y <- task$Y
      cutoffs <- as.vector(quantile(Y, seq(0, 1, length.out = num_bins)))
      print(length(cutoffs))
      Y <- as.factor(cutoffs[findInterval(Y, cutoffs, left.open = self$params$type != "left-continuous", all.inside
 = T)])
      print(length(unique(Y)))
      column_names <- task$add_columns(data.table(Y = Y))
      task <- task$next_in_chain(outcome = "Y", column_names = column_names)

      pooled_task <- pooled_hazard_task(task)
      lrnr <- lrnr$train(pooled_task)
      return(list(lrnr = lrnr, cutoffs = cutoffs))
    },
    .predict = function(task = NULL) {
      args <- self$params
      cutoffs <- self$fit_object$cutoffs
      orig_cutoffs <- cutoffs
      #cutoffs <- c(-Inf,cutoffs)
      threshs <- args$threshs
      lrnr <- self$fit_object$lrnr
      outcome_type <- self$get_outcome_type(task)
      data <- task$data
      Y <- task$Y
      Y <- as.factor(cutoffs[findInterval(Y, cutoffs,
                                          left.open = self$params$type != "left-continuous",
                                          all.inside = T)])

      column_names <- task$add_columns(data.table(Y = Y))
      task <- task$next_in_chain(outcome = "Y", column_names = column_names)
      pooled_task <- pooled_hazard_task(task, trim = F)

      predictions <- matrix(lrnr$predict(pooled_task), nrow = task$nrow)
      predictions <- cbind(rep(0, task$nrow), 1 - t(apply(1-predictions, 1, cumprod)))
      print(dim(predictions))
      print(length(cutoffs))
      #Interpolate to get values at desired cutoffs
      predictions <- rbindlist(lapply(1:nrow(predictions), function(i){
        approx(orig_cutoffs,predictions[i,], xout = threshs, rule = 2, yright = 1, yleft = 0 )
      }))
      return(as.vector(predictions[,2][[1]]))
    }
  )
)



