
#' @export
threshold_likelihood <- function(tmle_task, learner_list, cutoffs, bins = 10) {
  # covariates
  W_factor <- define_lf(LF_emp, "W")

  # treatment (bound likelihood away from 0 (and 1 if binary))
  A_type <- tmle_task$npsem[["A"]]$variable_type
  if (A_type$type == "continous") {
    A_bound <- c(1 / tmle_task$nrow, Inf)
  } else if (A_type$type %in% c("binomial", "categorical")) {
    A_bound <- 0.025
  } else {
    A_bound <- NULL
  }


  A_factor <- define_lf(LF_fit, "A", learner = Lrnr_CDF$new(learner_list[["A"]], bins, cutoffs), type = "mean", bound = A_bound)

  # outcome
  Y_factor <- LF_fit$new("Y", Lrnr_thresh$new(learner_list[["Y"]], tmle_task$npsem[["A"]]$variables, cutoffs =cutoffs ), type = "mean")


  # construct and train likelihood
  factor_list <- list(W_factor, A_factor, Y_factor)

  # add outcome censoring factor if necessary
  if (!is.null(tmle_task$npsem[["Y"]]$censoring_node)) {
    if (is.null(learner_list[["delta_Y"]])) {
      stop("Y is subject to censoring, but no learner was specified for censoring mechanism delta_Y")
    }

    delta_Y_factor <- define_lf(LF_fit, "delta_Y", learner = learner_list[["delta_Y"]], type = "mean", bound = c(0.025, 1))
    factor_list <- c(factor_list, delta_Y_factor)
  }

  if (!is.null(tmle_task$npsem[["A"]]$censoring_node)) {
    stop("A is subject to censoring, this isn't supported yet")
  }

  if (!is.null(tmle_task$npsem[["W"]]$censoring_node)) {
    stop("W is subject to censoring, this isn't supported yet")
  }

  likelihood_def <- Likelihood$new(factor_list)
  likelihood <- likelihood_def$train(tmle_task)
  return(likelihood)
}



#Construct lrnr that chains the S column
#' @export
Lrnr_thresh <- R6::R6Class(
  classname = "Lrnr_thresh", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(lrnr = make_learner(Lrnr_glm), strata_variable, cutoffs,
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
#' @export
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


      data <- task$data
      Y <- task$Y
      cutoffs <- as.vector(quantile(Y, seq(0, 1, length.out = num_bins)))

      Y <- as.factor(cutoffs[findInterval(Y, cutoffs, left.open = self$params$type != "left-continuous", all.inside
 = T)])

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

      data <- task$data
      Y <- task$Y

      Y <- as.factor(cutoffs[findInterval(Y, cutoffs,
                                          left.open = self$params$type != "left-continuous",
                                          all.inside = T)])

      column_names <- task$add_columns(data.table(Y = Y))

      task <- task$clone()
      nodes <- task$nodes
      nodes$outcome <- "Y"
      task$initialize(
        task$internal_data,
        nodes = nodes,
        folds = task$folds,
        column_names = column_names,
        row_index = task$row_index,
        outcome_type = "categorical",
        outcome_levels = as.factor(cutoffs[-1]))

      #task <- task$next_in_chain(outcome = "Y", column_names = column_names, outcome_type = "categorical", outcome_levels = as.factor(cutoffs))

      pooled_task <- pooled_hazard_task(task, trim = F)

      predictions <- matrix(lrnr$predict(pooled_task), nrow = task$nrow)

      predictions <- cbind(rep(0, task$nrow), 1 - t(apply(1-predictions, 1, cumprod)))

      #Interpolate to get values at desired cutoffs
      predictions <- rbindlist(lapply(1:nrow(predictions), function(i){
        approx(orig_cutoffs,predictions[i,], xout = threshs, rule = 2, yright = 1, yleft = 0 )
      }))
      return(as.vector(predictions[,2][[1]]))
    }
  )
)



