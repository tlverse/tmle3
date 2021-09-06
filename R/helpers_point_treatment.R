#' Helper Functions for Point Treatment
#'
#' Handles the common W (covariates), A (treatment/intervention), Y (outcome) data structure
#'
#' @param data a \code{data.frame}, or \code{data.table} containing data for use in estimation
#' @param node_list a list of character vectors, listing the variables that comprise each node
#' @param variable_types a list of variable types, one for each node. If missing, variable types will be guessed
#' @param tmle_task a \code{\link{tmle3_Task}} as constructed via \code{point_tx_task}
#' @param learner_list a list of sl3 learners, one for A and one for Y to be used for likelihood estimation
#' @param ... extra arguments.
#' @export
#' @rdname point_tx
point_tx_npsem <- function(node_list, variable_types = NULL, scale_outcome = TRUE, include_variance_node = FALSE) {
  # make tmle_task
  npsem <- list(
    define_node("W", node_list$W, variable_type = variable_types$W),
    define_node("A", node_list$A, c("W"), variable_type = variable_types$A),
    define_node("Y", node_list$Y, c("A", "W"), variable_type = variable_types$Y, scale = scale_outcome)
  )
  if (include_variance_node) {
    npsem$var_Y <- define_node("var_Y", node_list$Y, c("A", "W"), variable_type = variable_types$var_Y, scale = FALSE)
  }

  return(npsem)
}

#' @export
#' @rdname point_tx
point_tx_task <- function(data, node_list, variable_types = NULL, scale_outcome = TRUE, include_variance_node = FALSE, ...) {
  setDT(data)

  npsem <- point_tx_npsem(node_list, variable_types, scale_outcome, include_variance_node)

  if (!is.null(node_list$id)) {
    tmle_task <- tmle3_Task$new(data, npsem = npsem, id = node_list$id, ...)
  } else {
    tmle_task <- tmle3_Task$new(data, npsem = npsem, ...)
  }

  return(tmle_task)
}

#' @export
#' @rdname point_tx
point_tx_likelihood <- function(tmle_task, learner_list) {
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

  A_factor <- define_lf(LF_fit, "A", learner = learner_list[["A"]], bound = A_bound)

  # outcome
  Y_factor <- define_lf(LF_fit, "Y", learner = learner_list[["Y"]], type = "mean")


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

  # If conditional variance needs to be estimated, do so.
  if ("var_Y" %in% names(tmle_task$npsem)) {
    if (is.null(learner_list[["var_Y"]])) {
      learner_list[["var_Y"]] <- Lrnr_glmnet$new(family = "poisson")
      warning("Node var_Y is in npsem but no learner is provided in `learner_list`. Defaulting to glmnet with `poisson` family.")
    }
    if (tmle_task$npsem[["Y"]]$variable_type$type == "binomial") {
      mean_fun <- function(task, likelihood, tmle_task) {
        EY <-  likelihood$get_likelihood(tmle_task, "Y")
        return(EY * (1 - EY))
      }
      LF_var_Y <- LF_known$new("var_Y", mean_fun = mean_fun, base_likelihood = likelihood, type = "mean")
    } else {
      task_generator <- function(tmle_task, base_likelihood) {
        EY <- base_likelihood$get_likelihood(tmle_task, "Y")
        Y <- tmle_task$get_tmle_node("Y", format = TRUE)[[1]]
        outcome <- (Y - EY)^2
        task <- tmle_task$get_regression_task("Y")
        column_names <- task$add_columns(data.table("var_Y" = outcome))
        task <- task$next_in_chain(outcome = "var_Y", column_names = column_names)
      }
      LF_var_Y <- LF_derived$new("var_Y", learner_list[["var_Y"]], likelihood, task_generator = task_generator, type = "mean")
    }
    likelihood$add_factors(LF_var_Y)
  }

  return(likelihood)
}
