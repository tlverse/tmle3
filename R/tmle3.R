#' Fit Targeted Maximum Likelihood Estimator
#'
#' @param likelihood ...
#' @param task ...
#' @param param ...
#'
#' @export
#
fit_tmle_likelihood <- function(likelihood, task, param) {
  # use GLM to fit a fluctuation regression through submodels
  tmle_fluc <- Lrnr_glm_fast$new(intercept = FALSE, transform_offset = TRUE)
  # tmle_fluc <- make_learner(Lrnr_optim, submodel_logit, loss_loglik_binomial,
  #                           init_0 = TRUE)

  # make pipe that generates fluctuated predictions
  # wrap likelihood in custom chain that defines quantities for fluctuation
  fluc_chain <- function(likelihood, task) {
    Y <- task$get_regression_task(param$outcome_node)$Y
    EY <- likelihood$get_predictions(task, param$outcome_node)
    HA <- param$HA(likelihood, task)

    tmle_data <- data.table(HA = HA, EY = EY, Y = Y)
    fluc_task <- sl3_Task$new(tmle_data, outcome = "Y", offset = "EY",
                              covariates = "HA")
    return(fluc_task)
  }
  fluc_likelihood <- customize_chain(likelihood, fluc_chain)
  fluc_task <- fluc_chain(likelihood, task)
  fluc_fit <- tmle_fluc$train(fluc_task)

  # put in pipe with learner to estimate fluctuation
  tmle_pipe <- Pipeline$new(fluc_likelihood, fluc_fit)

  # update likelihood with fluctuation
  lf_y_updated <- LF_fit$new("Y", tmle_pipe, expects_tmle_task = TRUE)
  tmle_likelihood <- likelihood$modify_factors(list(lf_y_updated))
  return(tmle_likelihood)
}

