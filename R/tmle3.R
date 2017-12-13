# default_lrnr_submodel <- make_learner(Lrnr_glm_fast, intercept = FALSE, transform_offset = TRUE)

#' Fit Targeted Maximum Likelihood Estimator
#'
#' @param likelihood ...
#' @param task ...
#' @param param ...
#' @param lrnr_submodel ...
#'
#' @export
#' @import sl3
#' @import data.table
fit_tmle_likelihood <- function(likelihood, task, param, lrnr_submodel) {

  # prepare learner and task for tmle update
  lrnr_tmle_update <- make_learner(Lrnr_tmle_update, param, lrnr_submodel)
  fit_for_submodel <- sl3::customize_chain(likelihood, lrnr_tmle_update$prepare_task)
  submodel_task <- fit_for_submodel$chain()

  # fit update
  update_fit <- lrnr_tmle_update$train(submodel_task)

  # updated fit is a pipeline with the intial fit and then the updated fit
  tmle_pipe <- sl3::Pipeline$new(fit_for_submodel, update_fit)

  # modify likelihood with fluctuation
  lf_y_updated <- LF_fit$new("Y", tmle_pipe, expects_tmle_task = TRUE)
  tmle_likelihood <- likelihood$modify_factors(list(lf_y_updated))

  return(tmle_likelihood)
}
