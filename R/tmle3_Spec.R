#' Defines a tmle (minus the data)
#'
#' Current limitations:
#' pretty much tailored to Param_TSM
#' see todos for places generalization can be added
#' @importFrom R6 R6Class
#'
#' @export
#
tmle3_Spec <- R6Class(
  classname = "tmle3_Spec",
  portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(dag, default_learner_list = NULL) {
      private$.dag <- dag
      private$.default_learner_list <- default_learner_list
    },
    make_tmle_task = function(data, node_list) {
      # todo: generalize
      tmle_nodes <- list(
        define_node("W", node_list$W),
        define_node("A", node_list$A, c("W")),
        define_node("Y", node_list$Y, c("A", "W"))
      )

      tmle_task <- tmle3_Task$new(data, tmle_nodes = tmle_nodes)
      return(tmle_task)
    },
    make_likelihood = function(tmle_task, learner_list = NULL) {
      combined_learner_list <- self$default_learner_list
      combined_learner_list[names(learner_list)] <- learner_list[names(learner_list)]

      # todo: generalize
      factor_list <- list(
        define_lf(LF_np, "W", NA),
        define_lf(LF_fit, "A", combined_learner_list[["A"]]),
        define_lf(LF_fit, "Y", combined_learner_list[["Y"]], type = "mean")
      )

      likelihood_def <- Likelihood$new(factor_list)

      # fit_likelihood
      likelihood <- likelihood_def$train(tmle_task)
      return(likelihood)
    },
    make_params = function(tmle_task, likelihood) {
      # todo: generalize

      # todo: export and use sl3:::get_levels
      A_levels <- levels(tmle_task$get_tmle_node("A"))
      A_levels <- factor(A_levels, A_levels)
      A_level <- A_levels[[1]]
      tmle_params <- lapply(A_levels, function(A_level) {
        intervention <- define_cf(define_lf(LF_static, "A", value = A_level))
        tsm <- Param_TSM$new(intervention)
        # todo: maybe put setup in constructor
        tsm$setup(likelihood)
        return(tsm)
      })

      return(tmle_params)
    },
    make_updater = function(likelihood, tmle_params) {
      # todo: generalize
      updater <- tmle3_Update$new(tmle_params)
      likelihood$update_list <- updater
      return(updater)
    },
    tmle3 = function(data, node_list, learner_list = NULL) {
      tmle_task <- self$make_tmle_task(data, node_list)
      likelihood <- self$make_likelihood(tmle_task, learner_list)
      tmle_params <- self$make_params(tmle_task, likelihood)
      updater <- self$make_updater(likelihood, tmle_params)
      fit <- fit_tmle3(tmle_task, likelihood, tmle_params, updater)
      return(fit)
    }
  ),
  active = list(
    dag = function() {
      return(private$.dag)
    },
    default_learner_list = function() {
      return(private$.default_learner_list)
    }
  ),
  private = list(
    .dag = NULL,
    .default_learner_list = NULL
  )
)

#' All Treatment Specific Means
#'
#' O=(W,A,Y)
#' W=Covariates
#' A=Treatment (binary or categorical)
#' Y=Outcome (binary or bounded continuous)
#' @importFrom sl3 make_learner Lrnr_mean
#' @export
tmle_tsm_all <- function(){
  #todo: unclear why this has to be in a factory function
  tmle3_Spec$new(dag = list(Y = c("A", "W"), A = c("W"), W = c()),
                 default_learner_list = list(A = make_learner(Lrnr_mean), Y = make_learner(Lrnr_mean)))
}
