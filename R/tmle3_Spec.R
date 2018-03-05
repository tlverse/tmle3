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

      # bound Y if continuous
      Y_node <- node_list$Y
      Y_vals <- unlist(data[, Y_node, with = FALSE])
      Y_variable_type <- variable_type(x = Y_vals)
      if (Y_variable_type$type == "continuous") {
        min_Y <- min(Y_vals)
        max_Y <- max(Y_vals)
        range <- max_Y - min_Y
        lower <- min_Y - 0.1 * range
        upper <- max_Y + 0.1 * range
        Y_variable_type <- variable_type(type = "continuous", bounds = c(lower, upper))
      }

      # make tmle_task
      tmle_nodes <- list(
        define_node("W", node_list$W),
        define_node("A", node_list$A, c("W")),
        define_node("Y", node_list$Y, c("A", "W"), Y_variable_type)
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
        define_lf(LF_fit, "A", type = "density", learner = combined_learner_list[["A"]]),
        define_lf(LF_fit, "Y", type = "mean", learner = combined_learner_list[["Y"]])
      )

      likelihood_def <- Likelihood$new(factor_list)

      # fit_likelihood
      likelihood <- likelihood_def$train(tmle_task)
      return(likelihood)
    },
    make_params = function(tmle_task, likelihood) {
      # todo: generalize
      browser()
      # todo: export and use sl3:::get_levels
      A_levels <- levels(tmle_task$get_tmle_node("A"))
      A_levels <- factor(A_levels, A_levels)
      A_level <- A_levels[[1]]
      tmle_params <- lapply(A_levels, function(A_level) {
        intervention <- define_lf(LF_static, "A", value = A_level)
        tsm <- Param_TSM$new(likelihood, "Y", intervention)
        return(tsm)
      })

      return(tmle_params)
    },
    make_updater = function(likelihood, tmle_params) {
      # todo: generalize
      updater <- tmle3_Update$new(tmle_params)
      likelihood$update_list <- updater
      return(updater)
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
tmle_tsm_all <- function() {
  # todo: unclear why this has to be in a factory function
  tmle3_Spec$new(
    dag = list(Y = c("A", "W"), A = c("W"), W = c()),
    default_learner_list = list(A = make_learner(Lrnr_mean), Y = make_learner(Lrnr_mean))
  )
}
