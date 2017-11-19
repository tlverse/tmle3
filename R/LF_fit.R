#' Fitting Likelihood Functionals
#'
#' @importFrom R6 R6Class
#'
#' @export
#
LF_fit <- R6Class(classname = "LF_fit",
                             portable = TRUE,
                             class = TRUE,
                             inherit = LF_base,
  public = list(
    initialize = function(name, learner, expects_tmle_task = FALSE, ...) {
      private$.name <- name
      private$.learner <- learner
      private$.expects_tmle_task <- expects_tmle_task
    },

    get_learner_task <- function(task) {
      if (self$expects_tmle_task) {
        return(task)
      } else {
        return(task$get_regression_task(self$name))
      }
    },
    fit_learner <- function(task) {
      learner <- self$learner
      if (!learner$is_trained) {
        learner_task <- self$get_learner_task(task)
        private$.learner <- learner$train(learner_task)
      }
    },
    get_prediction <- function(task) {
      learner_task <- self$get_learner_task(task)
      learner <- self$learner
      preds <- learner$predict(learner_task)
      return(preds)
    },
    get_likelihood <- function(task, only_observed = FALSE) {
      learner_task <- self$get_learner_task(task)
      preds <- self$get_prediction(task)
      # variable type should come from node type
      # always store fit tmle_task (or at least fit node_list) for this info
      outcome_type <- self$learner$training_task$outcome_type
      if (only_observed) {
        observed <- outcome_type$format(learner_task$Y)
        if (outcome_type$type == "binomial") {
          likelihood <- ifelse(observed == 1, preds, 1 - preds)
        } else{
          stop("currently, only binomial likelihoods are supported")
        }
      } else {
        if (outcome_type$type == "binomial") {
          likelihood <- cbind(1-preds, preds)
        } else {
          stop("currently, only binomial likelihoods are supported")
        }
      }
      return(likelihood)
    }
  ),
  active = list(
    learner <- function() {
      return(private$.learner)
    },
    expects_tmle_task <- function() {
      return(private$.expects_tmle_task)
    }
  ),
  private = list(
    .name = NULL,
    .learner = NULL,
    .expects_tmle_task = FALSE
  )
)

