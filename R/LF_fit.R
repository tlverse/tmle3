#' Fitting Likelihoods
#'
#' @importFrom R6 R6Class
#'
#' @export
#
LF_fit <- R6Class(
  classname = "LF_fit",
  portable = TRUE,
  class = TRUE,
  inherit = LF_base,
  public = list(
    initialize = function(name, learner, ...) {
      private$.name <- name
      private$.learner <- learner
    },

    get_learner_task = function(task) {
        task$get_regression_task(self$name)
    },
    train = function(tmle_task) {
      super$train(tmle_task)
      learner <- self$learner
      if (!learner$is_trained) {
        learner_task <- self$get_learner_task(tmle_task)
        private$.learner <- learner$train(learner_task)
      }
    },
    get_prediction = function(tmle_task) {
      learner_task <- self$get_learner_task(tmle_task)
      learner <- self$learner
      preds <- learner$predict(learner_task)
      return(preds)
    },
    get_likelihood = function(tmle_task, only_observed = FALSE) {
      learner_task <- self$get_learner_task(tmle_task)
      preds <- self$get_prediction(tmle_task)
      #todo: add support for multinomial
      outcome_type <- self$learner$training_task$outcome_type
      if (only_observed) {
        observed <- outcome_type$format(learner_task$Y)
        if (outcome_type$type == "binomial") {
          likelihood <- ifelse(observed == 1, preds, 1 - preds)
        } else if(outcome_type$type =="categorical"){
          unpacked <- sl3::unpack_predictions(preds)  
          index_mat <- cbind(seq_along(observed), observed)
          likelihood <- unpacked[index_mat]
        } else {
          stop("currently, only binomial and multinomial likelihoods are supported")
        }
      } else {
        if (outcome_type$type == "binomial") {
          likelihood <- cbind(1 - preds, preds)
        } else if(outcome_type$type =="categorical"){
          likelihood <- sl3::unpack_predictions(preds)  
        } else {
          stop("currently, only binomial and multinomial likelihoods are supported")
        }
      }
      return(likelihood)
    }
  ),
  active = list(
    learner = function() {
      return(private$.learner)
    }
  ),
  private = list(
    .name = NULL,
    .learner = NULL
  )
)
