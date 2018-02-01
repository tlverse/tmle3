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
    initialize = function(name, learner, type="density", ...) {
      private$.name <- name
      private$.learner <- learner
      private$.type <- type
    },
    get_learner_task = function(task) {
      task$get_regression_task(self$name)
    },
    delayed_train = function(tmle_task){
      # just return prefit learner if that's what we have
      # otherwise, make a delayed fit and return that
      if (self$learner$is_trained) {
        return(self$learner)
      }
      
      outcome_node <- self$name
      # todo: we don't support delayed tmle_tasks
      # this relates to the issue of cross_validation for delayed tasks
      # and not being able to generate delayed calls on-the-fly
      get_learner_task <- function(tmle_task, outcome_node){
        tmle_task$get_regression_task(outcome_node)
      }
      
      # delayed_get_learner_task <- delayed_fun(get_learner_task)
      learner_task <- get_learner_task(tmle_task, outcome_node)
      learner_fit <- delayed_learner_train(self$learner, learner_task)
      return(learner_fit)
    },
    train = function(tmle_task, learner_fit) {
      super$train(tmle_task)
      private$.learner <- learner_fit
    },
    get_prediction = function(tmle_task) {
      if (self$memoize_predictions) {
        uuid <- tmle_task$uuid
        preds <- private$.memoized[[uuid]]
        if (is.null(preds)) {
          learner_task <- self$get_learner_task(tmle_task)
          learner <- self$learner
          preds <- learner$predict(learner_task)
          private$.memoized[[uuid]] <- preds
        }
      } else {
        learner_task <- self$get_learner_task(tmle_task)
        learner <- self$learner
        preds <- learner$predict(learner_task)
      }
      return(preds)
    },
    get_likelihood = function(tmle_task, only_observed = FALSE) {
      learner_task <- self$get_learner_task(tmle_task)
      preds <- self$get_prediction(tmle_task)
      # todo: add support for multinomial
      outcome_type <- self$learner$training_task$outcome_type
      if (only_observed) {
        observed <- outcome_type$format(learner_task$Y)
        if (outcome_type$type == "binomial") {
          likelihood <- ifelse(observed == 1, preds, 1 - preds)
        } else if (outcome_type$type == "categorical") {
          unpacked <- sl3::unpack_predictions(preds)
          index_mat <- cbind(seq_along(observed), observed)
          likelihood <- unpacked[index_mat]
        } else {
          stop("currently, only binomial and multinomial likelihoods are supported")
        }
      } else {
        if (outcome_type$type == "binomial") {
          likelihood <- cbind(1 - preds, preds)
        } else if (outcome_type$type == "categorical") {
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
    },
    memoize_predictions = function() {
      return(private$.memoize_predictions)
    }
  ),
  private = list(
    .name = NULL,
    .learner = NULL,
    .memoize_predictions = TRUE,
    .memoized = list()
  )
)
