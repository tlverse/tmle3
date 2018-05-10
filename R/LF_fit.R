#' Likelihood Factor Estimated from Data using sl3.
#'
#' Uses an \code{sl3} learner to estimate a likelihood factor from data.
#' Inherits from \code{\link{LF_base}}; see that page for documentation on likelihood factors in general.
#' Currently, predictions are memoized to speed up multiple calls to \code{get_likelihood} and \code{get_mean} for the same data.
#'
#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @family Likelihood objects
#' @keywords data
#'
#' @return \code{LF_base} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @section Constructor:
#'   \code{define_lf(LF_fit, name, learner, ..., type = "density")}
#'
#'   \describe{
#'     \item{\code{name}}{character, the name of the factor. Should match a node name in the nodes specified by \code{\link{tmle3_Task}$npsem}
#'     }
#'     \item{\code{learner}}{An sl3 learner to be used to estimate the factor
#'     }
#'     \item{\code{...}}{Not currently used.
#'     }
#'     \item{\code{type}}{character, either "density", for conditional density or, "mean" for conditional mean
#'     }
#'     }
#'
#' @section Fields:
#' \describe{
#'     \item{\code{learner}}{The learner or learner fit object}
#'     }
#'
#' @export
LF_fit <- R6Class(
  classname = "LF_fit",
  portable = TRUE,
  class = TRUE,
  inherit = LF_base,
  public = list(
    initialize = function(name, learner, ..., type="density") {
      super$initialize(name, ..., type = type)
      private$.learner <- learner
    },
    delayed_train = function(tmle_task) {
      # just return prefit learner if that's what we have
      # otherwise, make a delayed fit and return that
      if (self$learner$is_trained) {
        return(self$learner)
      }

      outcome_node <- self$name
      # todo: we don't support delayed tmle_tasks
      # this relates to the issue of cross_validation for delayed tasks
      # and not being able to generate delayed calls on-the-fly
      learner_task <- tmle_task$get_regression_task(outcome_node)
      learner_fit <- delayed_learner_train(self$learner, learner_task)
      return(learner_fit)
    },
    train = function(tmle_task, learner_fit) {
      super$train(tmle_task)
      private$.learner <- learner_fit
    },
    get_mean = function(tmle_task) {
      if (self$memoize_predictions) {
        uuid <- tmle_task$uuid
        preds <- private$.memoized[[uuid]]
        if (is.null(preds)) {
          learner_task <- tmle_task$get_regression_task(self$name)
          learner <- self$learner
          preds <- learner$predict(learner_task)
          private$.memoized[[uuid]] <- preds
        }
      } else {
        learner_task <- tmle_task$get_regression_task(self$name)
        learner <- self$learner
        preds <- learner$predict(learner_task)
      }
      return(preds)
    },
    get_likelihood = function(tmle_task) {
      learner_task <- tmle_task$get_regression_task(self$name)
      preds <- self$learner$predict(learner_task)
      outcome_type <- self$learner$training_task$outcome_type
      observed <- outcome_type$format(learner_task$Y)
      if (outcome_type$type == "binomial") {
        likelihood <- ifelse(observed == 1, preds, 1 - preds)
      } else if (outcome_type$type == "categorical") {
        unpacked <- sl3::unpack_predictions(preds)
        index_mat <- cbind(seq_along(observed), observed)
        likelihood <- unpacked[index_mat]
      } else if (outcome_type$type == "continuous") {
        likelihood <- unlist(preds)
      } else {
        stop(sprintf("unsupported outcome_type: %s", outcome_type$type))
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
