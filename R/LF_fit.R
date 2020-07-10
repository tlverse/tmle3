#' Likelihood Factor Estimated from Data using sl3.
#'
#' Uses an \code{sl3} learner to estimate a likelihood factor from data.
#' Inherits from \code{\link{LF_base}}; see that page for documentation on likelihood factors in general.
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
    initialize = function(name, learner, is_time_variant = FALSE, ..., type = "density") {
      super$initialize(name, ..., type = type)
      private$.learner <- learner
      # TODO: add parameter is_time_variant
      private$.is_time_variant <- is_time_variant
    },
    delayed_train = function(tmle_task) {
      # just return prefit learner if that's what we have
      # otherwise, make a delayed fit and return that
      if (self$learner$is_trained) {
        return(self$learner)
      }

      outcome_node <- self$name

      # fit scaled task for bounded continuous
      learner_task <- tmle_task$get_regression_task(outcome_node, scale = TRUE, drop_censored=TRUE, 
        is_time_variant = self$is_time_variant)
      learner_fit <- delayed_learner_train(self$learner, learner_task)
      return(learner_fit)
    },
    train = function(tmle_task, learner_fit) {
      super$train(tmle_task)
      private$.learner <- learner_fit
    },
    get_mean = function(tmle_task, fold_number) {
      # TODO: prediction is made on all data, so is_time_variant is set to TRUE
      learner_task <- tmle_task$get_regression_task(self$name, is_time_variant = TRUE)
      learner <- self$learner
      preds <- learner$predict_fold(learner_task, fold_number)

      # unscale preds (to handle bounded continuous)
      preds_unscaled <- tmle_task$unscale(preds, self$name)
      return(preds_unscaled)
    },
    get_density = function(tmle_task, fold_number) {
      # TODO: prediction is made on all data, so is_time_variant is set to TRUE
      learner_task <- tmle_task$get_regression_task(self$name, is_time_variant = TRUE)
      learner <- self$learner
      preds <- learner$predict_fold(learner_task, fold_number)

      outcome_type <- self$learner$training_task$outcome_type
      observed <- outcome_type$format(learner_task$Y)
      if (outcome_type$type == "binomial") {
        likelihood <- ifelse(observed == 1, preds, 1 - preds)
      } else if (outcome_type$type == "categorical") {
        unpacked <- sl3::unpack_predictions(as.vector(preds))
        index_mat <- cbind(seq_along(observed), observed)
        likelihood <- unpacked[index_mat]
      } else if (outcome_type$type == "continuous") {
        likelihood <- unlist(preds)
      } else {
        stop(sprintf("unsupported outcome_type: %s", outcome_type$type))
      }
      return(likelihood)
    },
    sample = function(tmle_task, n_samples = NULL, fold_number = "full") {
      # TODO: fold
      if (is.null(tmle_task)) {
        tmle_task <- self$training_task
      }
      if (is.null(n_samples)) {
        return(tmle_task)
      }

      learner_task <- tmle_task$get_regression_task(self$name)
      learner <- self$learner

      outcome_type <- learner$training_task$outcome_type

      if (outcome_type$type == "binomial") {
        # TODO: option to return task
        # TODO: think about how folds should be structured on resample
        # need to keep ids the same
        # probably also predict using training set fits
        preds <- learner$predict_fold(learner_task, "full")
        values <- sapply(preds, function(p) rbinom(n_samples, 1, p))
      } else if (outcome_type$type == "categorical") {
        preds <- learner$predict_fold(learner_task, "full")
        unpacked <- sl3::unpack_predictions(as.vector(preds))
        values <- apply(
          unpacked, 1,
          function(probs) {
            apply(
              rmultinom(n_samples, 1, probs) == 1, 2,
              function(onehots) outcome_type$levels[which(onehots)]
            )
          }
        )
      } else if (outcome_type$type == "continuous") {
        if ("sampling" %in% learner$properties) {
          values <- learner$sample(learner_task, n_samples, "full")
        } else {
          values <- matrix(nrow = n_samples, ncol = tmle_task$nrow)
          for (i in 1:tmle_task$nrow) {
            subject <- tmle_task[i]
            f_X <- function(a) {
              cf_data <- data.table(a)
              setnames(cf_data, names(cf_data), self$name)
              subject_a <- subject$generate_counterfactual_task(UUIDgenerate(), cf_data)

              pred <- learner$predict_fold(subject_a$get_regression_task(self$name), "full")
              likelihood <- unlist(pred)

              return(likelihood)
            }
            samples <- AR.Sim(n_samples, f_X,
              xlim = c(min(learner$training_task$Y), max(learner$training_task$Y))
            )
            values[, i] <- samples
          }
        }
      } else {
        stop(sprintf("unsupported outcome_type: %s", outcome_type$type))
      }
      values <- t(values)
      return(values)
    }
  ),
  active = list(
    learner = function() {
      return(private$.learner)
    },
    is_time_variant = function() {
      return(private$.is_time_variant)
    }
  ),
  private = list(
    .name = NULL,
    .learner = NULL,
    .is_time_variant = NULL
  )
)
