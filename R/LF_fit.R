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
    initialize = function(name, learner, is_time_variant = FALSE, is_level_variant = F, include_degeneracy = T,  ..., type = "density") {
      if(length(name) > 1){
        if(!(is_time_variant | is_level_variant)) {
          warning("You have specified that this a pooled regression but did not specify time or level variance.")
        }
        # Is time variant is needed to allow predictions for a single node (of pooled regression)
        #by using the individual nodes rgeression task (with time).
      }
      super$initialize(name, ..., type = type)
      private$.learner <- learner
      # TODO: add parameter is_time_variant
      private$.is_time_variant <- is_time_variant
      private$.is_level_variant <- is_level_variant
      private$.include_degeneracy <- include_degeneracy
    },
    delayed_train = function(tmle_task) {
      # just return prefit learner if that's what we have
      # otherwise, make a delayed fit and return that
      if (self$learner$is_trained) {
        return(self$learner)
      }

      outcome_node <- self$name

      # fit scaled task for bounded continuous
      learner_task <- tmle_task$get_regression_task(outcome_node,
        scale = TRUE, drop_censored = TRUE, expand = F)
      learner_fit <- delayed_learner_train(self$learner, learner_task)
      return(learner_fit)
    },
    train = function(tmle_task, learner_fit) {
      tmle_nodes <- lapply(self$names, function(node) tmle_task$npsem[[node]])
      private$.variable_type <- lapply(tmle_nodes, function(node) node$variable_type)
      private$.training_task <- tmle_task
      private$.learner <- learner_fit
    },
    get_likelihood = function(tmle_task, fold_number = "full", expand = T, node = NULL) {
      if(is.null(node)) node <- self$name
      if (self$type == "mean") {
        values <- self$get_mean(tmle_task, fold_number, expand = expand, node = node)
      } else {
        values <- self$get_density(tmle_task, fold_number, expand = expand, node = node)
      }
      if (!is.null(self$bound)) {
        values <- bound(values, self$bound)
        #values[, (node) := bound(values[[node]], self$bound)]
      }

      return(values)
    },
    get_mean = function(tmle_task, fold_number, node = NULL,  check_at_risk = T, to_wide = F, drop_id = F, drop_time = F, drop = T, expand = T, quick_pred = F) {
      # TODO: prediction is made on all data, so is_time_variant is set to TRUE
      #
      if(is.null(node)) node <- self$name
      learner_task <- tmle_task$get_regression_task(node, expand = expand, include_bins = self$is_level_variant, bin_num = match(node, self$name), is_time_variant = self$is_time_variant)
      learner <- self$learner
      preds <- learner$predict_fold(learner_task, fold_number)
      if(quick_pred){
        return(data.table(preds))
      }
      # unscale preds (to handle bounded continuous)
      preds <- tmle_task$unscale(preds, node)
      #preds <- as.data.table(preds)

      # data <-  learner_task$get_data()
      #
      # # For conditional means, degenerate value is exactly the value to be predicted
      # if(check_at_risk & "at_risk" %in% colnames(data) & self$include_degeneracy) {
      #   # By setting check_at_risk = F then one can obtain the counterfactual predictions
      #   # conditioned on everyone being at risk.
      #   names_degen_val <- paste0("degeneracy_value_", learner_task$nodes$outcome)
      #   # TODO Support multivariate outcome
      #   assertthat::assert_that(all(names_degen_val %in% colnames(data)), msg = "If at_risk is a column then last_val must be as well.")
      #   not_at_risk <- which(data$at_risk == 0)
      #   if(length(not_at_risk)>0){
      #     degen_val <- data[not_at_risk, names_degen_val, with = F]
      #     set(preds, not_at_risk, names(preds) ,  degen_val)
      #   }
      #
      # }

      # preds$id <- rep(learner_task$id, nrow(preds)/length(learner_task$id))
      # preds$t <- rep(learner_task$time, nrow(preds)/length(learner_task$time))
      # setnames(preds, c(node, "id", "t"))
      # if(to_wide){
      #   preds <- reshape(preds, idvar = "id", timevar = "t", direction = "wide")
      #   if(length(node) + 1  == ncol(preds)){
      #     setnames(preds, c("id", node))
      #   }
      #
      # }
      # if(drop_id & "id" %in% colnames(preds)) preds$id <- NULL
      # if(drop_time & "t" %in% colnames(preds)) preds$t <- NULL
      # if(drop & ncol(preds) == 1) preds <- unlist(preds, use.names = F)

      return(preds)
    },
    get_density = function(tmle_task, fold_number, node = NULL,  check_at_risk = T, to_wide = F, drop_id = F, drop_time = F, drop = T, expand = T, quick_pred = F) {
      # TODO: prediction is made on all data, so is_time_variant is set to TRUE
      if(is.null(node)) node <- self$name


      learner_task <- tmle_task$get_regression_task(node, expand = expand, include_bins = self$is_level_variant, bin_num = match(node, self$name), is_time_variant = self$is_time_variant)

      learner <- self$learner
      preds <- learner$predict_fold(learner_task, fold_number)
      if(quick_pred){
        return(preds)
      }

      outcome_type <- self$learner$training_task$outcome_type
      observed <- outcome_type$format(learner_task$Y)

      data <-  learner_task$get_data()

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
      # if(check_at_risk & "at_risk" %in% colnames(data) & self$include_degeneracy ) {
      #   # By setting check_at_risk = F then one can obtain the counterfactual predictions
      #   # conditioned on everyone being at risk.
      #   names_degen_val <- paste0("degeneracy_value_", learner_task$nodes$outcome)
      #   # TODO Support multivariate outcome
      #   assertthat::assert_that(all(names_degen_val %in% colnames(data)), msg = "If at_risk is a column then last_val must be as well.")
      #   not_at_risk <- which(data$at_risk == 0)
      #   if(length(not_at_risk)>0){
      #     degen_val <- data[not_at_risk, names_degen_val, with = F]
      #     #If not at risk then equal to degenerate value with prob 1
      #     #TODO  we should probably not compute the likelihood for all the not at risk people?
      #     #Since we change their value anyway?
      #     likelihood[not_at_risk] <- as.numeric(observed[not_at_risk] == degen_val)
      #   }
      #
      # }
      # likelihood <- data.table(likelihood)
      # likelihood$id <- rep(learner_task$data$id, length(preds)/length(learner_task$data$id))
      #
      #
      # likelihood$t <- learner_task$data$t
      #
      # setnames(likelihood, c(paste0(node, collapse = "%"), "id", "t"))
      # if(to_wide){
      #   likelihood <- reshape(likelihood, idvar = "id", timevar = "t", direction = "wide")
      #   if(length(node) + 1  == ncol(preds)){
      #     setnames(preds, c("id", node))
      #   }
      # }
      # if(drop_id & "id" %in% colnames(likelihood)) likelihood$id <- NULL
      # if(drop_time & "t" %in% colnames(likelihood)) likelihood$t <- NULL
      # if(drop & ncol(likelihood) == 1) likelihood <- unlist(likelihood, use.names = F)
      #

      return(likelihood)
    },
    sample = function(tmle_task, n_samples = NULL, fold_number = "full", node = NULL, expand = T) {
      # TODO: fold
      # TODO when we have pooled tasks
      if(is.null(node)) node <- self$name

      if (is.null(tmle_task)) {
        tmle_task <- self$training_task
      }
      if (is.null(n_samples)) {
        return(tmle_task)
      }
      #TODO sampling will be messed up with degenerate likelihoods for some people
      #Need to sample the at_risk people and the not at_risk people separately.
      learner_task <- tmle_task$get_regression_task(node, expand = expand)
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
    },
    is_level_variant = function(){
      return(private$.is_level_variant)
    },
    include_degeneracy = function() {
      return(private$.include_degeneracy)
    }
  ),
  private = list(
    .name = NULL,
    .learner = NULL,
    .is_time_variant = NULL,
    .is_level_variant = NULL,
    .include_degeneracy = NULL
  )
)
