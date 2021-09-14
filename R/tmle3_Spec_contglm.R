#' Defines a TML Estimator (except for the data)
#'
#' Current limitations: pretty much tailored to \code{Param_TSM}
#'
#' @importFrom R6 R6Class
#'
#' @export
#
tmle3_Spec_contglm <- R6Class(
  classname = "tmle3_Spec_contglm",
  portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(formula_continuous, formula_binary = formula_continuous, fA = function(x){x}, fA0 = fA, dfA = function(x){1}, dfA0 = dfA, include_A_binary = TRUE, submodel = NULL,
                          likelihood_override = NULL,
                          variable_types = NULL, delta_epsilon = 0.1, ...) {
      private$.options <- list(
        formula_continuous = formula_continuous, formula_binary = formula_binary,
        fA = fA, fA0 = fA0, dfA = dfA, dfA0 = dfA0,
        include_A_binary = include_A_binary,
        submodel = submodel,
        delta_epsilon = delta_epsilon,
        likelihood_override = likelihood_override,
        variable_types = variable_types, ...
      )
    },
    make_tmle_task = function(data, node_list, ...) {
      variable_types <- self$options$variable_types
      include_variance_node <- FALSE
      scale_outcome <- self$options$scale_outcome
      family <- self$options$submodel
      include_A_binary <- self$options$include_A_binary
      if(!is.null(family) && is.null(scale_outcome)) {
        scale_outcome <- FALSE
        if(family == "binomial") {
          scale_outcome <- TRUE
        }
      }

      Y <- data[[node_list$Y]]

      if (is.null(family)) {
        if (all(Y %in% c(0, 1))) {
          family <- "binomial"
          scale_outcome <- FALSE
        } else if (all(Y >= 0)) {
          family <- "poisson"
          scale_outcome <- FALSE
        } else {
          family <- "gaussian"
          scale_outcome <- FALSE
        }
      }
      private$.options$submodel <- family


      tmle_task <- point_tx_continuous_task(data, node_list, variable_types, scale_outcome = scale_outcome, include_A_binary = include_A_binary)

      return(tmle_task)
    },
    make_initial_likelihood = function(tmle_task, learner_list = NULL) {
      # Wrap baseline learner in semiparametric learner

      # produce trained likelihood when likelihood_def provided
      if (!is.null(self$options$likelihood_override)) {
        likelihood <- self$options$likelihood_override$train(tmle_task)
      } else {
        likelihood <- point_tx_continuous_likelihood(tmle_task, learner_list)
      }

      return(likelihood)
    },
    make_updater = function(convergence_type = "sample_size", verbose = TRUE, ...) {
      delta_epsilon <- self$options$delta_epsilon
      if (!is.null(self$options$verbose)) {
        verbose <- self$options$verbose
      }

      bounds <- list(Y = c(-Inf, Inf), A = 0.005)
      if (self$options$submodel == "poisson") {
        bounds <- list(Y = c(0.0025, Inf), A = 0.005)
      } else if (self$options$submodel == "binomial") {
        bounds <- list(Y = 0.0025, A = 0.005)
      }
      updater <- tmle3_Update$new(maxit = 100, one_dimensional = TRUE, delta_epsilon = 0.05, verbose = verbose, constrain_step = TRUE, bounds = bounds, ...)

      return(updater)
    },
    make_targeted_likelihood = function(likelihood, updater) {
      targeted_likelihood <- Targeted_Likelihood$new(likelihood, updater)
      return(targeted_likelihood)
    },
    make_params = function(tmle_task, targeted_likelihood) {
      formula_continuous <- self$options$formula_continuous
      formula_binary <- self$options$formula_binary
      family <- self$options$submodel

      param <- Param_contglm$new(targeted_likelihood, formula_CATE_binary = formula_binary, formula_CATE_continuous = formula_continuous,
                                 fA = self$options$fA,
                                 fA0 = self$options$fA0,
                                 dfA = self$options$dfA,
                                 dfA0 = self$options$dfA0,
                                 submodel = family)

      return(list(param))
    }
  ),
  active = list(
    options = function() {
      return(private$.options)
    }
  ),
  private = list(
    .options = NULL
  )
)
