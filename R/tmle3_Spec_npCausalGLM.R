#' Defines a TML Estimator (except for the data)
#'
#' Current limitations: pretty much tailored to \code{Param_TSM}
#'
#' @importFrom R6 R6Class
#'
#' @export
#
tmle3_Spec_npCausalGLM <- R6Class(
  classname = "tmle3_Spec_npCausalGLM",
  portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(formula, estimand = c("CATE", "CATT", "TSM", "OR", "RR"), treatment_level = 1, control_level = 0, submodel = NULL,
                          likelihood_override = NULL,
                          variable_types = NULL, delta_epsilon = 0.025, ...) {
      estimand <- match.arg(estimand)
      private$.options <- list(
        estimand = estimand, formula = formula, submodel = submodel,
        treatment_level = treatment_level, control_level = control_level, delta_epsilon = delta_epsilon,
        likelihood_override = likelihood_override,
        variable_types = variable_types, ...
      )
    },
    make_tmle_task = function(data, node_list, ...) {
      variable_types <- self$options$variable_types
      include_variance_node <- FALSE
      scale_outcome <- TRUE
      Y <- data[[node_list$Y]]
      family <- self$options$submodel

      if (is.null(family) && self$options$estimand %in% c("CATE", "CATT", "TSM")) {
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
      } else if (is.null(family) && self$options$estimand == "RR") {
        if (all(Y %in% c(0, 1))) {
          family <- "binomial"
        } else {
          family <- "poisson"
          scale_outcome <- FALSE
        }
      } else if (!is.null(family)) {
        if (family == "binomial") {
          scale_outcome <- TRUE
        } else {
          scale_outcome <- FALSE
        }
      }
      private$.options$submodel <- family
      binary_outcome <- all(data[[node_list$Y]] %in% c(0, 1))
      private$.options$binary_outcome <- binary_outcome
      if (self$options$estimand == "RR") {
        if (binary_outcome) {
          type <- "binomial"
        } else {
          type <- "continuous"
        }
        variable_types <- list(Y = variable_type(type))
        # scale_outcome <- binary_outcome
      } else if (self$options$estimand == "OR") {
        variable_types <- list(Y = variable_type("binomial"))
      }

      tmle_task <- point_tx_task(data, node_list, variable_types, scale_outcome = scale_outcome, include_variance_node = include_variance_node)

      return(tmle_task)
    },
    make_initial_likelihood = function(tmle_task, learner_list = NULL) {
      # Wrap baseline learner in semiparametric learner

      # produce trained likelihood when likelihood_def provided
      if (!is.null(self$options$likelihood_override)) {
        likelihood <- self$options$likelihood_override$train(tmle_task)
      } else {
        likelihood <- point_tx_likelihood(tmle_task, learner_list)
      }

      return(likelihood)
    },
    make_updater = function(convergence_type = "sample_size", verbose = TRUE, ...) {
      delta_epsilon <- self$options$delta_epsilon
      if (!is.null(self$options$verbose)) {
        verbose <- self$options$verbose
      }
      if (self$options$estimand == "CATE" || self$options$estimand == "CATT" || self$options$estimand == "TSM") {
        updater <- tmle3_Update$new(maxit = 100, one_dimensional = FALSE, delta_epsilon = 1, verbose = verbose, constrain_step = FALSE, bounds = c(-Inf, Inf), ...)
      } else if (self$options$estimand == "OR") {
        updater <- tmle3_Update$new(maxit = 200, one_dimensional = TRUE, convergence_type = convergence_type, verbose = verbose, delta_epsilon = delta_epsilon, constrain_step = TRUE, bounds = 0.0025, ...)
      } else if (self$options$estimand == "RR") {
        if (self$options$submodel == "poisson") {
          bounds <- list(Y = c(0.0025, Inf), A = 0.005)
        } else {
          bounds <- list(Y = 0.0025, A = 0.005)
        }


        updater <- tmle3_Update$new(maxit = 200, one_dimensional = TRUE, convergence_type = convergence_type, verbose = verbose, delta_epsilon = delta_epsilon, constrain_step = TRUE, bounds = bounds, ...)
      }
      return(updater)
    },
    make_targeted_likelihood = function(likelihood, updater) {
      targeted_likelihood <- Targeted_Likelihood$new(likelihood, updater)
      return(targeted_likelihood)
    },
    make_params = function(tmle_task, targeted_likelihood) {
      treatment_value <- self$options$treatment_level
      control_value <- self$options$control_level
      formula <- self$options$formula
      family <- self$options$submodel
      A_levels <- tmle_task$npsem[["A"]]$variable_type$levels
      if (!is.null(A_levels)) {
        treatment_value <- factor(treatment_value, levels = A_levels)
        control_value <- factor(control_value, levels = A_levels)
      }



      if (self$options$estimand == "TSM") {
        # If TSM generate params for all levels
        param <- lapply(union(treatment_value, control_value), function(value) {
          treatment <- define_lf(LF_static, "A", value = value)
          return(Param_npTSM$new(targeted_likelihood, formula, treatment, submodel = family))
        })
        return(param)
      } else {
        treatment <- define_lf(LF_static, "A", value = treatment_value)
        control <- define_lf(LF_static, "A", value = control_value)
      }

      if (self$options$estimand == "CATE") {
        param <- Param_npCATE$new(targeted_likelihood, formula, treatment, control, submodel = family)
      } else if (self$options$estimand == "CATT") {
        param <- Param_npCATT$new(targeted_likelihood, formula, treatment, control, submodel = family)
      } else if (self$options$estimand == "OR") {
        param <- Param_npOR$new(targeted_likelihood, formula, treatment, control)
      } else if (self$options$estimand == "RR") {
        param <- Param_npRR$new(targeted_likelihood, formula, treatment, control, binary_outcome = self$options$binary_outcome, submodel = family)
      }
      return(list(param))
    }
  ),
  active = list(
    options = function() {
      return(private$.options)
    },
    family = function() {
      return(private$.families[[self$options$estimand]])
    }
  ),
  private = list(
    .options = NULL,
    .families = list("CATE" = gaussian(), "RR" = poisson(), "OR" = binomial())
  )
)
