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
    initialize = function(formula, estimand = c("CATE", "CATT", "OR", "RR"), treatment_level = 1, control_level =0,
                          likelihood_override = NULL,
                          variable_types = NULL, ...) {
      estimand <- match.arg(estimand)
      private$.options <- list(estimand = estimand, formula = formula,
                               treatment_level = treatment_level, control_level = control_level,
                               likelihood_override = likelihood_override,
                               variable_types = variable_types, ...
      )
    },
    make_tmle_task = function(data, node_list, ...) {
      variable_types <- self$options$variable_types
      include_variance_node <- FALSE

      tmle_task <- point_tx_task(data, node_list, variable_types, scale_outcome = FALSE, include_variance_node = include_variance_node)

      return(tmle_task)
    },
    make_initial_likelihood = function(tmle_task, learner_list = NULL ) {
      #Wrap baseline learner in semiparametric learner

      # produce trained likelihood when likelihood_def provided
      if (!is.null(self$options$likelihood_override)) {
        likelihood <- self$options$likelihood_override$train(tmle_task)
      } else {
        likelihood <- point_tx_likelihood(tmle_task, learner_list)
      }

      return(likelihood)
    },
    make_updater = function(convergence_type = "sample_size", verbose = TRUE,...) {
      if(self$options$estimand == "CATE" || self$options$estimand == "CATT"){
        updater <- tmle3_Update$new(maxit=100,one_dimensional = FALSE,   verbose = verbose, constrain_step = FALSE, bounds = c(-Inf, Inf), ...)
      } else if (self$options$estimand == "OR"){
        updater <- tmle3_Update$new(maxit = 200, one_dimensional = TRUE, convergence_type = convergence_type, verbose = verbose,delta_epsilon = 0.001, constrain_step = TRUE, bounds = 0.0025, ...)
      } else if (self$options$estimand == "RR"){
        updater <- tmle3_Update$new(maxit = 200, one_dimensional = TRUE,  convergence_type = convergence_type, verbose = verbose, delta_epsilon = 0.001, constrain_step = TRUE, bounds = c(0.0025, Inf), ...)
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
      A_levels <- tmle_task$npsem[["A"]]$variable_type$levels
      if (!is.null(A_levels)) {
        treatment_value <- factor(treatment_value, levels = A_levels)
        control_value <- factor(control_value, levels = A_levels)
      }
      treatment <- define_lf(LF_static, "A", value = treatment_value)
      control <- define_lf(LF_static, "A", value = control_value)
      formula <- self$options$formula
      if(self$options$estimand == "CATE"){
        param <- Param_npCATE$new(targeted_likelihood,formula, treatment, control)
      } else if(self$options$estimand == "CATT"){
        param <- Param_npCATT$new(targeted_likelihood,formula, treatment, control)
      } else if (self$options$estimand == "OR"){
        param <- Param_npOR$new(targeted_likelihood,formula, treatment, control)
      } else if (self$options$estimand == "RR"){
        param <- Param_npRR$new(targeted_likelihood, formula, treatment, control)
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
