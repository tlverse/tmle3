#' Defines a TML Estimator for the Marginal Mean Outcome under a Candidate Adaptive Design
#'
#' @importFrom R6 R6Class
#' @importFrom tmle3 tmle3_Spec define_lf tmle3_Update Targeted_Likelihood
#'  Param_TSM point_tx_likelihood
#'
#' @export
tmle3_Spec_marginal_ADSM <- R6::R6Class(
  classname = "tmle3_Spec_marginal_ADSM",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Spec,
  public = list(
    initialize = function(treatment_level, control_level, g_treat, g_adapt, ...) {
      super$initialize(
        treatment_level = treatment_level,
        control_level = control_level,
        g_treat = g_treat,
        g_adapt = g_adapt, ...
      )
    },
    make_params = function(tmle_task, likelihood, ...) {
      g_treat <-self$options$g_treat
      g_adapt <- self$options$g_adapt
      if (!(is.vector(g_treat) &
            is.vector(g_adapt) &
            tmle_task$nrow == length(g_treat) &
            tmle_task$nrow == length(g_adapt))) {
        msg <- paste("`g_treat` and `g_adapt` must be vectors",
                     "with a length of `tmle_task$nrow`")
        stop(msg)
      }

      treatment_value <- self$options$treatment_level
      control_value <- self$options$control_level
      A_levels <- tmle_task$npsem[["A"]]$variable_type$levels
      if (!is.null(A_levels)) {
        treatment_value <- factor(treatment_value, levels = A_levels)
        control_value <- factor(control_value, levels = A_levels)
      }
      treatment <- define_lf(LF_static, "A", value = treatment_value)
      control <- define_lf(LF_static, "A", value = control_value)
      adsm <- Param_marginal_ADSM$new(likelihood, treatment, control, g_treat, g_adapt)
      tmle_params <- list(adsm)
      return(tmle_params)
    },
    make_updater = function() {
      updater <- tmle3_Update$new(cvtmle = TRUE)
    }
  ),
  active = list(),
  private = list()
)

################################################################################

#' Mean Outcome under a Candidate Adaptive Design (Marginalized Over Covaraite Distribution)
#'
#' O = (W, A, Y)
#' W = Covariates
#' A = Treatment (binary or categorical)
#' Y = Outcome (binary or bounded continuous)
#'
#' @importFrom sl3 make_learner Lrnr_mean
#' @param treatment_level the level of A that corresponds to treatment
#' @param control_level the level of A that corresponds to a control or reference level
#' @param g_treat the actual probability of A that corresponds to treatment
#' @param g_adapt the probability of A that corresponds to treatment under a candidate adaptive design
#' @export
tmle_marginal_ADSM <- function(treatment_level, control_level, g_treat, g_adapt) {
  # TODO: unclear why this has to be in a factory function
  tmle3_Spec_marginal_ADSM$new(treatment_level, control_level, g_treat, g_adapt)
}
