#' Defines a Stratified TML Estimator with MSM (except for the data)
#'
#' @importFrom R6 R6Class
#' @importFrom assertthat assert_that
#'
#' @export
#
tmle3_Spec_MSM <- R6Class(
  classname = "tmle3_Spec_MSM",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Spec,
  public = list(
    initialize = function(mass = NULL, ...) {
      super$initialize(
        mass = mass, ...
      )
    },
    make_tmle_task = function(data, node_list, ...) {
      assert_that(data.class(data[[node_list$V]]) == "numeric",
                  msg = "Stratified variable should be numeric.")
      
      private$.strata_variable = node_list$V
      # Initial estimate should include V as covariate
      # Note: Upon current structure, the easiest way is to add V into W.
      #       This won't cause problem in calculation, 
      #       but nodes dependency graph will be off.
      node_list$W = c(node_list$W, node_list$V)
      
      super$make_tmle_task(data, node_list, ...)
    },
    make_params = function(tmle_task, targeted_likelihood) {
      treatment_type <- variable_type(x=tmle_task$get_tmle_node("A"))$type
      
      if (treatment_type == "continuous") {
        tmle_params <- define_param(Param_MSM, targeted_likelihood, self$strata_variable, 
                                    continuous_treatment = TRUE, n_samples = 30, 
                                    mass = self$options$mass)
      } else {
        A_vals <- tmle_task$get_tmle_node("A")
        if (is.factor(A_vals)) {
          A_levels <- levels(A_vals)
          A_levels <- factor(A_levels, A_levels)
        } else {
          A_levels <- sort(unique(A_vals))
        }
        tmle_params <- define_param(Param_MSM, targeted_likelihood, self$strata_variable, 
                                    continuous_treatment = FALSE, treatment_values = A_levels, 
                                    mass = self$options$mass)
      }
      return(tmle_params)
    }
  ),
  active = list(
    strata_variable = function() {
      return(private$.strata_variable)
    }
  ),
  private = list(
    .strata_variable = NULL
  )
)

#' Make MSM version of Stratified TML estimator class
#'
#' O=(W,A,Y)
#' W=Covariates
#' A=Treatment (binary or categorical)
#' Y=Outcome (binary or bounded continuous)
#'
#' @importFrom sl3 make_learner Lrnr_mean
#'
#' @param mass h(A, V)
#'
#' @export
tmle_MSM <- function(mass = NULL) {
  # TODO: unclear why this has to be in a factory function
  tmle3_Spec_MSM$new(mass = mass)
}
