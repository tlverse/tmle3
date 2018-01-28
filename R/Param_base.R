#' Base Class for Parameter Estimation
#'
#' @importFrom R6 R6Class
#'
#' @export
#
Param_base <- R6Class(
  classname = "Param_base",
  portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(outcome_node) {
      private$.outcome_node <- outcome_node
    },
    clever_covariates = function(likelihood, task) {
      stop("Param_base is a base class")
    },
    estimates = function(likelihood, task) {
      stop("Param_base is a base class")
    },
    default_submodel = NULL,
    default_loss = NULL,
    default_dag = NULL
  ),
  active = list(
    outcome_node = function() {
      return(private$.outcome_node)
    },
    dag = function() {
      return(private$.dag)
    },
    submodel = function() {
      return(private$.submodel)
    },
    loss = function() {
      return(private$.loss)
    }
  ),
  private = list(
    .outcome_node = NULL,
    .dag = NULL,
    .submodel = NULL,
    .loss = NULL
  )
)
