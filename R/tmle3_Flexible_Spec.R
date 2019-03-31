#' Defines a TML Estimator (except for the data)
#'
#' Tries to do this a bit better than tmle3_Spec
#'
#' @importFrom R6 R6Class
#'
#' @export
#
tmle3_Flexible_Spec <- R6Class(
  classname = "tmle3_Flexible_Spec",
  portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(make_tmle_task = NULL, 
                          make_initial_likelihood = NULL, 
                          make_updater = NULL, 
                          make_targeted_likelihood = NULL,
                          make_params = NULL, ...) {
      private$.override_list = list(make_tmle_task = make_tmle_task,
                           make_initial_likelihood = make_initial_likelihood,
                           make_updater = make_updater,
                           make_targeted_likelihood = make_targeted_likelihood,
                           make_params = make_params,
                           ...)
    },
    make = function(component_name){
      fun_name <- sprintf("make_%s", component_name)
      step_fun <- private$.override_list[[step_name]]
      if(is.null(override_fun)){
        step_fun <- self[[step_name]])
      }
      result <- step_fun(self)
      private$.component_list[[component_name]] <- result
      return(result)
    },
    do_spec = function(...){
      tmle_task <- self$make("tmle_task")
      initial_likelihood <- self$make("initial_likelihood")
      updater <- self$make("updater")
      targeted_likelihood <- self$make("targeted_likelihood")
      tmle_params <- tmle_spec$make("params")
      
      fit <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)
      return(fit)
    },
    make_tmle_task = function(spec) {
      tmle_task <- point_tx_task(data, node_list)
      return(tmle_task)
    },
    make_initial_likelihood = function(spec) {
      # produce trained likelihood when likelihood_def provided

      if (!is.null(self$options$likelihood_override)) {
        likelihood <- self$options$likelihood_override$train(tmle_task)
      } else {
        likelihood <- point_tx_likelihood(tmle_task, learner_list)
      }
      
      return(likelihood)
    },
    make_updater = function(spec) {
      updater <- tmle3_Update$new()
    },
    make_targeted_likelihood = function(spec) {
      targeted_likelihood <- Targeted_Likelihood$new(likelihood, updater)
      return(targeted_likelihood)
    },
    make_params = function(spec) {
      stop("this is a base class, try tsm_Spec_TSM_all")
    },
    component = function(component_name, component_step=NULL){
      component <- private$.component_list[[component_name]]
      if(is.null(component)){
        if(is.null(component_step)){
          component_step = sprintf("make_%s",component_name)
        }
        stop(sprintf("%s not available. Do %s step first", component_name, component_step))
      }
      
      return(component)
    }
  ),
  active = list(
    options = function() {
      return(private$.options)
    },
    tmle_task = function(){
      return(self$component(tmle_task))
    },
    initial_likelihood = function(){
      return(self$component(initial_likelihood))      
    },
    updater = function(){
      return(self$component(updater))            
    },
    targeted_likelihood = function(){
      return(self$component(targeted_likelihood))
    },
    params = function(){
      return(self$component(params))
    }
  ),
  private = list(
    .override_list = list(),
    .component_list = list()
  )
)
