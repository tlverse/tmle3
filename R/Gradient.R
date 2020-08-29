#' Class representing a (possibly non-canonical) gradient for some parameter
#'
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#' @importFrom digest digest
#' @import data.table
#' @param Likelihood A trained likelihood object from which to compute the gradient from.
#' @param projection_task_generator A function that takes a task, likelihood, target parameter, and target_node
#' and returns a task with the outcome being the evaluated gradient,
#' and covariates being whatever variables the conditional density of the target node likelihood factor depends on.
#' (E.g. covariates + outcome in the regression task)
#' @export
#'
#' @keywords data
#'
#' @return \code{Gradient} object
#'
#' @format \code{\link{R6Class}} object.
#'
#'
#' @export
#
Gradient <- R6Class(
  classname = "Gradient",
  portable = TRUE,
  class = TRUE,
  inherit = Lrnr_base,
  public = list(
    initialize = function(Likelihood, projection_task_generator, target_param){
      params <- sl3::args_to_list()
      params$target_nodes <- target_param$update_nodes
      private$.params <- params
    },
    generate_task = function(tmle_task, node, include_outcome = T){
      self$projection_task_generator(task, self$Likelihood, self$target_param, node, outcome = include_outcome)
    },
    compute_EIC = function(tmle_task, node){
      fit_obj <- private$.component_fits[[node]]
      task <- self$generate_task(tmle_task, node, include_outcome = F)
      X <- task$X
      self$Likelihood$get_likelihood(tmle_task, node)

    },
    set_basis = function(){
      #Converts squashed basis to R functions of tmle3_tasks
      basis <- list()
      for(node in self$target_nodes){
        fit_obj <- private$.component_fits[[node]]
        basis_list <- fit_obj
        lst_to_func <- function(basis, node) {
          # Computes value of basis for 1 in task
          f <- function(task) {
            task <- self$generate_task(task, node, F)
            index <- which(names(task$X) == self$Likelihood$npsem[[node]]$variables)
            value <- as.numeric(all(unlist(task$X[, basis$cols, with = F], use.names = F) > basis$cutoffs))
            if(value == 0){ return(NULL)}
            return(basis$cutoffs[index])
            #return(as.vector(apply(task$X[, basis$cols, with = F], 1 , function(v){as.numeric(all(v > basis$cutoffs))})))
          }
          eval_lik <- function(y, fold_number = "full", task){
            #Generates likelihoods for outcome value y for task
            cf_data <- data.table(y)
            names(cf_data) <- node
            cf_task <- task$generate_counterfactual_task(UUIDgenerate(), cf_data)
            self$Likelihood$factor_list[[node]]$get_density(cf_task, fold_number, quick_pred = T, check_at_risk = F, expand = T)
          }
          eval_lik <- Vectorize(eval_lik, vectorize.args = "y")

          final_f <- function(task, fold_number){
           cutoff <- f(task)
           if(is.null(task)){
             return(0)
           }
           result <- 1 - integrate(eval_lik, lower = cutoff, upper = Inf, fold_number = fold_number, task = task)
          }
          return(final_f)
        }
        basis[node] <- lapply(basis_list, lst_to_func, node = node)
      }
      private$.basis <- basis
    },
    base_train = function(task, pretrain) {
      fit_object <- private$.train(task, pretrain)
      new_object <- self$clone() # copy parameters, and whatever else
      new_object$set_train(fit_object, task)
      return(new_object)
    }
  ),
  active = list(
    params = function(){
      private$.params
    },
    Likelihood = function(){
      private$.params$Likelihood
    },
    target_param = function(){
      private$.params$target_param
    },
    target_nodes = function(){
      private$.params$target_nodes
    },
    projection_task_generator = function(){
      private$.params$projection_task_generator
    },
    learner = function(){
      private$.learner
    },
    hal_args = function(args_to_add = NULL){
      if(!is.null(args_to_add)){
        for(name in names(args_to_add)){
          arg_value <- args_to_add[[name]]
          private$.learner_args[name] <- arg_value
        }
        args <- private$.learner_args
        args$learner_class <- Lrnr_hal9001
        private$.learner <- sl3:::call_with_args(make_learner, args)
      }
      return(private$.learner_args)
    }
  ),
  private = list(
    .train_sublearners = function(tmle_task){
      nodes <- self$target_nodes
      projected_fits <- lapply(nodes, function(node){
        task <- self$generate_task(task, node)
        lrnr <- self$learner$clone()
        return(lrnr$delayed_train(lrnr, task))
      })
      projected_fits <- bundle_delayed(projected_fits)
      return(projected_fits)
    },
    .train = function(tmle_task, projected_fits){
      #Store hal_fits
      component_fits <- lapply(projected_fits, `[[`, "fit_object")
      component_fits <- lapply(component_fits, squash_hal_fit)
      names(component_fits) <- self$target_nodes

      private$.component_fits <- component_fits
      self$set_basis()
    },
    .predict = function(tmle_task) {
      stop("This gradient has nothing to predict.")
    },
    .chain = function(tmle_task) {
      stop("This gradient has nothing to chain")
    },
    .params = NULL,
    .learner = make_learner(Lrnr_hal9001, max_degree = 2, family = "gaussian"),
    .learner_args = list(max_degree = 2, family = "gaussian"),
    .component_fits = list(),
    .basis = NULL
  )
)


Gradient_component =
