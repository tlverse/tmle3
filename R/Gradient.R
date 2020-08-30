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
      private$.cache <- new.env()
    },
    generate_task = function(tmle_task, node, include_outcome = T){
      self$projection_task_generator(tmle_task, self$Likelihood, self$target_param, node, outcome = include_outcome)
    },
    expand_task = function(tmle_task, node, force = F){
      #Computes expanded task where observations are repeated (with fake ids) for all levels of node
      # TODO these expanded tasks should only be targetted for this node.
      if(tmle_task$uuid %in% names(private$.uuid_expanded_history)){
        if(private$.uuid_expanded_history[[tmle_task$uuid]] != node){
          stop("This expanded task does not match its node. You shouldn't be targeting this node for this task ")
        }
        return(tmle_task)
      }

      key <- paste0(tmle_task$uuid, node, sep = "%")

      cached_task <- get0(key, self$cache, inherits = FALSE)
      if(!is.null(cached_task)){
        return(cached_task)
      }
      variables <- tmle_task$npsem[[node]]$variables
      if(length(variables) >1) stop("Multivariate nodes not supported")
      data <- tmle_task$data
      data$trueid <- data$id
      levels <- sort(unique(unlist(data[, variables, with = F])))
      if(!force & length(levels) > 100){
        stop("Too many levels in node.")
      }
      long_data <- rbindlist(lapply(levels, function(level) {
        data <- copy(data)
        set(data ,, variables, level)
        return(data)
      }))
      long_data$id <-  paste(long_data$trueid,long_data[, variables, with = F][[1]], sep = "_")
      long_task <- tmle3_Task$new(long_data, tmle_task$npsem, id = "id", time = "t", force_at_risk = tmle_task$force_at_risk, summary_measure_columns = c(tmle_task$summary_measure_columns, "trueid"))

      assign(key, long_task, self$cache)
      setattr(long_task, "target_nodes", node)
      private$.uuid_expanded_history[[long_task$uiid]] <- node
      return(long_task)
    },
    compute_component = function(tmle_task, node, fold_number = "full"){
      self$assert_trained()
      #Converts squashed basis to R functions of tmle3_tasks

      fit_obj <- private$.component_fits[[node]]
      task <- self$generate_task(tmle_task, node, include_outcome = F)
      col_index <- which(colnames(task$X) == tmle_task$npsem[[node]]$variables )
      preds <- self$Likelihood$factor_list[[node]]$get_density(tmle_task, fold_number, quick_pred = T)
      type <- tmle_task$npsem[[node]]$variable_type$type
      if(type == "binomial"){
        preds <- data.table(cbind(1-preds, preds))
        setnames(preds, c("1", "2"))
        levels <- c(0,1)
      } else if (type == "categorical") {
        levels <- tmle_task$npsem[[node]]$variable_type$levels
        preds <- data.table(sl3::unpack_predictions(as.vector(preds)))
        setnames(preds, as.character(seq_along(levels)))
      } else{
        stop("Type not supported")
      }

      cdf <- data.table(t(apply(preds, 1, cumsum)))
      basis_list <- fit_obj$basis_list
      coefs <- as.vector(fit_obj$coefs)
      X <- as.matrix(task$X)
      design <- as.data.table(as.matrix(hal9001::make_design_matrix(X, basis_list)))

      diff_map <- unlist(lapply(seq_along(basis_list), function(i) {
        basis <- basis_list[[i]]
        if(!(col_index %in% basis$cols)){
          return(NULL)
        }
        result <- (list(which(levels == basis$cutoffs[which(basis$cols == col_index)])))
        names(result) = i
        return(result)
      }))

      center_basis <- lapply((names(diff_map)), function(i){
        col_index <- diff_map[[i]]
        diff <- design[[as.integer(i)]] - 1 + cdf[[col_index]]
        set(design, , as.integer(i), diff)
      })

      clean_basis <- function(basis){
        if(!(col_index %in% basis$cols)){
          return(basis)
        }
        index = which(basis$cols == col_index)
        basis$cutoffs[index] <- min(task$X[[node]]) - 1
        return(basis)
      }
      clean_list = lapply(basis_list, clean_basis)

      clean_design <- hal9001::make_design_matrix(X, clean_list)
      clean_design <- data.table(as.matrix(clean_design))

      #TODO only do this for basis functions containing y

      mid_result <- as.matrix(design * clean_design)
      result = coefs[1] + mid_result %*% coefs[-1]
      out = list(mat = mid_result, EIC = result)
      return(out)


    },
    base_train = function(task, pretrain) {
      fit_object <- private$.train(task, pretrain)
      new_object <- self$clone() # copy parameters, and whatever else
      new_object$set_train(fit_object, task)
      private$.training_task <- task
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
    basis = function(){
      private$.basis
    },
    training_task = function(){
      private$.training_task
    },
    cache = function(){
      private$.cache
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
        return(delayed_learner_train(lrnr, task))
      })
      projected_fits <- bundle_delayed(projected_fits)
      return(projected_fits)
    },
    .train = function(tmle_task, projected_fits){
      #Store hal_fits
      component_fits <- lapply(projected_fits, `[[`, "fit_object")
      component_fits <- lapply(component_fits, hal9001::squash_hal_fit)
      names(component_fits) <- self$target_nodes

      private$.component_fits <- component_fits
      return(component_fits)

    },
    .predict = function(tmle_task) {
      stop("This gradient has nothing to predict.")
    },
    .chain = function(tmle_task) {
      stop("This gradient has nothing to chain")
    },
    .params = NULL,
    .learner = make_learner(Lrnr_hal9001, max_degree = 1, family = "gaussian"),
    .learner_args = list(max_degree = 1, family = "gaussian"),
    .component_fits = list(),
    .basis = NULL,
    .training_task = NULL,
    .cache = NULL,
    .uuid_expanded_history = list()
  )
)



