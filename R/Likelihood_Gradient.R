#' @export
ipw_gen <- function(Y, A, g, W){
  Y*A/g  - Y*(1-A)/g
}

#' @export
generator_ate <-function(tmle_task, lik = NULL, target_param = NULL, node, outcome = T){
  task <- tmle_task$get_regression_task(node)
  A <- task$X$A
  Y <- task$Y
  W <- task$X$W

  g <- lik$get_likelihood(tmle_task, "A")

  IC <- ipw_gen(Y,A,g, W)

  cols <- task$add_columns(data.table(IC = IC))
  task <- task$clone()
  nodes <- task$nodes
  nodes$outcome <- "IC"
  nodes$covariates <- c(nodes$covariates, node)
  task$initialize(
    task$internal_data,
    nodes = nodes,
    folds = task$folds,
    column_names = cols,
    row_index = task$row_index,
    outcome_type = "continuous"
  )
  return(task)
  # task$next_in_chain(column_names = cols, covariates  = c(task$nodes$covariates, task$nodes$outcome), outcome = "IC")
}


#' Class representing a (possibly non-canonical) gradient for some parameter.
#' Currently, this gradient object takes a IPW/unobserved model gradient
#' and numerically projects it onto the tangent space.
#'
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#' @importFrom digest digest
#' @import data.table
#' @import sl3
#' @param likelihood A trained likelihood object from which to compute the gradient from.
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
    initialize = function(likelihood, ipw_args = NULL, projection_task_generator, target_nodes = "Y"){
      params <- sl3::args_to_list()
      params$target_nodes <- target_nodes
      private$.params <- params
      private$.cache <- new.env()
      private$.learner <- Lrnr_hal9001a$new(max_degree = 3, family = "gaussian")
    },
    train_projections = function(tmle_task, fold_number = private$.fold_number){
      # Based on the data of tmle_task, the L2(P_n) projection of the IPW gradient is computed onto the variable and past variabes for each node
      # This corresponds to a projection for each node to be targeted.
      private$.fold_number <- fold_number
      self$train(tmle_task)
    },
    generate_task = function(tmle_task, node, include_outcome = T, fold_number = "full"){
      #This generates the IPW projection task
      self$projection_task_generator(tmle_task, self$likelihood, node, include_outcome = include_outcome, self$params$ipw_args, fold_number = fold_number)
    },
    expand_task = function(tmle_task, node, force = F){
      #Computes expanded task where observations are repeated (with fake ids) for all levels of node
      # TODO these expanded tasks should only be targetted for this node.


      #  Cache's task for efficiency
      if(tmle_task$uuid %in% names(private$.uuid_expanded_history)){
        if(private$.uuid_expanded_history[[tmle_task$uuid]] != node){
          stop("This expanded task does not match its node. You shouldn't be targeting this node for this task ")
        }
        return(tmle_task)
      }

      key <- paste0(tmle_task$uuid, node, sep = "%")

      cached_task <- get0(key, self$cache, inherits = FALSE)
      if(!is.null(cached_task)){
        if(is.null(attr(cached_task, "target_nodes"))) {
          print(node)
          stop("wrong")
        }
        return(cached_task)
      }
      variables <- tmle_task$npsem[[node]]$variables

      if(length(variables) >1) stop("Multivariate nodes not supported")
      # Expand task by levels of node
      # TODO extend for continuous
      data <- tmle_task$data
      data$trueid <- data$id
      time <- tmle_task$npsem[[node]]$time
      levels <- sort(unique(unlist(tmle_task$get_tmle_node(node))))  #data[, variables, with = F])))

      if(!force & length(levels) > 100){
        stop("Too many levels in node.")
      }
      long_data <- rbindlist(lapply(levels, function(level) {
        data <- copy(data)
        set(data , which(data$t == time), variables, level)
        data$levelcopy <- level
        return(data)
      }))

      long_data$id <-  paste(long_data$trueid,long_data$levelcopy, sep = "_")

      suppressWarnings(long_task <- tmle3_Task$new(long_data, tmle_task$npsem, id = "id", time = "t", force_at_risk = tmle_task$force_at_risk, summary_measure_columns = c(tmle_task$summary_measure_columns, "trueid")))

      # By setting the "target_nodes" attribute to node, we ensure that this task is only targeted for this node
      #This is done for efficiency. If this likleihood values for this long_task is needed at a different node then this will lead to errors.
      setattr(long_task, "target_nodes", node)
      if(is.null(attr(long_task, "target_nodes"))) {
        print(node)
        stop("wrong")
      }
      assign(key, long_task, self$cache)

      private$.uuid_expanded_history[[long_task$uuid]] <- node
      return(long_task)
    },
    compute_component = function(tmle_task, node, fold_number = "full"){
      print("Gradient")
      print(node)
      print(fold_number)

      # Here, we compute the EIC component by first getting the IPW projection onto fewer variables from the hal fit and the nintegrating out.
      time <- tmle_task$npsem[[node]]$time

      self$assert_trained()
      #Converts squashed basis to R functions of tmle3_tasks

      fit_obj <- private$.component_fits[[node]]

      long_task <- self$expand_task(tmle_task, node)

      IC_task <- self$generate_task(tmle_task, node, include_outcome = F, fold_number = fold_number)


      col_index <- which(colnames(IC_task$X) == node )

      long_preds <- NULL

      # This gets the likelihood values of the long_task for node.
      # This trycatch stuff is only needed if we target tasks after already having targeted the likelihood.

      tryCatch({
        value_step1 <- self$likelihood$cache$get_update_step(self$likelihood$factor_list[[node]], tmle_task, fold_number, node = node)
        value_step2 <- self$likelihood$cache$get_update_step(self$likelihood$factor_list[[node]], long_task, fold_number, node = node)
        if(is.null(value_step1)){
          value_step1 <- 0
        }
        if(is.null(value_step2)){
          value_step2 <- 0
        }
        if(value_step1!=value_step2) {
          stop("Long_task and tmle_task out of sync.")
        }
        long_preds <- self$likelihood$get_likelihood(long_task, node, fold_number = fold_number, drop_id = T, drop_time = T, drop = T  )},
      error = function(e){
        #long_task is probably out of sync with tmle_task
        #Update it to the same level
        value_step <- self$likelihood$cache$get_update_step(self$likelihood$factor_list[[node]], tmle_task, fold_number, node = node)

        self$likelihood$sync_task(long_task, fold_number = fold_number, check = F, max_step = value_step)
        value_step1 <- self$likelihood$cache$get_update_step(self$likelihood$factor_list[[node]], tmle_task, fold_number, node = node)
        value_step2 <- self$likelihood$cache$get_update_step(self$likelihood$factor_list[[node]], long_task, fold_number, node = node)
        if(is.null(value_step1)){
          value_step1 <- 0
        }
        if(is.null(value_step2)){
          value_step2 <- 0
        }
        if(value_step1!=value_step2) {
          stop("Long_task and tmle_task out of sync.")
        }
        long_preds <<- self$likelihood$get_likelihood(long_task, node, fold_number = fold_number, drop_id = T, drop_time = T, drop = T, check_sync = F  )

      })


      data <- IC_task$data

      # TODO

      variables <- node
      #TODO check id order

      data <- data.table(cbind(merge(long_task$data[,c("id", "trueid", "t")], cbind(long_task$get_tmle_node(node, format = T, include_id = T, include_time = T), long_preds), by = c("id", "t"))
                               ))

      idkey <- data$trueid
      data$t <- NULL
      data$id <- NULL

      setnames(data, c("id", node, "pred"))

      setkeyv(data, cols = c("id", node))

      #TODO handle updating of expanded task
      #This is done by stacking copies of the cdf

      data <- dcast(data, as.formula(paste0("id ~ ", node)), value.var = "pred")
      id <- data$id
      data$id <- NULL
      levels <- as.numeric(colnames(data))

      # We need cdf to center the indicator basis functions.
      cdf <- as.data.table(t(apply(data, 1, cumsum)))
      setnames(cdf, as.character(levels))
      #print(cdf)

      if(long_task$uuid == tmle_task$uuid){
        #if expanded task is tmle_task then obtain then expand cdf to match
        #This ensures we dont have any recursion errors by expanding an expanded task

        match_index <- match(idkey, id)
        cdf <- cdf[match_index]
      }


      fit_obj <- private$.component_fits[[node]]
      basis_list <- fit_obj$basis_list
      coefs <- fit_obj$coefs
      col_index <- which(colnames(IC_task$X) == node )


      keep <- sapply(basis_list, function(b){
        col_index %in% b$cols
      })

      basis_list <- basis_list[keep]
      coefs <- coefs[c(T, keep)]

      #Should already be sorted
      X <- as.matrix(IC_task$X)
      # This is the design matrix generated by the basis which spans the initial IPW projection
      design <- as.data.table(as.matrix(hal9001::make_design_matrix(X, basis_list)))

      # For all levels of the node, we need to know what cutoff value for the node is used for each bsis function.
      # E.g if basis is 1(node > v) then this will give the index in "levels" which has the value v.

      basis_node_cutoff_map <- sapply(seq_along(basis_list), function(i) {
        basis <- basis_list[[i]]
        result <- (list(which(levels == basis$cutoffs[which(basis$cols == col_index)])))

        return(result)
      })

      # Based on the basis_node_cutoff_map index, for each basis function, we center it by subtracting its mean.
      # In the simple case of 1(node>v), we just need to do 1(node > v) - P(node > v), which is just the CDF.
      # This centers the design matrix in place.
      center_basis <- lapply(seq_along(basis_node_cutoff_map), function(i){
        col_index <- basis_node_cutoff_map[[i]]
        diff <- design[[as.integer(i)]] - 1 + cdf[[col_index]]
        set(design, , as.integer(i), diff)
      })
      # This is to make sure that basis functions which are zero (not because of the value of node) remain zero.
      # For example 1(a > v) 1(node >v) is centered as 1(a >v) ( 1 (node >v) - P(node >v|...)).
      # The center_basis actually computes 1(a >v)1 (node >v) - P(node >v|...), which is not correct.
      #We therefore need to multiply 1(a >v)1 (node >v) - P(node >v|...) by 1(a >v) to get 1(a >v) ( 1 (node >v) - P(node >v|...))
      # This is done by doing design * design_without_node

      min_val <- min(IC_task$X[[node]]) - 5
      get_basis_without_node <- function(basis){
        index = which(basis$cols == col_index)
        basis$cutoffs[index] <- min_val
        return(basis)
      }
      clean_list = lapply(basis_list, get_basis_without_node)
      #print(table(unlist( lapply(basis_list, `[[`, "cols"))))
      design_without_node <- hal9001::make_design_matrix(X, clean_list)
      design_without_node <- data.table(as.matrix(design_without_node))

      #print(as.data.table(design_without_node))
      #print(as.data.table(coefs))

      centered_design <- as.matrix(design * design_without_node)
      # Once we have our centered design matrix (centered_design), we can just multiply this by the hal basis coefficients.
      result =  centered_design %*% coefs[-1]
      out = list(col_index = col_index,Y = IC_task$Y, cdf = cdf,design = design,  centered_design = centered_design, coefs = coefs[-1], EIC = result)
      return(out)

    },
    compute_component_initial = function(tmle_task, node, fold_number = "full"){
      self$assert_trained()
      #Converts squashed basis to R functions of tmle3_tasks
      stop("no")
      fit_obj <- private$.component_fits[[node]]
      task <- self$generate_task(tmle_task, node, include_outcome = F)
      col_index <- which(colnames(task$X) == tmle_task$npsem[[node]]$variables )
      preds <- self$likelihood$factor_list[[node]]$get_density(tmle_task, fold_number, quick_pred = T)
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

      basis_node_cutoff_map <- unlist(lapply(seq_along(basis_list), function(i) {
        basis <- basis_list[[i]]
        if(!(col_index %in% basis$cols)){
          return(NULL)
        }
        result <- (list(which(levels == basis$cutoffs[which(basis$cols == col_index)])))
        names(result) = i
        return(result)
      }))

      center_basis <- lapply((names(basis_node_cutoff_map)), function(i){
        col_index <- basis_node_cutoff_map[[i]]
        diff <- design[[as.integer(i)]] - 1 + cdf[[col_index]]
        set(design, , as.integer(i), diff)
      })

      get_basis_without_node <- function(basis){
        if(!(col_index %in% basis$cols)){
          return(basis)
        }
        index = which(basis$cols == col_index)
        basis$cutoffs[index] <- min(task$X[[node]]) - 1
        return(basis)
      }
      clean_list = lapply(basis_list, get_basis_without_node)

      design_without_node <- hal9001::make_design_matrix(X, clean_list)
      design_without_node <- data.table(as.matrix(design_without_node))

      #TODO only do this for basis functions containing y

      centered_design <- as.matrix(design * design_without_node)
      result =  centered_design %*% coefs[-1]
      out = list(mat = centered_design, EIC = result)
      return(out)


    },
    base_train = function(task, pretrain) {
      fit_object <- private$.train(task, pretrain)
      #new_object <- self$clone() # copy parameters, and whatever else
      self$set_train(fit_object, task)
      private$.training_task <- task

      return(self)
    }
  ),
  active = list(
    params = function(){
      private$.params
    },
    likelihood = function(){
      private$.params$likelihood
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
    fold_number = function(){
      private$.fold_number
    },
    component_fits = function(){
      private$.component_fits
    },
    hal_args = function(args_to_add = NULL){
      if(!is.null(args_to_add)){
        for(name in names(args_to_add)){
          arg_value <- args_to_add[[name]]
          private$.learner_args[name] <- arg_value
        }
        args <- private$.learner_args
        #args$learner_class <- Lrnr_hal9001a
        private$.learner <- do.call(Lrnr_hal9001a$new, args)
      }
      return(private$.learner_args)
    }
  ),
  private = list(
    .train_sublearners = function(tmle_task){
      nodes <- c(self$target_nodes)

      projected_fits <- lapply(nodes, function(node){
        task <- self$generate_task(tmle_task, node, fold_number = self$fold_number)
        lrnr <- self$learner$clone()

        return(delayed_learner_train(lrnr, task))
      })

      projected_fits <- bundle_delayed(projected_fits)
      return(projected_fits)
    },
    .train = function(tmle_task, projected_fits){
      #Store hal_fits
      names(projected_fits) <- self$target_nodes

      component_fits <- lapply(projected_fits, `[[`, "fit_object")
      private$.fit_object <- projected_fits

      component_fits <- lapply(component_fits, function(fit){

        basis_list <- fit$basis_list[as.numeric(names(fit$copy_map))]
        coefs <- fit$coefs

        keep <- coefs[-1]!=0

        basis_list <- basis_list[keep]
        coefs_new <- coefs[c(T, keep)]
        if(sum(coefs_new) != sum(coefs)) stop("squash went wrong")
        if(length(coefs_new) != length(basis_list)+1) stop("squash went wrong")
        return(list(basis_list = basis_list, coefs = coefs_new))
      })
      names(component_fits) <- self$target_nodes

      private$.component_fits <- component_fits
      return(private$.fit_object)

    },
    .predict = function(tmle_task) {
      stop("This gradient has nothing to predict.")
    },
    .chain = function(tmle_task) {
      stop("This gradient has nothing to chain")
    },
    .params = NULL,
    .learner = NULL,
    .learner_args = list(max_degree = 3, family = "gaussian"),
    .component_fits = list(),
    .basis = NULL,
    .training_task = NULL,
    .cache = NULL,
    .uuid_expanded_history = list(),
    .fold_number = "full"
  )
)




#' @docType class
#' @importFrom R6 R6Class
#' @importFrom digest digest
#' @import hal9001
#' @import data.table
#' @import sl3
#' @export
Lrnr_hal9001a <- R6::R6Class(
  classname = "Lrnr_hal9001a", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(max_degree = 3,
                          fit_type = "glmnet",
                          n_folds = 10,
                          use_min = TRUE,
                          reduce_basis = NULL,
                          return_lasso = TRUE,
                          return_x_basis = FALSE,
                          basis_list = NULL,
                          cv_select = TRUE,
                          ...) {
      params <- args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("continuous", "binomial"),

    .train = function(task) {
      args <- self$params

      outcome_type <- self$get_outcome_type(task)

      if (is.null(args$family)) {
        args$family <- args$family <- outcome_type$glm_family()
      }

      args$X <- as.matrix(task$X)
      args$Y <- outcome_type$format(task$Y)
      args$yolo <- FALSE

      if (task$has_node("weights")) {
        args$weights <- task$weights
      }

      if (task$has_node("offset")) {
        args$offset <- task$offset
      }
      args$standardize <- F
      fit_object <- sl3:::call_with_args(hal9001::fit_hal, args)
      args$lambda <- fit_object$lambda_star
      args$lambda <- args$lambda * 0.15
      args$cv_select <- F
      if(args$fit_type == "glmnet"){
        fit_object <- sl3:::call_with_args(hal9001::fit_hal, args)
      }

      return(fit_object)
    },
    .predict = function(task = NULL) {
      predictions <- predict(self$fit_object, new_data = as.matrix(task$X))
      return(predictions)
    },
    .required_packages = c("hal9001")
  )
)




