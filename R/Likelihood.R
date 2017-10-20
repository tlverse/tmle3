#' @importFrom R6 R6Class
#' @importFrom sl3 Lrnr_base args_to_list
#' @importFrom assertthat assert_that
#' @importFrom delayed bundle_delayed
Likelihood <- R6Class(classname = "Likelihood",
                      portable = TRUE,
                      class = TRUE,
                      inherit = Lrnr_base,
                      public = list(
                        initialize = function(learner_list, ...){
                          params <- args_to_list()
                          super$initialize(params=params, ...)
                        }
                      ),
                      active = list(
                        
                      ),
                      private = list(
                        .pretrain = function(task){
                          learner_list <- self$params$learner_list
                          npsem <- task$npsem
                          
                          nodes <- names(learner_list)
                          assert_that(all(nodes%in%names(npsem)))
                          
                          delayed_learner_fits <- lapply(nodes, function(node){
                            node_task <- task$get_regression_task(node)
                            delayed_learner_train(learner_list[[node]], node_task)
                          })
                          
                          return(bundle_delayed(delayed_learner_fits))
                        },
                        .train = function(task, learner_fits){
                          names(learner_fits) <- names(self$params$learner_list)
                          return(learner_fits)
                        },
                        .predict = function(task){
                          learner_fits <- private$.fit_object
                          nodes <- names(learner_fits)
                          n_to_pred <- task$nrow
                          n_learners <- length(nodes)
                          
                          ## Cannot use := to add columns to a null data.table (no columns),
                          ## hence we have to first seed an initial column, then delete it later
                          learner_preds <- data.table::data.table(init_seed_preds_to_delete = rep(NA_real_,n_to_pred))
                          
                          for(node in nodes) {
                            current_fit  <- learner_fits[[node]]
                            node_task <- task$get_regression_task(node)
                            current_preds <- current_fit$predict(node_task)
                            current_names <- node
                            if (!is.na(safe_dim(current_preds)[2]) && safe_dim(current_preds)[2] > 1) {
                              current_names <- paste0(learner_names[i], "_", names(current_preds))
                              stopifnot(length(current_names) == safe_dim(current_preds)[2])
                            }
                            data.table::set(learner_preds, j = current_names, value = current_preds)
                            invisible(NULL)
                          }
                          
                          ## remove the initial seeded column by reference
                          data.table::set(learner_preds, j = "init_seed_preds_to_delete", value = NULL)
                          return(learner_preds)
                        }
                      )
)