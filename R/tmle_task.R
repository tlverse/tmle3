#' @importFrom R6 R6Class
#' @importFrom sl3 sl3_Task
#' @export
tmle_Task <- R6Class(classname = "tmle_Task",
                    portable = TRUE,
                    class = TRUE,
                    inherit = sl3_Task,
                    public = list(
                      initialize=function(data, tmle_nodes, ...){
                        super$initialize(data, covariates = c(), outcome=NULL, ...)
                        
                        node_names <- sapply(tmle_nodes, `[[`, "name")
                        names(tmle_nodes) <- node_names
                        for(node_name in node_names){
                          variables <- tmle_nodes[[node_name]]$variables
                          variable_data <- super$get_data(,variables)
                          if(ncol(variable_data)==1){
                            variable_data <- unlist(variable_data, use.names = FALSE)
                          } 
                          if(is.null(tmle_nodes[[node_name]]$variable_type)){
                            tmle_nodes[[node_name]]$guess_variable_type(variable_data)
                          }
                        }
                        
                        private$.tmle_nodes <- tmle_nodes
                        y_task <- self$get_regression_task("Y", data)
                        
                      },
                      get_regression_task=function(target_node, data = NULL){
                        nodes <- self$tmle_nodes
                        target_node <- nodes[[target_node]]
                        outcome <- target_node$variables
                        parent_names <- target_node$parents
                        parent_nodes <- nodes[parent_names]
                        covariates <- unlist(lapply(parent_nodes,`[[`, "variables"))
                        if(is.null(data)){
                          data <- self$raw_data
                        }
                        # todo: transfer goodies like weights and ids and folds over
                        # todo: fix next_in_chain and maybe use that
                        # todo: subset based on censoring nodes, etc
                        return(sl3_Task$new(data, covariates=covariates, outcome=outcome, 
                                            outcome_type=target_node$variable_type,
                                            column_names = self$column_names))
                      },
                      generate_counterfactual_task=function(uuid, new_data){
                        new_task <- self$clone()
                        new_column_names <- new_task$add_columns(uuid, new_data)
                        new_task$initialize(self$raw_data, self$tmle_nodes, column_names = new_column_names)
                        return(new_task)
                      },
                      next_in_chain = function(...){
                        return(super$next_in_chain(tmle_nodes = self$tmle_nodes, ...))
                      }
                      
                    ),
                    active = list(
                      tmle_nodes=function(){
                        return(private$.tmle_nodes)
                      }
                    ),
                    private = list(
                      .tmle_nodes = NULL
                    )
)

