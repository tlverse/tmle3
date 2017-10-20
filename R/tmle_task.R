#' @importFrom R6 R6Class
#' @importFrom sl3 sl3_Task
tmle_Task <- R6Class(classname = "tmle_Task",
                    portable = TRUE,
                    class = TRUE,
                    inherit = sl3_Task,
                    public = list(
                      initialize=function(data, tmle_nodes, npsem = NULL, ...){
                        private$.tmle_nodes <- tmle_nodes
                        
                        if(is.null(npsem)){
                          npsem <- list(W=c(),
                                        A=c("W"),
                                        Y=c("A","W"))
                          
                        }
                        private$.npsem <- npsem
                        
                        super$initialize(data, covariates = unlist(tmle_nodes), outcome=tmle_nodes$Y, ...)
                      },
                      get_regression_task=function(pred_node){
                        nodes <- self$tmle_nodes
                        npsem <- self$npsem
                        outcome <- nodes[[pred_node]]
                        parent_nodes <- npsem[[pred_node]]
                        covariates <- unlist(nodes[parent_nodes])
                        
                        # todo: transfer goodies like weights and ids and folds over
                        # todo: fix next_in_chain and maybe use that
                        return(sl3_Task$new(self$raw_data, covariates=covariates, outcome=outcome))
                      }
                    ),
                    active = list(
                      tmle_nodes=function(){
                        return(private$.tmle_nodes)
                      },
                      npsem = function(){
                        return(private$.npsem)
                      }
                        
                    ),
                    private = list(
                      .tmle_nodes = NULL,
                      .npsem = NULL
                    )
)
