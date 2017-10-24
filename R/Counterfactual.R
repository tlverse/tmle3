#' @importFrom R6 R6Class
#' @importFrom sl3 Lrnr_base args_to_list
#' @export
Counterfactual <- R6Class(classname = "Counterfactual",
                          portable = TRUE,
                          class = TRUE,
                          inherit = Lrnr_base,
                          public = list(
                            initialize = function(intervention_list, name = NULL){
                              if(inherits(intervention_list, "LF_base")){
                                intervention_list <- list(intervention_list)
                              }
                              
                              intervention_nodes <- sapply(intervention_list, `[[`, "name")
                              names(intervention_list) <- intervention_nodes
                              
                              private$.name = name
                              private$.intervention_list = intervention_list
                            },
                            cf_likelihood = function(likelihood){
                              return(likelihood$modify_factors(self$intervention_list))
                            },
                            cf_task = function(task){
                              cf_task <- task

                              for(current_factor in self$intervention_list){
                                
                                node_name <- current_factor$name
                                
                                if(!current_factor$is_degenerate){
                                  stop("all intervention nodes must be degenerate (static or dynamic rule nodes, no stochastic nodes)")
                                }
                                # for current_factor, generate counterfactual values
                                node_variable <- cf_task$tmle_nodes[[node_name]]$variables
                                node_values <- current_factor$get_values(cf_task)
                                
                                # create new tmle task with counterfactual values
                                new_data <- data.table(node_values)
                                setnames(new_data, names(new_data), node_variable)
                                cf_task <- cf_task$generate_counterfactual_task(self$fit_uuid, new_data)  
                              }
                              
                              return(cf_task)
                            }
                            
                          ),
                          active = list(
                            name = function(){
                              return(private$.name)
                            },
                            intervention_list = function(){
                              return(private$.intervention_list)
                            },
                            intervention_nodes = function(){
                              return(names(self$intervention_list))
                            }
                          ),
                          private = list(
                            .name = NULL,
                            .intervention_list = NULL
                          )
)

#' @export
define_cf <- function(intervention_list, name = NULL){
  return(Counterfactual$new(intervention_list, name))
}