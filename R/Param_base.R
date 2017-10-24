#' @importFrom R6 R6Class
#' @export
Param_base <- R6Class(classname = "Param_base",
                          portable = TRUE,
                          class = TRUE,
                          public = list(
                            initialize = function(outcome_node){
                              private$.outcome_node <- outcome_node
                            },
                            HA = function(likelihood, task){
                              stop("Param_base is a base class")
                            },
                            estimates = function(likelihood, task){
                              stop("Param_base is a base class")
                            }
                          ),
                          active = list(
                            outcome_node = function(){
                              return(private$.outcome_node)
                            }
                          ),
                          private = list(
                            .outcome_node = NULL
                          )
)