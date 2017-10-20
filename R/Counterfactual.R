#' @importFrom R6 R6Class
Counterfactual <- R6Class(classname = "Counterfactual",
                          portable = TRUE,
                          class = TRUE,
                          public = list(
                            initialize = function(npsem){
                              private$.npsem=npsem
                            }),
                          active = list(
                            npsem = function(){
                              return(private$.npsem)
                            }
                          ),
                          private = list(
                            .npsem = NULL
                          )
)