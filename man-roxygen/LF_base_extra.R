#' @section Constructor:
#'   \code{define_lf(LF_base, name, ..., type = "density")}
#'                      
#'   \describe{
#'     \item{\code{name}}{character, the name of the factor. Should match a node name in the nodes specified by \code{\link{tmle3_Task}$npsem}
#'     }
#'     \item{\code{...}}{Not currently used.
#'     }
#'     \item{\code{type}}{character, either "density", for conditional density or, "mean" for conditional mean
#'     }
#'     }
#'     
#' @section Methods:
#'
#' \describe{
#' \item{\code{get_density(tmle_task)}}{
#'   Get conditional density values for for the observations in \code{tmle_task}.
#'
#'   \itemize{
#'     \item{\code{tmle_task}: \code{\link{tmle3_Task}} to get likelihood values for
#'     }
#'   }
#'   }
#'\item{\code{get_mean(tmle_task)}}{
#'   Get conditional mean values for for the observations in \code{tmle_task}.
#'   \itemize{
#'     \item{\code{tmle_task}: \code{\link{tmle3_Task}} to get likelihood values for
#'     }
#'   }
#'   }
#'}
#'
#' @section Fields:
#' \describe{
#'     \item{\code{name}}{character, the name of the factor. Should match a node name in the nodes specified by \code{\link{tmle3_Task}$npsem}
#'     }
#'     \item{\code{type}}{character, either "density", for conditional density or, "mean" for conditional mean
#'     }
#'     \item{\code{variable_type}}{\code{\link[sl3]{variable_type}} object, specifying the data type of the outcome variable. Only available after Likelihood training.
#'     }
#'     \item{\code{values}}{Possible values of the outcome variable, retrivied from the \code{variable_type} object. Only available after Likelihood training.
#'     }
#'}
