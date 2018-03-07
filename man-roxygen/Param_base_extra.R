#' @section Constructor:
#'   \code{define_param(Param_base, observed_likelihood, ..., outcome_node)}
#'                      
#'   \describe{
#'     \item{\code{observed_likelihood}}{A \code{\link{Likelihood}} corresponding to the observed likelihood
#'     }
#'     \item{\code{...}}{Not currently used.
#'     }
#'     \item{\code{outcome_node}}{character, the name of the node that should be treated as the outcome
#'     }
#'     }
#'     
#' @section Methods:
#'
#' \describe{
#' \item{\code{clever_covariates(tmle_task = NULL)}}{
#'   Get the clever covariates for an TMLE update step.
#'
#'   \itemize{
#'     \item{\code{tmle_task}: \code{\link{tmle3_Task}} to get clever covariate values for. 
#'     If NULL, the \code{tmle_task} used to train the \code{observed likelihood will be used}
#'     }
#'   }
#'   }
#' \item{\code{estimates(tmle_task = NULL)}}{
#'   Get the parameter estimates and influence curve values.
#'
#'   \itemize{
#'     \item{\code{tmle_task}: \code{\link{tmle3_Task}} to get clever covariate values for. 
#'     If NULL, the \code{tmle_task} used to train the \code{observed likelihood will be used}
#'     }
#'   }
#'   }
#'}
#'
#' @section Fields:
#' \describe{
#'     \item{\code{observed_likelihood}}{the observed likelihood
#'     }
#'     \item{\code{outcome_node}}{character, the name of the outcome node
#'     }
#'}
