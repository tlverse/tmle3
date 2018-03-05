#' Shifted Likelihood Factor
#'
#' Shifts a likelihood factor according to a \code{shift_function}. 
#' In effect, \code{get_likelihood(tmle_task)} will instead the likelihood from the \code{original_lf},
#' but for shifted value \eqn{A'=}\code{shift_function}\eqn{(A,W)}
#' 
#' @references 
#' Díaz, Iván, and Mark J van der Laan. 2017. “Stochastic Treatment Regimes.” In Targeted Learning in Data Science: Causal Inference for Complex Longitudinal Studies, 167–80. Springer Science & Business Media.
#' Muñoz, Iván Díaz, and Mark J van der Laan. 2012. “Population Intervention Causal Effects Based on Stochastic Interventions.” Biometrics 68 (2). Wiley Online Library: 541–49. 
#'
#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @family Likelihood objects
#' @keywords data
#'
#' @return \code{LF_base} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @section Constructor:
#'   \code{define_lf(LF_shift, name, type = "density", original_lf, shift_function, ...)}
#'
#'   \describe{
#'     \item{\code{name}}{character, the name of the factor. Should match a node name in the nodes specified by \code{\link{tmle3_Task}$tmle_nodes}
#'     }
#'     \item{\code{type}}{character, either "density", for conditional density or, "mean" for conditional mean
#'     }
#'     \item{\code{original_lf}}{\code{\link{LF_base}} object, the likelihood factor to shift
#'     }
#'     \item{\code{shift_function}}{\code{function}, defines the shift
#'     }
#'     \item{\code{shift_inverse}}{\code{function}, the inverse of shift_function
#'     }
#'     \item{\code{...}}{Not currently used.
#'     }
#'     }
#'
#' @section Fields:
#' \describe{
#'     \item{\code{original_lf}}{\code{\link{LF_base}} object, the likelihood factor to shift
#'     }
#'     \item{\code{shift_function}}{\code{function}, defines the shift
#'     }
#'     \item{\code{shift_inverse}}{\code{function}, the inverse of shift_function
#'     }
#'     }
#' @export
LF_shift <- R6Class(
  classname = "LF_shift",
  portable = TRUE,
  class = TRUE,
  inherit = LF_base,
  public = list(
    initialize = function(name, type="density", original_lf, shift_function, shift_inverse, ...) {
      super$initialize(name, type)
      private$.original_lf <- original_lf
      private$.shift_function <- shift_function
      private$.shift_inverse <- shift_inverse
    },
    get_mean = function(tmle_task) {
      stop("get_mean not supported for LF_shift")
    },
    get_likelihood = function(tmle_task) {
      #get shifted data
      shifted_values <- self$shift_inverse(tmle_task)
      
      #generate cf_task data
      cf_data <- data.table(shifted_values)
      setnames(cf_data, self$name)
      
      cf_task <- tmle_task$generate_counterfactual_task(UUIDgenerate(), cf_data)
      
      #get original likelihood for shifted data
      cf_likelihood <- self$original_lf$get_likelihood(cf_task)

      return(cf_likelihood)
    },
    cf_values = function(tmle_task) {
      cf_values <- self$shift_function(tmle_task)
      return(cf_values)
    }
  ),
  active = list(
    original_lf = function() {
      return(private$.original_lf)
    },
    shift_function = function() {
      return(private$.shift_function)
    },
    shift_inverse = function() {
      return(private$.shift_inverse)
    }
    
  ),
  private = list(
    .name = NULL,
    .original_lf = NULL,
    .shift_function = NULL,
    .shift_inverse = NULL
  )
)
