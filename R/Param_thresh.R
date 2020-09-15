#' Estimates the threshold function E_W E[Y | A >=c, W] for a range of given values threshold values c
#'
#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @importFrom dplyr near
#' @family Parameters
#' @keywords data
#'
#' @return \code{Param_base} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @section Constructor:
#'   \code{define_param(Param_ATT, observed_likelihood, intervention_list, ..., outcome_node)}
#'
#'   \describe{
#'     \item{\code{observed_likelihood}}{A \code{\link{Likelihood}} corresponding to the observed likelihood
#'     }
#'     \item{\code{intervention_list_treatment}}{A list of objects inheriting from \code{\link{LF_base}}, representing the treatment intervention.
#'     }
#'     \item{\code{intervention_list_control}}{A list of objects inheriting from \code{\link{LF_base}}, representing the control intervention.
#'     }
#'     \item{\code{...}}{Not currently used.
#'     }
#'     \item{\code{outcome_node}}{character, the name of the node that should be treated as the outcome
#'     }
#'     }
#'

#' @section Fields:
#' \describe{
#'     \item{\code{cf_likelihood_treatment}}{the counterfactual likelihood for the treatment
#'     }
#'     \item{\code{cf_likelihood_control}}{the counterfactual likelihood for the control
#'     }
#'     \item{\code{intervention_list_treatment}}{A list of objects inheriting from \code{\link{LF_base}}, representing the treatment intervention
#'     }
#'     \item{\code{intervention_list_control}}{A list of objects inheriting from \code{\link{LF_base}}, representing the control intervention
#'     }
#' }
#' @export
Param_thresh <- R6Class(
  classname = "Param_thresh",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, cutoffs, thresh_node = "A", outcome_node = "Y") {
      super$initialize(observed_likelihood, list(), outcome_node = outcome_node)
      cf_task <- observed_likelihood$training_task
      # cf_data where everyone is the maximum level of of the node so that they are above threshhold in every group
      cf_task$data
      cf_data <- data.table(rep(2*max(cutoffs), cf_task$nrow))

      setnames(cf_data, thresh_node)

      cf_data$id <- cf_task$id
      cf_data$t <- cf_task$time

      cf_task <- cf_task$generate_counterfactual_task(UUIDgenerate(), cf_data)
      #cache task


      observed_likelihood$get_likelihood(cf_task, "Y")

      private$.censoring_node <- (observed_likelihood$censoring_nodes[[outcome_node]])
      private$.thresh_node <- thresh_node
      private$.cutoffs <- cutoffs
      private$.strict_threshold <- F
      private$.cf_task <- cf_task

    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      thresh_node <- private$.thresh_node
      censoring_node <- private$.censoring_node
      cutoffs <- private$.cutoffs

      cdfS <- as.vector(self$observed_likelihood$get_likelihood(tmle_task, thresh_node, fold_number))

      cdfS <- bound(cdfS, c(0.0005, .9995))
      S <- tmle_task$get_tmle_node(thresh_node)

      if(!is.null(censoring_node)) {
        pCensoring <- self$observed_likelihood$get_likelihood(tmle_task, censoring_node, fold_number)
        uncensored <- as.vector(as.numeric(as.numeric(tmle_task$get_tmle_node(censoring_node))==1))
      } else {
        pCensoring <- 1
        uncensored <- 1
      }
      if(private$.strict_threshold) {
        indS <- as.vector(unlist(lapply(cutoffs, function(cutoff) {as.numeric(S >= cutoff)})))
      } else {
        indS <- as.vector(unlist(lapply(cutoffs, function(cutoff) {as.numeric(S >= cutoff)})))
      }
      #Uses

      HA <- indS * uncensored / (pCensoring * (1-cdfS))
      n = tmle_task$nrow
      k = length(cutoffs)
      H <- matrix(0, nrow = length(HA), ncol = k)
      for(i in 1:k){
        first <- (i-1)*n + 1
        last <- (i)*n
        H[first:last,i] <- HA[first:last]
      }

      if(any(!dplyr::near(rowSums(H),HA))){
        stop("oops")
      }
      if(length(indS)!= length(cdfS)) {
        stop("Uneven lengths in cdfS and indS")
      }
      return(list(Y = H))
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      thresh_node <- private$.thresh_node
      censoring_node <- private$.censoring_node
      cutoffs <- private$.cutoffs
      cf_task <- private$.cf_task
      intervention_nodes <- union(names(self$intervention_list_treatment), names(self$intervention_list_control))

      # clever_covariates happen here (for this param) only, but this is repeated computation
      HA <- matrix(rowSums(self$clever_covariates(tmle_task, fold_number)[[self$outcome_node]]), nrow = tmle_task$nrow)

      Y <- matrix(tmle_task$get_tmle_node(self$outcome_node, impute_censoring = TRUE), nrow = tmle_task$nrow)

      #get E[Y|A>=1, W]
      EY <- matrix(self$observed_likelihood$get_likelihood(tmle_task, self$outcome_node, fold_number), nrow = tmle_task$nrow)

      #get E[Y|A>=1, W]
      EY1 <- matrix(self$observed_likelihood$get_likelihood(cf_task, self$outcome_node, fold_number), nrow = tmle_task$nrow)

      psi <- colMeans(EY1)


      IC <- HA * (as.vector(Y) - EY)  + t((t(EY1)  - psi))
      weights <- tmle_task$get_regression_task(self$outcome_node)$weights
      result <- list(psi = psi, IC = IC * weights)
      return(result)
    }
  ),
  active = list(
    name = function() {
      param_form <- sprintf("ATE[%s_{%s}-%s_{%s}]", self$outcome_node, self$cf_likelihood_treatment$name, self$outcome_node, self$cf_likelihood_control$name)
      return(param_form)
    },
    cf_task = function() {
      return(private$.cf_task)
    },

    update_nodes = function() {
      return(c(self$outcome_node))
    }
  ),
  private = list(
    .type = "Threshold",
    .cf_likelihood_treatment = NULL,
    .cf_likelihood_control = NULL,
    .supports_outcome_censoring = TRUE,
    .submodel_type_supported = c("logistic"),
    .supports_weights = T,
    .censoring_node = NULL,
    .thresh_node = NULL,
    .cutoffs = NULL,
    .strict_threshold = NULL,
    .cf_task = NULL
  )
)
