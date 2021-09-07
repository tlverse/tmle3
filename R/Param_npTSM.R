#' Nonparametric inference for user-specified parametric working models for the conditional treatment effect.
#' The true conditional average treatment effect is projected onto a parametric working model using least-squares regression.
#' Unlike \code{Param_npCATT}, this function uses all observations to compute the projection.
#' This can be used to assess heterogeneity of the average treatment effect.
#' We note that `formula_TSM = ~ 1` gives an estimator of the nonparametric average treatment effect (ATE).
#'
#' Parameter definition for the Average Treatment Effect (ATE).
#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
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
#'     \item{\code{formula_TSM}}{...
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
Param_npTSM <- R6Class(
  classname = "Param_npTSM",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, formula_TSM = ~1, intervention_list, family_fluctuation = c("binomial", "gaussian", "poisson"), outcome_node = "Y") {
      family_fluctuation <- match.arg(family_fluctuation)
      private$.submodel <- list(Y = family_fluctuation)
      super$initialize(observed_likelihood, list(), outcome_node)
      training_task <- self$observed_likelihood$training_task
      W <- training_task <- self$observed_likelihood$training_task$get_tmle_node("W")
      V <- model.matrix(formula_TSM, as.data.frame(W))
      private$.formula_names <- colnames(V)
      private$.targeted <- rep(T, ncol(V))

      if (!is.null(observed_likelihood$censoring_nodes[[outcome_node]])) {
        # add delta_Y=0 to intervention lists
        outcome_censoring_node <- observed_likelihood$censoring_nodes[[outcome_node]]
        censoring_intervention <- define_lf(LF_static, outcome_censoring_node, value = 1)
        intervention_list <- c(intervention_list, censoring_intervention)
      }
      private$.formula_TSM <- formula_TSM
      private$.cf_likelihood <- CF_Likelihood$new(observed_likelihood, intervention_list)
    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full", is_training_task = TRUE) {
      training_task <- self$observed_likelihood$training_task
      if (is.null(tmle_task)) {
        tmle_task <- training_task
      }


      cf_task1 <- self$cf_likelihood$enumerate_cf_tasks(tmle_task)[[1]]

      W <- tmle_task$get_tmle_node("W")
      V <- model.matrix(self$formula_TSM, as.data.frame(W))
      A <- tmle_task$get_tmle_node("A", format = T)[[1]]
      Y <- tmle_task$get_tmle_node("Y")
      W_train <- training_task$get_tmle_node("W")
      V_train <- model.matrix(self$formula_TSM, as.data.frame(W_train))
      A_train <- training_task$get_tmle_node("A", format = TRUE)[[1]]
      Y_train <- training_task$get_tmle_node("Y", format = TRUE)[[1]]

      intervention_nodes <- names(self$intervention_list)
      pA <- self$observed_likelihood$get_likelihoods(tmle_task, intervention_nodes, fold_number)
      cf_pA <- self$cf_likelihood$get_likelihoods(tmle_task, intervention_nodes, fold_number)

      if (!is.null(ncol(pA)) && ncol(pA) > 1) {
        pA <- apply(pA, 1, prod)
      }
      if (!is.null(ncol(cf_pA)) && ncol(cf_pA) > 1) {
        cf_pA <- apply(cf_pA, 1, prod)
      }

      Q <- as.vector(self$observed_likelihood$get_likelihoods(tmle_task, "Y", fold_number))
      Q1 <- as.vector(self$cf_likelihood$get_likelihoods(cf_task1, "Y", fold_number))

      beta1 <- coef(glm.fit(as.matrix(V_train), Q1, family = gaussian(), weights = self$weights, intercept = F))
      Q1beta <- as.vector(V %*% beta1)

      H <- V * (cf_pA / pA)

      EIF_Y <- NULL
      # Store EIF component
      if (is_training_task) {
        scale <- apply(V, 2, function(v) {
          apply(self$weights * V * (v), 2, mean)
        })

        scaleinv <- solve(scale)
        EIF_Y <- self$weights * (H %*% scaleinv) * as.vector(Y - Q)
        EIF_WA <-
          apply(V, 2, function(v) {
            self$weights * (v * (Q1 - Q1beta) - mean(self$weights * v * (Q1 - Q1beta)))
          }) %*% scaleinv
      }


      return(list(Y = H, EIF = list(Y = EIF_Y, WA = EIF_WA)))
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      cf_task1 <- self$cf_likelihood$enumerate_cf_tasks(tmle_task)[[1]]

      W <- tmle_task$get_tmle_node("W")
      V <- model.matrix(self$formula_TSM, as.data.frame(W))
      A <- tmle_task$get_tmle_node("A", format = T)[[1]]
      Y <- tmle_task$get_tmle_node("Y")

      weights <- tmle_task$weights
      # clever_covariates happen here (for this param) only, but this is repeated computation
      EIF <- self$clever_covariates(tmle_task, fold_number, is_training_task = TRUE)$EIF
      EIF <- EIF$Y + EIF$WA
      Q <- self$observed_likelihood$get_likelihoods(tmle_task, "Y", fold_number)
      Q1 <- self$cf_likelihood$get_likelihoods(cf_task1, "Y", fold_number)
      beta1 <- coef(glm.fit(V, Q1, family = gaussian(), weights = weights, intercept = F))

      IC <- as.matrix(EIF)

      result <- list(psi = beta1, IC = IC)
      return(result)
    }
  ),
  active = list(
    name = function() {
      param_form <- unlist(paste0(sprintf("E[%s_{%s}]: ", self$outcome_node, self$cf_likelihood$name), private$.formula_names)) # sprintf("CATE[%s_{%s}-%s_{%s}]", self$outcome_node, self$cf_likelihood_treatment$name, self$outcome_node, self$cf_likelihood_control$name)
      return(param_form)
    },
    cf_likelihood = function() {
      return(private$.cf_likelihood)
    },
    intervention_list = function() {
      return(self$cf_likelihood$intervention_list)
    },
    update_nodes = function() {
      return(c(self$outcome_node))
    },
    formula_TSM = function() {
      return(private$.formula_TSM)
    }
  ),
  private = list(
    .type = "TSM",
    .cf_likelihood = NULL,
    .supports_outcome_censoring = TRUE,
    .formula_TSM = NULL,
    .submodel = list(Y = "binomial_identity"),
    .formula_names = NULL
  )
)
