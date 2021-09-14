#' Nonparametric inference for user-specified parametric working models for the conditional treatment effect.
#' The true conditional average treatment effect is projected onto a parametric working model using least-squares regression.
#' Unlike \code{Param_npCATT}, this function uses all observations to compute the projection.
#' This can be used to assess heterogeneity of the average treatment effect.
#' We note that `formula_CATE = ~ 1` gives an estimator of the nonparametric average treatment effect (ATE).
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
#'     \item{\code{formula_CATE}}{...
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
Param_npCATE <- R6Class(
  classname = "Param_npCATE",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, formula_CATE = ~1, intervention_list_treatment, intervention_list_control, submodel = c("binomial", "gaussian", "poisson"), outcome_node = "Y") {
      submodel <- match.arg(submodel)
      super$initialize(observed_likelihood, list(), outcome_node, submodel = submodel)
      training_task <- self$observed_likelihood$training_task
      W <- training_task <- self$observed_likelihood$training_task$get_tmle_node("W")
      V <- model.matrix(formula_CATE, as.data.frame(W))
      private$.formula_names <- colnames(V)
      private$.targeted <- rep(T, ncol(V))



      if (!is.null(observed_likelihood$censoring_nodes[[outcome_node]])) {
        # add delta_Y=0 to intervention lists
        outcome_censoring_node <- observed_likelihood$censoring_nodes[[outcome_node]]
        censoring_intervention <- define_lf(LF_static, outcome_censoring_node, value = 1)
        intervention_list_treatment <- c(intervention_list_treatment, censoring_intervention)
        intervention_list_control <- c(intervention_list_control, censoring_intervention)
      }
      private$.formula_CATE <- formula_CATE
      private$.cf_likelihood_treatment <- CF_Likelihood$new(observed_likelihood, intervention_list_treatment)
      private$.cf_likelihood_control <- CF_Likelihood$new(observed_likelihood, intervention_list_control)
    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full", is_training_task = F) {
      training_task <- self$observed_likelihood$training_task
      if (is.null(tmle_task)) {
        tmle_task <- training_task
      }


      cf_task1 <- self$cf_likelihood_treatment$enumerate_cf_tasks(tmle_task)[[1]]
      cf_task0 <- self$cf_likelihood_control$enumerate_cf_tasks(tmle_task)[[1]]
      intervention_nodes <- union(names(self$intervention_list_treatment), names(self$intervention_list_control))
      pA <- self$observed_likelihood$get_likelihoods(tmle_task, intervention_nodes, fold_number)
      cf_pA_treatment <- self$cf_likelihood_treatment$get_likelihoods(tmle_task, intervention_nodes, fold_number)
      cf_pA_control <- self$cf_likelihood_control$get_likelihoods(tmle_task, intervention_nodes, fold_number)
      # collapse across multiple intervention nodes
      if (!is.null(ncol(cf_pA_treatment)) && ncol(cf_pA_treatment) > 1) {
        cf_pA_treatment <- apply(cf_pA_treatment, 1, prod)
      }
      if (!is.null(ncol(cf_pA_control)) && ncol(cf_pA_control) > 1) {
        cf_pA_control <- apply(cf_pA_control, 1, prod)
      }
      if (!is.null(ncol(pA)) && ncol(pA) > 1) {
        pA <- apply(pA, 1, prod)
      }


      W <- tmle_task$get_tmle_node("W")
      V <- model.matrix(self$formula_CATE, as.data.frame(W))

      Y <- tmle_task$get_tmle_node("Y")
      W_train <- training_task$get_tmle_node("W")
      V_train <- model.matrix(self$formula_CATE, as.data.frame(W_train))

      # g <- self$observed_likelihood$get_likelihoods(tmle_task, "A", fold_number)
      # g1 <- ifelse(A == 1, g, 1 - g)
      # g0 <- 1 - g1

      Q <- as.vector(self$observed_likelihood$get_likelihoods(tmle_task, "Y", fold_number))
      Q0 <- as.vector(self$cf_likelihood_treatment$get_likelihoods(cf_task0, "Y", fold_number))
      Q1 <- as.vector(self$cf_likelihood_treatment$get_likelihoods(cf_task1, "Y", fold_number))
      beta <- coef(glm.fit(V_train, Q1 - Q0, family = gaussian(), weights = self$weights))
      CATE <- as.vector(V_train %*% beta)
      # var_Y <- self$cf_likelihood_treatment$get_likelihoods(tmle_task, "var_Y", fold_number)
      # var_Y0 <- self$cf_likelihood_treatment$get_likelihoods(cf_task0, "var_Y", fold_number)
      # var_Y1 <- self$cf_likelihood_treatment$get_likelihoods(cf_task1, "var_Y", fold_number)


      H <- V / pA * (cf_pA_treatment - cf_pA_control)

      EIF_Y <- NULL
      EIF_WA <- NULL
      # Store EIF component
      if (is_training_task) {
        scale <- apply(V, 2, function(v) {
          apply(self$weights * V * (v), 2, mean)
        })

        scaleinv <- solve(scale)
        EIF_Y <- self$weights * (H %*% scaleinv) * as.vector(Y - Q)
        EIF_WA <- apply(V, 2, function(v) {
          self$weights * (v * (Q1 - Q0 - CATE) - mean(v * self$weights * (Q1 - Q0 - CATE)))
        }) %*% scaleinv

        # print(dim(EIF_Y))
        # print(mean(EIF_Y))
      }


      return(list(Y = H, EIF = list(Y = EIF_Y, WA = EIF_WA)))
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      cf_task1 <- self$cf_likelihood_treatment$enumerate_cf_tasks(tmle_task)[[1]]
      cf_task0 <- self$cf_likelihood_control$enumerate_cf_tasks(tmle_task)[[1]]

      W <- tmle_task$get_tmle_node("W")
      V <- model.matrix(self$formula_CATE, as.data.frame(W))

      # A <- tmle_task$get_tmle_node("A", format = T)[[1]]
      Y <- tmle_task$get_tmle_node("Y")

      weights <- tmle_task$weights
      # clever_covariates happen here (for this param) only, but this is repeated computation
      EIF_list <- self$clever_covariates(tmle_task, fold_number, is_training_task = TRUE)$EIF
      EIF <- EIF_list
      EIF <- EIF$Y + EIF$WA
      Q <- self$observed_likelihood$get_likelihoods(tmle_task, "Y", fold_number)
      Q0 <- self$cf_likelihood_treatment$get_likelihoods(cf_task0, "Y", fold_number)
      Q1 <- self$cf_likelihood_treatment$get_likelihoods(cf_task1, "Y", fold_number)


      beta <- coef(glm.fit(V, Q1 - Q0, family = gaussian(), weights = self$weights))

      CATE <- Q1 - Q0

      IC <- as.matrix(EIF)

      result <- list(psi = beta, IC = IC, IC_list = EIF_list, CATE = CATE)
      return(result)
    }
  ),
  active = list(
    name = function() {
      param_form <- private$.formula_names # sprintf("CATE[%s_{%s}-%s_{%s}]", self$outcome_node, self$cf_likelihood_treatment$name, self$outcome_node, self$cf_likelihood_control$name)
      return(param_form)
    },
    cf_likelihood_treatment = function() {
      return(private$.cf_likelihood_treatment)
    },
    cf_likelihood_control = function() {
      return(private$.cf_likelihood_control)
    },
    intervention_list_treatment = function() {
      return(self$cf_likelihood_treatment$intervention_list)
    },
    intervention_list_control = function() {
      return(self$cf_likelihood_control$intervention_list)
    },
    update_nodes = function() {
      return(c(self$outcome_node))
    },
    formula_CATE = function() {
      return(private$.formula_CATE)
    }
  ),
  private = list(
    .type = "CATE",
    .cf_likelihood_treatment = NULL,
    .cf_likelihood_control = NULL,
    .supports_outcome_censoring = TRUE,
    .formula_CATE = NULL,
    .formula_names = NULL
  )
)
