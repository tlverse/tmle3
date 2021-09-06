#' Nonparametric inference for user-specified parametric working models for the conditional odds ratio between two binary variables
#' The true conditional odds ratio is projected onto a parametric working model using logistic-regression.
#' This can be used to assess heterogeneity of the odds ratio.
#' We note that `formula_logOR = ~ 1` gives an estimator of the nonparametric marginal odds ratio among the treated.
#' The parametric model is at the log-scale and therefore the coefficients returned code the linear predictor for the `log`-conditional odds ratio.
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
#'     \item{\code{formula_logOR}}{...
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
Param_npOR <- R6Class(
  classname = "Param_npOR",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, formula_logOR = ~1, intervention_list_treatment, intervention_list_control, outcome_node = "Y") {
      super$initialize(observed_likelihood, list(), outcome_node)
      if (!is.null(observed_likelihood$censoring_nodes[[outcome_node]])) {
        # add delta_Y=0 to intervention lists
        outcome_censoring_node <- observed_likelihood$censoring_nodes[[outcome_node]]
        censoring_intervention <- define_lf(LF_static, outcome_censoring_node, value = 1)
        intervention_list_treatment <- c(intervention_list_treatment, censoring_intervention)
        intervention_list_control <- c(intervention_list_control, censoring_intervention)
      }
      private$.formula_logOR <- formula_logOR
      private$.cf_likelihood_treatment <- CF_Likelihood$new(observed_likelihood, intervention_list_treatment)
      private$.cf_likelihood_control <- CF_Likelihood$new(observed_likelihood, intervention_list_control)
    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full", is_training_task = TRUE) {
      training_task <- self$observed_likelihood$training_task
      if (is.null(tmle_task)) {
        tmle_task <- training_task
      }


      cf_task1 <- self$cf_likelihood_treatment$enumerate_cf_tasks(tmle_task)[[1]]
      cf_task0 <- self$cf_likelihood_control$enumerate_cf_tasks(tmle_task)[[1]]
      intervention_nodes <- union(names(self$intervention_list_treatment), names(self$intervention_list_control))

      W <- tmle_task$get_tmle_node("W")
      V <- model.matrix(self$formula_logOR, as.data.frame(W))
      A <- tmle_task$get_tmle_node("A", format = TRUE)[[1]]
      Y <- tmle_task$get_tmle_node("Y", format = TRUE)[[1]]
      W_train <- training_task$get_tmle_node("W")
      V_train <- model.matrix(self$formula_logOR, as.data.frame(W_train))
      A_train <- training_task$get_tmle_node("A", format = TRUE)[[1]]
      Y_train <- training_task$get_tmle_node("Y", format = TRUE)[[1]]

      g <- self$observed_likelihood$get_likelihoods(tmle_task, "A", fold_number)
      g1 <- ifelse(A == 1, g, 1 - g)
      g0 <- 1 - g1

      Q <- as.vector(self$observed_likelihood$get_likelihoods(tmle_task, "Y", fold_number))
      Q0 <- as.vector(self$cf_likelihood_treatment$get_likelihoods(cf_task0, "Y", fold_number))
      Q1 <- as.vector(self$cf_likelihood_treatment$get_likelihoods(cf_task1, "Y", fold_number))
      Qorig <- Q
      Q0 <- bound(Q0, 0.005)
      Q1 <- bound(Q1, 0.005)
      beta <- get_beta(W_train, A_train, self$formula_logOR, Q1, Q0, family = binomial(), weights = self$weights)

      Q1beta <- as.vector(plogis(qlogis(Q0) + V %*% beta))

      sigma_rel <- Q1beta * (1 - Q1beta) / (Q0 * (1 - Q0))


      omega <- as.vector((g0 + g1 * sigma_rel) / (g0))

      h_star <- -1 * as.vector((g1 * sigma_rel) / (g1 * sigma_rel + (1 - g1)))
      H <- as.matrix(omega * V * (A + h_star))

      # Store EIF component
      EIF_Y <- NULL
      EIFWA <- NULL
      if (is_training_task) {
        tryCatch(
          {
            scale <- apply(V, 2, function(v) {
              apply(self$weights * (A * Q1beta * (1 - Q1beta) * v * V), 2, mean)
            })
            scaleinv <- solve(scale)
            EIF_Y <- self$weights * (H %*% scaleinv) * as.vector(Y - Q)
            EIFWA <- apply(V, 2, function(v) {
              (self$weights * (A * v * (Q1 - Q1beta)) - mean(self$weights * (A * v * (Q1 - Q1beta))))
            })

            EIFWA <- EIFWA %*% scaleinv
          },
          error = function(...) {}
        )
      }

      return(list(Y = H, EIF = list(Y = EIF_Y, WA = EIFWA)))
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      cf_task1 <- self$cf_likelihood_treatment$enumerate_cf_tasks(tmle_task)[[1]]
      cf_task0 <- self$cf_likelihood_control$enumerate_cf_tasks(tmle_task)[[1]]

      W <- tmle_task$get_tmle_node("W")
      A <- tmle_task$get_tmle_node("A", format = TRUE)[[1]]
      Y <- tmle_task$get_tmle_node("Y", format = TRUE)[[1]]
      weights <- tmle_task$weights
      # clever_covariates happen here (for this param) only, but this is repeated computation
      EIF <- self$clever_covariates(tmle_task, fold_number, is_training_task = TRUE)$EIF
      EIF <- EIF$Y + EIF$WA
      Q <- self$observed_likelihood$get_likelihoods(tmle_task, "Y", fold_number)
      Q0 <- self$cf_likelihood_treatment$get_likelihoods(cf_task0, "Y", fold_number)
      Q1 <- self$cf_likelihood_treatment$get_likelihoods(cf_task1, "Y", fold_number)
      Qtest <- ifelse(A == 1, Q1, Q0)
      if (!all(Qtest - Q == 0)) {
        stop("Q and Q1,Q0 dont match")
      }
      # Q_packed <- sl3::unpack_predictions(self$observed_likelihood$get_likelihoods(tmle_task, "Y", fold_number))
      # Q0 <- Q_packed[[1]]
      # Q1 <- Q_packed[[2]]
      # Q <- Q_packed[[3]]
      Q0 <- bound(Q0, 0.0005)
      Q1 <- bound(Q1, 0.0005)
      beta <- get_beta(W, A, self$formula_logOR, Q1, Q0, family = binomial(), weights = weights)
      V <- model.matrix(self$formula_logOR, as.data.frame(W))
      OR <- exp(V %*% beta)

      IC <- EIF

      result <- list(psi = beta, IC = IC, OR = OR)
      return(result)
    }
  ),
  active = list(
    name = function() {
      param_form <- sprintf("log({P(Y=1|A=1,W)/P(Y=0|A=1,W)}/{P(Y=1|A=0,W)/P(Y=0|A=0,W)})")
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
    formula_logOR = function() {
      return(private$.formula_logOR)
    }
  ),
  private = list(
    .type = "OR",
    .cf_likelihood_treatment = NULL,
    .cf_likelihood_control = NULL,
    .supports_outcome_censoring = TRUE,
    .formula_logOR = NULL,
    .submodel = list(Y = "gaussian_identity")
  )
)
