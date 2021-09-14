#' Nonparametric inference for user-specified linear working models for user-specified conditional continuous treatment effect contrasts.
#' Treatment must be nonnegative with `A=0` being the control assignment.
#' Specifically, the user-specified model is of the form `fA(E[Y|A,W]) - fA0(E[Y|A=0,W])  = 1(A >0) * formula_binary(W) + A * (formula_binary)`
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
Param_contglm <- R6Class(
  classname = "Param_contglm",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, formula_CATE_binary = ~1, formula_CATE_continuous = ~1, fA, fA0, dfA, dfA0, transform = NULL,  submodel = c("binomial", "gaussian", "poisson"), outcome_node = "Y") {
      submodel <- match.arg(submodel)
      super$initialize(observed_likelihood, list(), outcome_node, submodel = submodel)
      training_task <- self$observed_likelihood$training_task
      W <- training_task <- self$observed_likelihood$training_task$get_tmle_node("W")
      Vbin <- model.matrix(formula_CATE_binary, as.data.frame(W))
      Vcont <- model.matrix(formula_CATE_continuous, as.data.frame(W))
      private$.formula_names <- list(bin = colnames(Vbin), cont = colnames(Vcont))
      private$.targeted <- rep(T, ncol(Vbin))
      private$.custom_functions <- list(fA = fA, fA0 = fA0, dfA = dfA, dfA0 = dfA0, transform = transform)

      intervention_list <- list()
      if (!is.null(observed_likelihood$censoring_nodes[[outcome_node]])) {
        # add delta_Y=0 to intervention lists
        outcome_censoring_node <- observed_likelihood$censoring_nodes[[outcome_node]]
        censoring_intervention <- define_lf(LF_static, outcome_censoring_node, value = 1)
        intervention_list <- c(censoring_intervention)
      }
      private$.formula_CATE_binary <- formula_CATE_binary
      private$.formula_CATE_continuous <- formula_CATE_continuous
      private$.cf_likelihood <- CF_Likelihood$new(observed_likelihood, intervention_list)
    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full", is_training_task = FALSE) {
      training_task <- self$observed_likelihood$training_task
      if (is.null(tmle_task)) {
        tmle_task <- training_task
      }



      intervention_nodes <- names(self$intervention_list) # Just censoring
      pA1 <- self$observed_likelihood$get_likelihoods(tmle_task, "A_binary", fold_number)
      pA0 <- 1 - pA1

      EA <- self$observed_likelihood$get_likelihoods(tmle_task, "A", fold_number)
      if (length(intervention_nodes) > 0) {
        G <- self$observed_likelihood$get_likelihoods(tmle_task, intervention_nodes[[1]], fold_number)
        Delta <- self$cf_likelihood$get_likelihoods(tmle_task, intervention_nodes[[1]], fold_number)
      } else {
        G <- 1
        Delta <- 1
      }




      W <- tmle_task$get_tmle_node("W")

      V_binary <- model.matrix(self$formula_CATE_binary, as.data.frame(W))
      V_cont <- model.matrix(self$formula_CATE_continuous, as.data.frame(W))
      A_binary <- tmle_task$get_tmle_node("A_binary")

      A <- tmle_task$get_tmle_node("A")
      Y <- tmle_task$get_tmle_node("Y")
      W_train <- training_task$get_tmle_node("W")
      A_binary_train <- training_task$get_tmle_node("A_binary")
      A_train <- training_task$get_tmle_node("A")
      V_train_binary <- model.matrix(self$formula_CATE_binary, as.data.frame(W_train))
      V_train_cont <- model.matrix(self$formula_CATE_continuous, as.data.frame(W_train))

      cftask0 <- tmle_task$generate_counterfactual_task(UUIDgenerate(), data.table(A = rep(0, tmle_task$nrow), A_binary = rep(0, tmle_task$nrow)))

      Q <- as.vector(self$observed_likelihood$get_likelihoods(tmle_task, "Y", fold_number))
      Q0 <- as.vector(self$observed_likelihood$get_likelihoods(cftask0, "Y", fold_number))

      V_full <- cbind(A_binary_train * V_train_binary, A_binary_train * A_train * V_train_cont)
      custom_funs <- self$custom_functions

      # var_Y <- self$cf_likelihood_treatment$get_likelihoods(tmle_task, "var_Y", fold_number)
      # var_Y0 <- self$cf_likelihood_treatment$get_likelihoods(cf_task0, "var_Y", fold_number)
      # var_Y1 <- self$cf_likelihood_treatment$get_likelihoods(cf_task1, "var_Y", fold_number)


      H_binary <- V_binary * (A_binary * custom_funs$dfA(Q) - (1 - A_binary) * pA1 / pA0 * custom_funs$dfA0(Q0) )
      H_cont <- V_cont * (A_binary * A * custom_funs$dfA(Q)  - (1 - A_binary) * EA / pA0 * custom_funs$dfA0(Q0))
      H <- Delta / G * cbind(H_binary, H_cont)
      EIF_Y <- NULL
      EIF_WA <- NULL
      # Store EIF component
      if (is_training_task) {
        true_contrast <- custom_funs$fA(Q) - custom_funs$fA0(Q0)
        beta <- coef(glm.fit(V_full, true_contrast, family = gaussian(), weights = self$weights, intercept = FALSE))
        contrast <- as.vector(V_full %*% beta)

        scale <- apply(V_full, 2, function(v) {
          apply(self$weights * V_full * (v), 2, mean)
        })

        scaleinv <- solve(scale)
        EIF_Y <- self$weights * (H %*% scaleinv) * as.vector(Y - Q)
        EIF_WA <- apply(V_full, 2, function(v) {
          self$weights * (v * (true_contrast - contrast) - mean(v * self$weights * (true_contrast - contrast)))
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

      custom_funs <- self$custom_functions
      W <- tmle_task$get_tmle_node("W")
      V_binary <- model.matrix(self$formula_CATE_binary, as.data.frame(W))
      V_cont <- model.matrix(self$formula_CATE_continuous, as.data.frame(W))
      A_binary <- tmle_task$get_tmle_node("A_binary")
      A <- tmle_task$get_tmle_node("A")



      weights <- tmle_task$weights
      # clever_covariates happen here (for this param) only, but this is repeated computation
      IC_list <- self$clever_covariates(tmle_task, fold_number, is_training_task = TRUE)$EIF
      EIF <- IC_list
      EIF <- EIF$Y + EIF$WA
      cftask0 <- tmle_task$generate_counterfactual_task(UUIDgenerate(), data.table(A = rep(0, tmle_task$nrow), A_binary = rep(0, tmle_task$nrow)))
      Q <- as.vector(self$observed_likelihood$get_likelihoods(tmle_task, "Y", fold_number))
      Q0 <- as.vector(self$observed_likelihood$get_likelihoods(cftask0, "Y", fold_number))

      V_full <- cbind(A_binary * V_binary, A_binary * A * V_cont)
      true_contrast <- custom_funs$fA(Q) - custom_funs$fA0(Q0)
      beta <- coef(glm.fit(V_full, true_contrast, family = gaussian(), weights = self$weights, intercept = FALSE))
      contrast <- as.vector(V_full %*% beta)



      IC <- as.matrix(EIF)

      result <- list(psi = beta, IC = IC, IC_list = IC_list, transform = self$custom_functions$transform)
      return(result)
    }
  ),
  active = list(
    name = function() {
      param_form <- c(
        paste0("1(A>0)*", private$.formula_names$bin),
        paste0("A*", private$.formula_names$cont)
      )
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
    formula_CATE_continuous = function() {
      return(private$.formula_CATE_continuous)
    },
    formula_CATE_binary = function() {
      return(private$.formula_CATE_binary)
    },
    custom_functions = function(){
      return(private$.custom_functions)
    }
  ),
  private = list(
    .type = "contCATE",
    .cf_likelihood = NULL,
    .supports_outcome_censoring = TRUE,
    .formula_CATE_binary = NULL,
    .formula_CATE_continuous = NULL,
    .formula_names = NULL,
    .custom_functions = NULL
  )
)
