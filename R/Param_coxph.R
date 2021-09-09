#' Nonparametric inference for user-specified parametric working models for the conditional treatment effect.
#' The true conditional average treatment effect is projected onto a parametric working model using least-squares regression.
#' Unlike \code{Param_npCATT}, this function uses all observations to compute the projection.
#' This can be used to assess heterogeneity of the average treatment effect.
#' We note that `formula_coxph = ~ 1` gives an estimator of the nonparametric average treatment effect (ATE).
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
#'     \item{\code{observed_likelihood}}{A \code{\link{Likelihood}} cocoxphesponding to the observed likelihood
#'     }
#'     \item{\code{formula_coxph}}{...
#'     }
#'     \item{\code{intervention_list_treatment}}{A list of objects inheriting from \code{\link{LF_base}}, representing the treatment intervention.
#'     }
#'     \item{\code{intervention_list_control}}{A list of objects inheriting from \code{\link{LF_base}}, representing the control intervention.
#'     }
#'     \item{\code{...}}{Not cucoxphently used.
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
Param_coxph <- R6Class(
  classname = "Param_coxph",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, formula_coxph = ~1, intervention_list_treatment, intervention_list_control, family_fluctuation = c("binomial"), outcome_node = "N") {

      super$initialize(observed_likelihood, list(), outcome_node = outcome_node)
      family_fluctuation <- match.arg(family_fluctuation)
      training_task <- self$observed_likelihood$training_task
      W <- training_task$get_regression_task("W", is_time_variant = TRUE)$Y


      V <- model.matrix(formula_coxph, as.data.frame(W))
      private$.formula_names <- colnames(V)
      private$.targeted <- rep(T, ncol(V))
      private$.submodel <- list(N = family_fluctuation)


      if (!is.null(observed_likelihood$censoring_nodes[[outcome_node]])) {
        # add delta_Y=0 to intervention lists
        outcome_censoring_node <- observed_likelihood$censoring_nodes[[outcome_node]]
        censoring_intervention <- define_lf(LF_static, outcome_censoring_node, value = 1)
        intervention_list_treatment <- c(intervention_list_treatment, censoring_intervention)
        intervention_list_control <- c(intervention_list_control, censoring_intervention)
      }
      private$.formula_coxph <- formula_coxph
      private$.cf_likelihood_treatment <- CF_Likelihood$new(observed_likelihood, intervention_list_treatment)
      private$.cf_likelihood_control <- CF_Likelihood$new(observed_likelihood, intervention_list_control)
    },
    long_to_mat = function(x, id, time) {
      dt <- data.table(id = id, time = time, x = as.vector(x))
      wide <- dcast(dt, id ~ time, value.var = "x")
      mat <- as.matrix(wide[, -1, with = FALSE])
      return(mat)
    },
    hm_to_sm = function(hm) {
      # TODO: check
      sm <- t(apply(1 - hm, 1, cumprod))
      # sm <- cbind(1,sm[,-ncol(sm)])
      return(sm)
    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full", is_training_task = TRUE) {
      training_task <- self$observed_likelihood$training_task
      if (is.null(tmle_task)) {
        tmle_task <- training_task
      }


      cf_task1 <- self$cf_likelihood_treatment$enumerate_cf_tasks(tmle_task)[[1]]
      cf_task0 <- self$cf_likelihood_control$enumerate_cf_tasks(tmle_task)[[1]]
      intervention_nodes <- union(names(self$intervention_list_treatment), names(self$intervention_list_control))

      W <- tmle_task$get_regression_task("W", is_time_variant = TRUE)$Y
      A <- tmle_task$get_tmle_node("A", format = T)[[1]]
      dNt <- tmle_task$get_tmle_node("N", format = F)
      dCt <- tmle_task$get_tmle_node("A_c", format = F)
      prefailure <- as.numeric(tmle_task$get_tmle_node("pre_failure"))
      Vt <- model.matrix(self$formula_coxph, as.data.frame(W))

      g <- self$observed_likelihood$get_likelihoods(tmle_task, "A", fold_number)
      g1 <- ifelse(A == 1, g, 1 - g)
      g0 <- 1 - g1

      pN <- self$observed_likelihood$get_likelihoods(tmle_task, "N", fold_number)
      pC <- self$observed_likelihood$get_likelihoods(tmle_task, "A_c", fold_number)
      pN0 <- as.vector(self$observed_likelihood$get_likelihoods(cf_task0, "N", fold_number))
      pN1 <- self$observed_likelihood$get_likelihoods(cf_task1, "N", fold_number)

      time <- tmle_task$time
      id <- tmle_task$id
      long_order <- order(id, time)

      pC_mat <- self$long_to_mat(pC, id, time)
      S_censor_mat <- self$hm_to_sm(pC_mat)
      S_censor_mat <- cbind(1, S_censor_mat[, -ncol(S_censor_mat)])
      S_censor <-  pmax(as.vector(S_censor_mat), 0.005)# Back to long, CHECK
      pN_mat <- self$long_to_mat(pN, id, time)
      S_surv_mat <- self$hm_to_sm(pN_mat)
      S_surv_mat <- cbind(1, S_surv_mat[, -ncol(S_surv_mat)])

      S_surv <- pmax(as.vector(S_surv_mat), 0.005)

      beta <- suppressWarnings(coef(glm.fit(Vt, pN1, offset = log(pN0), family = poisson(), weights = self$weights)))
      HR <- as.vector(exp(Vt %*% beta))

      t_grid <- sort(unique(time))


      H <- as.matrix(Vt * (prefailure / S_censor / S_surv) * (A / g1 * HR - (1 - A) / g0))

      #print(quantile(H))

      EIF_N <- NULL

      # Store EIF component
      if (is_training_task) {
        scale <- apply(Vt, 2, function(v) {
          apply(self$weights * Vt * (v) * HR * pN0, 2, sum) / length(unique(id))
        })



        scaleinv <- solve(scale)
        EIF_N <- self$weights * (H %*% scaleinv) * as.vector(dNt - pN)
        EIF_WA <- apply(Vt, 2, function(v) {
          long_vec <- self$weights * (v * (HR * pN0 - pN1))
          wide_vec <- self$long_to_mat(long_vec, id, time)
          means <- colMeans(wide_vec)
          as.vector(t(t(wide_vec) - means))
        }) %*% scaleinv



      }





      return(list(N = H, EIF = list(N = EIF_N, WA = EIF_WA)))
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      cf_task1 <- self$cf_likelihood_treatment$enumerate_cf_tasks(tmle_task)[[1]]
      cf_task0 <- self$cf_likelihood_control$enumerate_cf_tasks(tmle_task)[[1]]

      W <- tmle_task$get_regression_task("W", is_time_variant = TRUE)$Y
      A <- tmle_task$get_tmle_node("A", format = T)[[1]]
      dNt <- tmle_task$get_tmle_node("N", format = F)
      dCt <- tmle_task$get_tmle_node("A_c", format = F)
      prefailure <- as.numeric(tmle_task$get_tmle_node("pre_failure"))
      Vt <- model.matrix(self$formula_coxph, as.data.frame(W))

      weights <- tmle_task$weights
      time <- tmle_task$time
      id <- tmle_task$id
      long_order <- order(id, time)
      # clever_covariates happen here (for this param) only, but this is repeated computation
      EIF <- self$clever_covariates(tmle_task, fold_number, is_training_task = TRUE)$EIF
      EIF <- EIF$N + EIF$WA

      EIF <- apply(EIF, 2, function(col) {
        rowSums(self$long_to_mat(col, id, time))
      })

      pN <- self$observed_likelihood$get_likelihoods(tmle_task, "N", fold_number)
      pC <- self$observed_likelihood$get_likelihoods(tmle_task, "A_c", fold_number)
      pN0 <- self$observed_likelihood$get_likelihoods(cf_task0, "N", fold_number)
      pN1 <- self$observed_likelihood$get_likelihoods(cf_task1, "N", fold_number)



      pC_mat <- self$long_to_mat(pC, id, time)
      S_censor_mat <- self$hm_to_sm(pC_mat)
      S_censor_mat <- cbind(1, S_censor_mat[, -ncol(S_censor_mat)])
      S_censor <- as.vector(S_censor_mat) # Back to long, CHECK



      beta <- suppressWarnings(coef(glm.fit(Vt, pN1, offset = log(pN0), family = poisson(), weights =  self$weights )))


      HR <- exp(Vt %*% beta)

      IC <- as.matrix(EIF)

      result <- list(psi = beta, IC = IC, HR = HR, transform = exp)
      return(result)
    }
  ),
  active = list(
    name = function() {
      param_form <- private$.formula_names # sprintf("coxph[%s_{%s}-%s_{%s}]", self$outcome_node, self$cf_likelihood_treatment$name, self$outcome_node, self$cf_likelihood_control$name)
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
    formula_coxph = function() {
      return(private$.formula_coxph)
    }
  ),
  private = list(
    .type = "coxph (HR)",
    .cf_likelihood_treatment = NULL,
    .cf_likelihood_control = NULL,
    .supports_outcome_censoring = TRUE,
    .formula_coxph = NULL,
    .submodel = list(N = "binomial_logit"),
    .formula_names = NULL
  )
)
