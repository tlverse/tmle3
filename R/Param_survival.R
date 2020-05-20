#' Survival Curve
#'
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
#'   \code{define_param(Param_survival, observed_likelihood, intervention_list, ..., outcome_node)}
#'
#'   \describe{
#'     \item{\code{observed_likelihood}}{A \code{\link{Likelihood}} corresponding to the observed likelihood
#'     }
#'     \item{\code{intervention_list}}{A list of objects inheriting from \code{\link{LF_base}}, representing the intervention.
#'     }
#'     \item{\code{...}}{Not currently used.
#'     }
#'     \item{\code{outcome_node}}{character, the name of the node that should be treated as the outcome
#'     }
#'     }
#'

#' @section Fields:
#' \describe{
#'     \item{\code{cf_likelihood}}{the counterfactual likelihood for this treatment
#'     }
#'     \item{\code{intervention_list}}{A list of objects inheriting from \code{\link{LF_base}}, representing the intervention
#'     }
#' }
#' @export
Param_survival <- R6Class(
  classname = "Param_survival",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, intervention_list, ..., outcome_node) {
      super$initialize(observed_likelihood, ..., outcome_node = outcome_node)
      private$.cf_likelihood <- make_CF_Likelihood(observed_likelihood, intervention_list)
    },
    hazards_to_survival = function(p_hazards, t_max) {
      # TODO: whether change the input data format to consecutive times for each observation
      n <- length(p_hazards) / t_max
      p_surv <- copy(p_hazards)
      for (i in 1:n) {
        p_surv[i + n * seq(0, t_max - 1)] <- cumprod(1 - p_hazards[i + n * seq(0, t_max - 1)])
      }
      return(p_surv)
    },
    get_pSN_at_time = function(pS_N, time, t_max, long_format = TRUE) {
      n <- length(pS_N) / t_max
      pS_N_time <- pS_N[seq(1 + (time - 1) * n, time * n)]
      if (long_format) {
        pS_N_time <- rep(pS_N_time, t_max)
      } 
      return(pS_N_time)
    },
    get_single_time_ht = function(time, pA, cf_pA, pS_N, pS_A_c, k_list, t_max) {
      I_k <- ifelse(k_list <= time, 1, 0)
      pS_N_time <- self$get_pSN_at_time(pS_N, time, t_max)
      ht <- -((cf_pA * I_k) / (pA * pS_A_c)) * (pS_N_time / pS_N)
      return(ht)
    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      intervention_nodes <- names(self$intervention_list)
      pA <- self$observed_likelihood$get_likelihoods(tmle_task, intervention_nodes, fold_number)
      # I(A=1)
      cf_pA <- self$cf_likelihood$get_likelihoods(tmle_task, intervention_nodes, fold_number)

      # TODO: whether modify LF_fit_hazards get_density without 1 - preds
      pN <- self$observed_likelihood$get_likelihoods(tmle_task, "N", fold_number)
      pA_c <- self$observed_likelihood$get_likelihoods(tmle_task, "A_c", fold_number)

      t_max <- max(tmle_task$get_tmle_node("T_tilde"))
      pS_N <- self$hazards_to_survival(pN, t_max)
      pS_A_c <- self$hazards_to_survival(pA_c, t_max)

      k_list <- tmle_task$get_tmle_node("t")
      all_ht <- lapply(seq(t_max), function(time) {
        self$get_single_time_ht(time, pA, cf_pA, pS_N, pS_A_c, k_list, t_max)
      })
      all_ht_dt <- as.data.table(all_ht)

      # TODO: return format
      HA <- all_ht_dt
    },
    get_psi = function(pS_N1, t_max) {
      n <- length(pS_N1) / t_max
      psi <- lapply(seq(t_max), function(t) {
        mean(pS_N1[seq(1 + (t - 1) * n, t * n)])
      })
      return(unlist(psi))
    },
    get_single_time_Dt = function(time, HA, pN1, pS_N1, psi, I1_mat, I2_mat, t_max) {
      n <- length(pS_N1) / t_max
      # TODO: initialize properly
      Dt <- NULL
      for (k in 1:time) {
        temp = HA[seq(1 + (k - 1) * n, k * n), time] * 
        (I1_mat[, k] - I2_mat[, k] * pN1[seq(1 + (k - 1) * n, k * n)])
        if (is.null(Dt)) {
          Dt = temp
        } else {
          Dt = Dt + temp
        }
      }
      Dt = Dt + self$get_pSN_at_time(pS_N1, time, t_max, long_format = FALSE) - psi[time]
      return(Dt)
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }

      cf_task <- self$cf_likelihood$enumerate_cf_tasks(tmle_task)[[1]]

      # TODO: return format
      HA <- self$clever_covariates(tmle_task, fold_number)

      t_max <- max(tmle_task$get_tmle_node("T_tilde"))
      pN1 <- self$observed_likelihood$get_likelihoods(cf_task, "N", fold_number)
      pS_N1 <- self$hazards_to_survival(pN1, t_max)
      psi <- self$get_psi(pS_N1, t_max)

      # TODO: create I1 and I2
    }
  ),
  active = list(
    name = function() {
      param_form <- sprintf("E[%s_{%s}]", self$outcome_node, self$cf_likelihood$name)
      return(param_form)
    },
    cf_likelihood = function() {
      return(private$.cf_likelihood)
    },
    intervention_list = function() {
      return(self$cf_likelihood$intervention_list)
    },
    update_nodes = function() {
      return(self$outcome_node)
    }
  ),
  private = list(
    .type = "survival",
    .cf_likelihood = NULL
  )
)
