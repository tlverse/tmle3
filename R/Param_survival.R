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
      # TODO: check outcome_node, current I(T<=t, delta=1), need I(T=t, delta=1)
      super$initialize(observed_likelihood, ..., outcome_node = outcome_node)
      private$.cf_likelihood <- make_CF_Likelihood(observed_likelihood, intervention_list)
    },
    reshape_long_data = function(long_data, t_max) {
      n <- length(long_data) / t_max
      # TODO: assume long_data is a list
      rs <- list()
      for (i in 1:t_max) {
        current <- long_data[seq(1 + (i - 1) * n, i * n)]
        rs <- c(rs, list(current))
      }
      rs <- do.call(cbind, rs)
      return(rs)
    },
    hazards_to_survival = function(p_hazards, t_max) {
      # TODO: whether change the input data format to consecutive times for each observation
      n <- length(p_hazards) / t_max
      p_surv <- copy(p_hazards)
      for (i in 1:n) {
        # TODO: make hazard at starting time 0
        temp <- p_hazards[i + n * seq(0, t_max - 1)]
        temp <- c(0, temp)
        p_surv[i + n * seq(0, t_max - 1)] <- cumprod(1 - temp[-length(temp)])
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
      # TODO: check
      HA <- as.matrix(HA)
      return(list(N = HA))
    },
    get_psi = function(pS_N1, t_max) {
      n <- length(pS_N1) / t_max
      psi <- lapply(seq(t_max), function(t) {
        mean(pS_N1[seq(1 + (t - 1) * n, t * n)])
      })
      return(unlist(psi))
    },
    get_single_time_Dt = function(time, HA, pN1, pS_N1, psi, T_tilde_data_short, Delta_data_short, t_max) {
      n <- length(pS_N1) / t_max
      # TODO: initialize properly
      Dt <- NULL
      for (k in 1:time) {
        # TODO: optimize I1, I2 creation
        I1 <- ifelse(T_tilde_data_short == k & Delta_data_short == 1, 1, 0)
        I2 <- ifelse(T_tilde_data_short >= k, 1, 0)
        # TODO: check
        temp = HA[seq(1 + (k - 1) * n, k * n), time] * 
        (I1 - I2 * pN1[seq(1 + (k - 1) * n, k * n)])
        if (is.null(Dt)) {
          Dt = temp
        } else {
          Dt = Dt + temp
        }
      }
      Dt = Dt + self$get_pSN_at_time(pS_N1, time, t_max, long_format = FALSE) - psi[time]
      return(Dt)
    },
    # get_all_Dt = function(HA, pN1, pS_N1, psi, T_tilde_data_short, Delta_data_short, t_max) {
    #   n <- length(pS_N1) / t_max

    #   all_Dt_mat <- matrix(rep(0, n * t_max), nrow = n, ncol = t_max, byrow = TRUE)
    #   all_Dt <- as.data.frame(all_Dt_mat)
    #   cum_sum1 <- NULL
    #   for (t in 1:t_max) { 
    #     I1 <- ifelse(T_tilde_data_short == t & Delta_data_short == 1, 1, 0)
    #     I2 <- ifelse(T_tilde_data_short >= t, 1, 0)
    #     # TODO: not work because ht differs for different t
    #     part1 <- HA[seq(1 + (t - 1) * n, t * n), t] * 
    #     (I1 - I2 * pN1[seq(1 + (t - 1) * n, t * n)])
    #     part2 <- self$get_pSN_at_time(pS_N1, t, t_max, long_format = FALSE) - psi[t]
    #     if (is.null(cum_sum1)) {
    #       cum_sum1 <- part1
    #     } else {
    #       cum_sum1 = cum_sum1 + part1
    #     }
    #     all_Dt[, t] <- cum_sum1 + part2
    #   }
    #   return(all_Dt)
    # },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }

      cf_task <- self$cf_likelihood$enumerate_cf_tasks(tmle_task)[[1]]

      # TODO: return format
      HA <- self$clever_covariates(tmle_task, fold_number)[["N"]]

      T_tilde_data <- tmle_task$get_tmle_node("T_tilde")
      Delta_data <- tmle_task$get_tmle_node("Delta")
      t_max <- max(T_tilde_data)
      n <- length(T_tilde_data) / t_max
      T_tilde_data_short <- T_tilde_data[seq(n)]
      Delta_data_short <- Delta_data[seq(n)]

      pN1 <- self$observed_likelihood$get_likelihoods(cf_task, "N", fold_number)
      pS_N1 <- self$hazards_to_survival(pN1, t_max)

      psi <- self$get_psi(pS_N1, t_max)

      # # TODO: create I1 and I2
      # I1_mat <- NULL
      # I2_mat <- NULL

      # TODO: as time gets large, slow
      all_Dt <- lapply(seq(t_max), function(time) {
        self$get_single_time_Dt(time, HA, pN1, pS_N1, psi, T_tilde_data_short, Delta_data_short, t_max)
      })
      all_Dt_table <- as.data.table(all_Dt)

      # TODO: return format
      IC <- all_Dt_table
      IC <- as.matrix(IC)
      result <- list(psi = psi, IC = IC)
      return(result)
    }
  ),
  active = list(
    # TODO: modify
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
