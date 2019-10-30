#' Stratified Parameter Estimates
#'
#'
#' @section Current Issues:
#' \itemize{
#'   \item clever covariates doesn't support updates; always uses initial (necessary for iterative TMLE, e.g. stochastic intervention)
#'   \item doesn't integrate over possible counterfactuals (necessary for stochastic intervention)
#'   \item clever covariate gets recalculated all the time (inefficient)
#' }
#' @importFrom R6 R6Class
#' @importFrom stats glm predict
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
#'   \code{define_param(Param_TSM, observed_likelihood, intervention_list, ..., outcome_node)}
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
Param_MSM <- R6Class(
  classname = "Param_MSM",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, treatment_values, strata_variable,
                          mass = NULL, ..., treatment_node = "A", outcome_node = "Y") {
      super$initialize(observed_likelihood, ..., outcome_node = "Y")
      private$.treatment_node <- treatment_node
      private$.type <- sprintf("stratified %s", param_base$type)
      # numeralize treatments and store mapping
      private$.treatment_values <- setNames(1:length(treatment_values), treatment_values)
      
      private$.strata_variable <- strata_variable
      
      V <- observed_likelihood$training_task$get_data(, strata_variable)
      strata <- V[, list(weight = observed_likelihood$training_task$nrow / .N), by = names(V)]
      set(strata, , "strata_i", factor(1:nrow(strata)))
      private$.strata <- strata
      # TODO: implement the case where mass is given
      # private$.mass <- mass
    },
    
    clever_covariates = function(tmle_task = NULL, fold_number = "full") {
      pA <- self$observed_likelihood$get_likelihoods(tmle_task, self$treatment_node, fold_number)
      base_covs <- 1. / pA
      if (is.null(private$.mass)) {
        strata_weights <- self$get_strata_weights(tmle_task)
      } else {
        # TODO: implement the case where mass is given
      }
      strata_covs <- base_covs * strata_weights
      
      treatment_indicators <- self$get_treatment_indicators(tmle_task)
      V <- as.matrix(tmle_task$get_data(, self$strata_variable))
      covs_matrix <- cbind(treatment_indicators, V)
      clever_covs <- strata_covs * covs_matrix
      
      return(clever_covs)
    },
    
    estimates = function(tmle_task = NULL, fold_number = "full") {
      qY <- self$observed_likelihood$get_likelihoods(tmle_task, self$outcome_node, fold_number)
      V <- as.vector(as.matrix(tmle_task$get_data(, self$strata_variable)))
      A <- diag(length(self$treatment_values))
      
      qY_ext <- rep(qY, times = length(self$treatment_values))
      V_ext <- rep(V, times = length(self$treatment_values))
      A_ext <- sapply(data.frame(A), rep, each=length(V))
      regress_mat = cbind(qY_ext, V_ext, A_ext)
      
      regress_table <- data.table(regress_mat)
      colnames(regress_table) <- c("Q", self$strata_variable, names(self$treatment_values))
      formula <- as.formula(paste0(c(
        "Q~",
        paste(c(self$strata_variable, names(self$treatment_values)), collapse = "+"),
        "-1"
      )))
      if (is.null(private$.mass)) {
        model <- glm(formula, family = binomial(link = "logit"), data = regress_table)
      } else {
        # TODO: implement the case where mass is given
      }
      psi <- model$coefficients
      
      residuals <- residuals(model, type = "response")
      residuals_A_mat <- matrix(residuals, nrow = length(V), byrow = FALSE)
      residuals_V_mat <- rowSums(residuals_A_mat) * V
      IC_term2 <- cbind(residuals_A_mat, residuals_V_mat)
      
      pA <- self$observed_likelihood$get_likelihoods(tmle_task, self$treatment_node, fold_number)
      Y <- tmle_task$get_tmle_node(self$outcome_node)
      base_res <- (Y - qY) / pA
      if (is.null(private$.mass)) {
        strata_weights <- self$get_strata_weights(tmle_task)
      } else {
        # TODO: implement the case where mass is given
      }
      strata_res <- base_res * strata_weights
      treatment_indicators <- self$get_treatment_indicators(tmle_task)
      covs_matrix <- cbind(treatment_indicators, V)
      IC_term1 <- strata_res * covs_matrix
      
      IC <- IC_term1 + IC_term2
      
      result <- list(psi = psi, IC = IC)
      return(result)
    },
    
    get_strata_weights = function(tmle_task) {
      V <- tmle_task$get_data(, self$strata_variable)
      strata <- self$strata
      combined <- merge(V, strata, by = self$strata_variable, sort = FALSE, all.x = TRUE)
      strata_weights <- combined$weight
      return(strata_weights)
    },
    
    get_treatment_indicators = function(tmle_task) {
      treatment_variable <- tmle_task$npsem[["A"]]$variables
      A <- tmle_task$get_data(, treatment_variable)
      treatment = factor(A[[treatment_variable]])
      levels(treatment) <- self$treatment_values[levels(treatment)]
      
      combined <- A[, "indicators":= 1][, "level":= treatment]
      combined[, index := .I]
      treatment_indicators_dt <- dcast(combined, index ~ level, value.var = "indicators", fill = 0, drop = FALSE)
      treatment_indicators <- as.matrix(strata_indicators_dt[, -1, with = FALSE])
      return(treatment_indicators)
    }
    
  ),
  
  active = list(
    name = function() {
      treatment_labels <- names(self$treatment_values)
      strata_label <- self$strata_variable
      param_form <- sprintf(
        "Param_MSM: %s",
        c(treatment_labels, strata_label)
      )
      return(param_form)
    },
    treatment_values = function() {
      return(private$.treatment_values)
    },
    strata_variable = function() {
      return(private$.strata_variable)
    },
    update_nodes = function() {
      return(self$outcome_node)
    },
    strata = function() {
      return(private$.strata)
    },
    treatment_node = function() {
      return(private$.treatment_node)
    }
  ),
  
  private = list(
    .type = NULL,
    .treatment_values = NULL,
    .strata_variable = NULL,
    .strata = NULL,
    .mass = NULL,
    .treatment_node = NULL
  )
)
