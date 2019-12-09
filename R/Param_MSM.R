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
    initialize = function(observed_likelihood, strata_variable, 
                          continuous_treatment = FALSE, treatment_values = NULL, n_samples = 30,
                          mass = NULL, ...,
                          covariate_node = "W", treatment_node = "A", outcome_node = "Y") {
      super$initialize(observed_likelihood, ..., outcome_node = "Y")
      private$.covariate_node <- covariate_node
      private$.treatment_node <- treatment_node
      private$.continuous_treatment <- continuous_treatment
      if (continuous_treatment) {
        private$.n_samples <- n_samples
        private$.treatment_values <- "A" # not needed but stored for simplicity
      } else {
        # numeralize treatments and store mapping
        private$.treatment_values <- setNames(paste0("A", 1:length(treatment_values)), treatment_values)
      }
      private$.strata_variable <- strata_variable
      '
      V <- observed_likelihood$training_task$get_data(, strata_variable)
      strata <- V[, list(weight = observed_likelihood$training_task$nrow / .N), by = names(V)]
      set(strata, , "strata_i", factor(1:nrow(strata)))
      private$.strata <- strata
      '
      if (is.null(private$.mass)) {
        private$.mass <- function(tmle_task = NULL, fold_number = "full") {
          self$observed_likelihood$get_likelihoods(tmle_task, self$treatment_node, fold_number)
        }
      } else {
        private$.mass <- function(tmle_task = NULL, fold_number = "full") {
          # TODO: fold?
          A <- tmle_task$get_tmle_node(self$treatment_node)
          V <- tmle_task$get_tmle_node(self$covariate_node)[[self$strata_variable]]
          sapply(1:length(A), function(i) mass(A[i], V[i]))
        }
      }
    },
    
    clever_covariates = function(tmle_task = NULL, fold_number = "full") {
      pA <- self$observed_likelihood$get_likelihoods(tmle_task, self$treatment_node, fold_number)
      h <- private$.mass(tmle_task, fold_number)
      
      base_covs <- 1. / pA
      strata_covs <- c(h * base_covs)
      
      A <- as.matrix(tmle_task$get_tmle_node(self$treatment_node))
      V <- as.matrix(tmle_task$get_data(, self$strata_variable))
      if (self$continuous_treatment) {
        clever_covs <- strata_covs * cbind(A, V)
      } else {
        covs_A_mat <- matrix(, nrow = length(h), ncol = length(self$treatment_values))
        for (i in 1:length(h)) {
          covs_A_mat[i, which(self$treatment_values == A[i,])] <- strata_covs[i]
        }
        clever_covs <- cbind(covs_A_mat, strata_covs * V)
      }
      
      return(list(Y = clever_covs))
    },
    
    estimates = function(tmle_task = NULL, fold_number = "full") {
      # TODO: estimate variance of MC integration
      Y <- as.matrix(tmle_task$get_tmle_node(self$outcome_node))
      qY <- self$observed_likelihood$get_likelihoods(tmle_task, self$outcome_node, fold_number)
      V <- as.vector(as.matrix(tmle_task$get_data(, self$strata_variable)))
      A <- as.matrix(tmle_task$get_tmle_node(self$treatment_node))
      pA <- self$observed_likelihood$get_likelihoods(tmle_task, self$treatment_node, fold_number)
      h <- private$.mass(tmle_task, fold_number)
      if (self$continuous_treatment) {
        A_range <- self$observed_likelihood$factor_list[[self$treatment_node]]$get_outcome_range(tmle_task, fold_number)
        U <- runif(100)
        A_vals <- t(apply(A_range, 1, function(r) r[1] + (r[2] - r[1]) * U))
        # Generate counterfactual tasks for each sample of A:
        cf_tasks <- lapply(A_samp, function(A_val) {
          newdata <- data.table(A = A_val)
          cf_task <- tmle_task$generate_counterfactual_task(UUIDgenerate(), new_data = newdata)
          return(cf_task)
        })
      } else {
        A_vals <- diag(length(self$treatment_values))
        # Generate counterfactual tasks for each value of A:
        cf_tasks <- lapply(names(self$treatment_values), function(A_val) {
          newdata <- data.table(A = A_val)
          cf_task <- tmle_task$generate_counterfactual_task(UUIDgenerate(), new_data = newdata)
          return(cf_task)
        })
      }
      
      qY_vals <- sapply(cf_tasks, self$observed_likelihood$get_likelihood, "Y", fold_number)
      h_vals <- sapply(cf_tasks, private$.mass, fold_number)
      
      #
      qY_ext <- c(qY_vals)
      if (self$continuous_treatment) {
        V_ext <- rep(V, times = self$n_samples)
        A_ext <- c(A_vals)
      } else {
        V_ext <- rep(V, times = length(self$treatment_values))
        A_ext <- sapply(data.frame(A_vals), rep, each=length(V))
      }
      h_ext <- c(h_vals)
      regress_mat = cbind(qY_ext, A_ext, V_ext)
      
      regress_table <- data.table(regress_mat)
      colnames(regress_table) <- c("Q", self$treatment_values, "V")
      formula <- as.formula(paste0("Q~", 
                                   paste(c(self$treatment_values, "V"), collapse = "+"),
                                   "-1"))
      if (length(levels(factor(Y))) == 2 && all.equal(levels(factor(Y)), c("0", "1"))) {
        model <- glm(formula, family = binomial, data = regress_table, weights = h_ext)
      } else {
        model <- glm(formula, family = gaussian, data = regress_table, weights = h_ext)
      }
      
      psi <- model$coefficients
      
      weighted_res <- residuals(model, type = "response") * h_ext
      if (self$continuous_treatment) {
        res_A_mat <- rowMeans(matrix(weighted_res * A_ext, nrow = length(V), byrow = FALSE))
        res_V_mat <- rowMeans(matrix(weighted_res, nrow = length(V), byrow = FALSE)) * V
      } else {
        res_A_mat <- matrix(weighted_res, nrow = length(V), byrow = FALSE)
        res_V_mat <- rowSums(res_A_mat) * V
      }
      IC_term2 <- cbind(res_A_mat, res_V_mat)
      
      base_res <- (Y - qY) / pA
      strata_res <- c(h * base_res)
      if (self$continuous_treatment) {
        IC_term1 <- strata_res * cbind(A, V)
      } else {
        ic_A_mat <- matrix(, nrow = length(h), ncol = length(self$treatment_values))
        for (i in 1:length(h)) {
          ic_A_mat[i, which(self$treatment_values == A[i,])] <- strata_res[i]
        }
        IC_term1 <- cbind(ic_A_mat, strata_res * V)
      }
      
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
      # deprecated
      treatment_variable <- tmle_task$npsem[["A"]]$variables
      A <- tmle_task$get_data(, treatment_variable)
      treatment <- factor(A[[treatment_variable]])
      levels(treatment) <- self$treatment_values[levels(treatment)]
      
      combined <- A[, "indicators":= 1][, "level":= treatment]
      treatment_indicators <- t(sapply(combined$level, function(a) as.numeric(a == self$treatment_values)))
      colnames(treatment_indicators) = self$treatment_values
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
    continuous_treatment = function() {
      return(private$.continuous_treatment)
    },
    treatment_values = function() {
      return(private$.treatment_values)
    },
    n_samples = function() {
      return(private$.n_samples)
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
    covariate_node = function() {
      return(private$.covariate_node)
    },
    treatment_node = function() {
      return(private$.treatment_node)
    }
  ),
  
  private = list(
    .type = "stratified MSM",
    .continuous_treatment = NULL,
    .treatment_values = NULL,
    .n_samples = NULL,
    .strata_variable = NULL,
    .strata = NULL,
    .mass = NULL,
    .covariate_node = NULL,
    .treatment_node = NULL
  )
)
