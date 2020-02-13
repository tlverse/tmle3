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
                          mass = NULL, mass_ub = 1/0.025, ...,
                          covariate_node = "W", treatment_node = "A", outcome_node = "Y") {
      super$initialize(observed_likelihood, ..., outcome_node = "Y")
      private$.covariate_node <- covariate_node
      private$.treatment_node <- treatment_node
      private$.continuous_treatment <- continuous_treatment
      private$.mass_ub <- mass_ub
      if (continuous_treatment) {
        private$.n_samples <- n_samples
        private$.U <- runif(n_samples)
        private$.treatment_values <- setNames(treatment_node, treatment_node) # not needed but stored for simplicity
      } else {
        if (is.null(treatment_values)) {
          treatment_values <- observed_likelihood$factor_list[[treatment_node]]$variable_type$levels
        }
        # numeralize treatments and store mapping
        private$.treatment_values <- setNames(treatment_values, paste0(treatment_node, "_", treatment_values))
        # 
        private$.cf_likelihoods <- lapply(self$treatment_values, function(A_level) {
          intervention_list <- define_lf(LF_static, "A", value = A_level)
          cf_likelihood <- make_CF_Likelihood(observed_likelihood, intervention_list)
          return(cf_likelihood)
        })
        names(private$.cf_likelihoods) <- names(self$treatment_values)
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
          m = self$observed_likelihood$get_likelihoods(tmle_task, self$treatment_node, fold_number)
          if (!is.null(self$mass_ub)) {
            m = pmin(m, self$mass_ub)
          }
          m
        }
      } else {
        private$.mass <- function(tmle_task = NULL, fold_number = "full") {
          # TODO: fold?
          A <- tmle_task$get_tmle_node(self$treatment_node)
          V <- tmle_task$get_tmle_node(self$covariate_node)[[self$strata_variable]]
          m = sapply(1:length(A), function(i) mass(A[i], V[i]))
          if (!is.null(self$mass_ub)) {
            m = pmin(m, self$mass_ub)
          }
          m
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
        A_mat <- self$get_treatment_indicators(tmle_task)
        clever_covs <- strata_covs * cbind(A_mat, V)
      }
      colnames(clever_covs) <- c(names(self$treatment_values), "V")
      
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
        A_range <- self$observed_likelihood$factor_list[[self$treatment_node]]$learner$get_outcome_range(tmle_task$get_regression_task(self$outcome_node), fold_number)
        A_vals <- t(apply(A_range, 1, function(r) r[1] + (r[2] - r[1]) * self$U))
        # Generate counterfactual tasks for each sample of A:
        cf_tasks <- apply(A_vals, 2, function(A_val) {
          newdata <- data.table(A = A_val)
          cf_task <- tmle_task$generate_counterfactual_task(UUIDgenerate(), new_data = newdata)
          return(cf_task)
        })
      } else {
        A_vals <- diag(length(self$treatment_values))
        # Generate counterfactual tasks for each value of A:
        cf_tasks <- lapply(self$cf_likelihoods, function(cf_likelihood) {
          cf_task <- cf_likelihood$enumerate_cf_tasks(tmle_task)[[1]]
          return(cf_task)
        })
      }
      
      qYA <- sapply(cf_tasks, self$observed_likelihood$get_likelihood, self$outcome_node, fold_number)
      hA <- sapply(cf_tasks, private$.mass, fold_number)
      
      # psi and ic
      qY_ext <- matrix(qYA, ncol = 1, byrow = FALSE)
      if (self$continuous_treatment) {
        V_ext <- rep(V, times = self$n_samples)
        A_ext <- matrix(A_vals, ncol = 1, byrow = FALSE)
      } else {
        V_ext <- rep(V, times = length(self$treatment_values))
        A_ext <- sapply(data.frame(A_vals), rep, each=length(V))
      }
      h_ext <- matrix(hA, ncol = 1, byrow = FALSE)
      
      regress_table <- data.table(cbind(qY_ext, A_ext, V_ext))
      colnames(regress_table) <- c("Q", names(self$treatment_values), "V")
      formula <- as.formula(paste0("Q~", 
                                   paste(c(names(self$treatment_values), "V"), collapse = "+"),
                                   "-1"))
      if (length(levels(factor(Y))) == 2 && all.equal(levels(factor(Y)), c("0", "1"))) {
        family = "binomial"
      } else {
        family = "gaussian"
      }
      suppressWarnings(model <- glm(formula, family = family, data = regress_table, weights = h_ext))
      
      psi <- model$coefficients
      names(psi) <- c(names(self$treatment_values), "V")
      
      weighted_res <- residuals(model, type = "response") * h_ext
      if (self$continuous_treatment) {
        res_A <- (A_range[, 2] - A_range[, 1]) * rowMeans(matrix(weighted_res * A_ext, nrow = length(V), byrow = FALSE))
        res_V <- (A_range[, 2] - A_range[, 1]) * rowMeans(matrix(weighted_res, nrow = length(V), byrow = FALSE)) * V
      } else {
        res_A <- matrix(weighted_res, nrow = length(V), byrow = FALSE)
        res_V <- rowSums(res_A) * V
      }
      H2 <- cbind(res_A, res_V)
      
      H1 <- self$clever_covariates(tmle_task, fold_number)[[self$outcome_node]]
      
      IC <- H1 * c(Y - qY) + H2
      
      # normalization
      if (self$continuous_treatment) {
        # TODO: continuous normalization
      } else {
        ## reference: tmle - Susan Gruber
        if (!any(is.na(IC))) {
          nterms <- ncol(IC)
          ntreats <- length(self$treatment_values)
          f <- function(x) {
            x[1:nterms] %*% t(x[(nterms + 1):(2 * nterms)])
          }
          if (family == "binomial") {
            mA <- matrix(predict(model, type = "response"), nrow = length(V), byrow = FALSE)
            derivFactor <- mA * (1 - mA)
          } else {
            derivFactor <- matrix(1, nrow = nrow(IC), ncol = ntreats)
          }
          deriv.terms <- rep(list(matrix(NA, nterms^2, length(V))), ntreats)
          for (i in 1:length(self$treatment_values)) {
            covar.MSMA <- cbind(matrix(A_ext[, i], nrow = length(V), byrow = FALSE), V)
            deriv.terms[[i]] <- apply(cbind(-hA[, i] * derivFactor[, i] * covar.MSMA, covar.MSMA), 1, f)
          }
          ddpsi.IC <- as.matrix(Reduce('+', deriv.terms))
          M <- -matrix(rowMeans(ddpsi.IC), nrow = nterms)
          Minv <- try(solve(M))
          if (identical(class(Minv), "try-error")) {
            warning("Inference unavailable: normalizing matrix not invertible. IC not normalized. \n")
          }
          else {
            IC <- t(Minv %*% t(IC))
          }
        }
      }
      
      colnames(IC) <- c(names(self$treatment_values), "V")
      
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
    
    get_treatment_indicators = function(tmle_task, treatment_node = self$treatment_node) {
      treatment_variable <- tmle_task$npsem[[treatment_node]]$variables
      A <- tmle_task$get_data(, treatment_variable)
      treatment <- factor(A[[treatment_variable]])
      
      combined <- A[, "indicators":= 1][, "level":= treatment]
      treatment_indicators <- t(sapply(combined$level, function(a) as.numeric(a == self$treatment_values)))
      colnames(treatment_indicators) <- names(self$treatment_values)
      return(treatment_indicators)
    }
    
  ),
  
  active = list(
    name = function() {
      if (self$continuous_treatment) {
        treatment_labels <- self$treatment_node
        treatment_names <- "Treatment Node"
      } else {
        treatment_labels <- as.character(self$treatment_values)
        treatment_names <- rep("Treatment Level", length(treatment_labels))
      }
      strata_label <- self$strata_variable
      
      param_form <- sprintf(
        "%s: %s",
        c(treatment_names, "Strata Variable"), 
        c(treatment_labels, strata_label)
      )
      
      return(param_form)
    },
    continuous_treatment = function() {
      return(private$.continuous_treatment)
    },
    mass_ub = function() {
      return(private$.mass_ub)
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
    },
    U = function() {
      return(private$.U)
    },
    cf_likelihoods = function() {
      return(private$.cf_likelihoods)
    }
  ),
  
  private = list(
    .type = "stratified MSM",
    .continuous_treatment = NULL,
    .mass_ub = NULL,
    .treatment_values = NULL,
    .n_samples = NULL,
    .strata_variable = NULL,
    .strata = NULL,
    .mass = NULL,
    .covariate_node = NULL,
    .treatment_node = NULL,
    .U = NULL,
    .cf_likelihoods = NULL
  )
)