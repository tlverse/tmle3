


context("Coxph does not error.")



if (require(simcausal)) {
  tmax <- 5
  print(i)
  D <- DAG.empty()
  D <- D + node("W1", distr = "runif", min = -1, max = 1) +
    node("Wa2", distr = "rbinom", size = 1, prob = 0.5) +
    node("W2", distr = "rconst", const = Wa2 - 0.5) +
    node("A", distr = "rbinom", size = 1, prob = plogis(W1 + W2)) +
    node("dNt", t = 1:tmax, EFU = TRUE, distr = "rbinom", size = 1, prob = exp(0.5 * A) * 0.15 * plogis(W1 + W2)) +
    node("dCt", t = 1:tmax, EFU = TRUE, distr = "rbinom", size = 1, prob = 0 * plogis(W1 + W2 + t))
  D <- set.DAG(D)
  data <- sim(D, n = 1000)
  data

  data_N <- data[, grep("[d][N].+", colnames(data)), drop = F]
  data_C <- data[, grep("[d][C].+", colnames(data)), drop = F]

  data_surv <- as.data.frame(do.call(rbind, lapply(1:nrow(data), function(i) {
    rowN <- data_N[i, ]
    rowC <- data_C[i, ]
    t <- which(rowN == 1)
    tc <- which(rowC == 1)
    if (length(t) == 0) {
      t <- tmax + 2
    }
    if (length(tc) == 0) {
      tc <- tmax + 1
    }

    Ttilde <- min(t, tc)
    Delta <- t <= tc
    return(matrix(c(Ttilde, Delta), nrow = 1))
  })))
  colnames(data_surv) <- c("Ttilde", "Delta")
  data$Ttilde <- data_surv$Ttilde
  data$Delta <- data_surv$Delta
  data <- data[, -grep("[d][C].+", colnames(data))]
  data <- data[, -grep("[d][N].+", colnames(data))]
  data



  doMC::registerDoMC(16)

  tmle_spec_np <- tmle3_Spec_coxph$new(formula = ~1, delta_epsilon = 0.05, verbose = T, treatment_level = 1, control_level = 0)
  learner_list <- list(A = Lrnr_glm$new(), N = Lrnr_glm$new(formula = ~ .^2), A_c = Lrnr_glm$new(formula = ~ .^2))
  node_list <- list(W = c("W1", "W2"), A = "A", T_tilde = "Ttilde", Delta = "Delta")

  tmle3_fit <- suppressMessages(suppressWarnings(tmle3(tmle_spec_np, data, node_list, learner_list)))
}
