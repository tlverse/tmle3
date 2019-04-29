if (FALSE) {
  library(tmle3)
  library(future)
  library(future.apply)
  library(data.table)
  library(sl3)
  sl3_debug_mode()

  one_sim <- function(iteration = -1) {
    setDTthreads(1, restore_after_fork = FALSE)
    n <- 1000
    W <- rnorm(n)
    A <- rbinom(n, 1, plogis(W))
    A_cont <- plogis(rnorm(n, mean = plogis(W)))
    Y <- rbinom(n, 1, plogis(A + 0.2 * W))
    data <- data.table(W = W, A = factor(A), Y = Y, A_cont = A_cont)
    tmle_spec <- tmle_TSM_all()
    node_list <- list(W = "W", A = "A", Y = "Y")
    Q_learner <- make_learner(Lrnr_hal9001)
    g_learner <- make_learner(Lrnr_glm)

    #
    # qlib <- make_learner_stack(
    #   "Lrnr_mean",
    #   "Lrnr_glm_fast"
    # )
    #
    # glib <- make_learner_stack(
    #   "Lrnr_mean",
    #   "Lrnr_glm_fast"
    # )
    #
    # logit_metalearner <- make_learner(Lrnr_solnp, metalearner_logistic_binomial, loss_loglik_binomial)
    # Q_learner <- make_learner(Lrnr_sl, qlib, logit_metalearner)
    # g_learner <- make_learner(Lrnr_sl, glib, logit_metalearner)
    #
    learner_list <- list(Y = Q_learner, A = g_learner)

    suppressWarnings({
      fit <- tmle3(tmle_spec, data, node_list, learner_list)
    })
    results <- as.list(fit$summary[1])
    results$iteration <- iteration

    return(results)
  }

  single_time <- system.time({
    single_result <- one_sim()
  })

  # multiprocess test
  n_sim <- 20
  plan(multiprocess, workers = 4)
  multiprocess_time <- system.time({
    mpresults <- future.apply::future_lapply(X = seq_len(n_sim), one_sim)
  })

  # sequential comparison
  plan(sequential)
  sequential_time <- system.time({
    seqresults <- future.apply::future_lapply(X = seq_len(n_sim), one_sim)
  })
}
