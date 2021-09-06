context("spCATE, npCATE, npCATT test")


passes <- c()
passes1 <- c()
passes2 <- c()

for (i in 1:1) {
  print(i)

  n <- 500
  W <- runif(n, -1, 1)
  A <- rbinom(n, size = 1, prob = plogis(W))
  Y <- rnorm(n, mean = A + W, sd = 0.3)
  data <- data.table(W, A, Y)
  lrnr_Y0W <- Lrnr_gam$new()
  lrnr_A <- Lrnr_gam$new()

  node_list <- list(W = "W", A = "A", Y = "Y")
  learner_list <- list(A = lrnr_A, Y = lrnr_Y0W, var_Y = Lrnr_mean$new())
  # spec_spCATE <- tmle3_Spec_spCausalGLM$new(~1, "CATE")
  # out <- tmle3(spec_spCATE, data, node_list, learner_list = learner_list)
  spec_spCATE <- tmle3_Spec_npCausalGLM$new(~1, "CATE")
  suppressWarnings(out <- tmle3(spec_spCATE, data, node_list, learner_list = learner_list))
  out <- out$summary
  passes <- c(passes, out$lower <= 1 & out$upper >= 1)


  spec_spCATE <- tmle3_Spec_npCausalGLM$new(~1, "CATT")
  suppressWarnings(out <- tmle3(spec_spCATE, data, node_list, learner_list = learner_list))
  out <- out$summary
  passes1 <- c(passes1, out$lower <= 1 & out$upper >= 1)


  spec_spCATE <- tmle3_Spec_spCausalGLM$new(~1, "CATE")
  suppressWarnings(out <- tmle3(spec_spCATE, data, node_list, learner_list = learner_list))
  out <- out$summary
  passes2 <- c(passes2, out$lower <= 1 & out$upper >= 1)

  print(mean(passes))
  print(mean(passes1))
  print(mean(passes2))
}
