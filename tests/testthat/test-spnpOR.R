context("spnpOR test")

passes <- c()
passes1 <- c()
for (i in 1:1) {
  print(i)
  library(sl3)
  n <- 500
  W <- runif(n, -1, 1)
  A <- rbinom(n, size = 1, prob = plogis(0))
  Y <- rbinom(n, size = 1, prob = plogis(A + W + A * W))
  quantile(plogis(1 + W) * (1 - plogis(1 + W)) / (plogis(W) * (1 - plogis(W))))
  data <- data.table(W, A, Y)
  lrnr_Y0W <- Lrnr_glm$new()
  lrnr_A <- Lrnr_glm$new()
  node_list <- list(W = "W", A = "A", Y = "Y")
  learner_list <- list(A = lrnr_A, Y = lrnr_Y0W)
  spec_spCATE <- tmle3_Spec_spCausalGLM$new(~ 1 + W, "OR")
  suppressWarnings(out <- tmle3(spec_spCATE, data, node_list, learner_list = learner_list))
  out <- out$summary
  passes <- cbind(passes, out$lower <= 1 & out$upper >= 1)


  spec_spCATE <- tmle3_Spec_npCausalGLM$new(~ 1 + W, "OR")
  suppressWarnings(out <- tmle3(spec_spCATE, data, node_list, learner_list = learner_list))
  out <- out$summary
  passes1 <- cbind(passes1, out$lower <= 1 & out$upper >= 1)

  print(rowMeans(passes))
  print(rowMeans(passes1))
}
