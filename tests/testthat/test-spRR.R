
context("spRR test")




passes <- c()
for(i in 1:1){
  print(i)
  n <- 500
  W <- runif(n, -1, 1)
  A <- rbinom(n, size = 1, prob = plogis(W))
  Y <-  rpois(n, exp(A  + A*W + W))
  data <- data.table(W,A,Y)
  data
  lrnr_Y0W <- Lrnr_glmnet$new(family = "poisson")
  lrnr_A <- Lrnr_glm$new()

  node_list <- list (W = "W", A = "A", Y= "Y")
  learner_list <- list(A  = lrnr_A, Y = lrnr_Y0W)
  spec_spCATE <- tmle3_Spec_spCausalGLM$new(~1 + W, "RR")
  out <- suppressWarnings(tmle3(spec_spCATE, data, node_list, learner_list = learner_list) )
  out <- out$summary
  passes <- cbind(passes , out$lower <= 1 & out$upper >= 1)
  print(rowMeans(passes))

}
