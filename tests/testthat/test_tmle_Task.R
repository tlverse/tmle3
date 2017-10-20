context("test_tmle_Task.R")

data(cpp)
cpp <- cpp[!is.na(cpp[, "haz"]), ]
cpp$parity3 <- cpp$parity
cpp$parity3[cpp$parity3>3] <- 3
cpp$parity01 <- as.numeric(cpp$parity>0)
cpp[is.na(cpp)] <- 0

W_nodes <- c("apgar1", "apgar5", "gagebrth", "mage", "meducyrs", "sexn")
A_node <- "parity01"
Y_node <- "haz"
tmle_nodes=list(W=W_nodes,
                A=A_node,
                Y=Y_node)

task = tmle_Task$new(cpp, tmle_nodes=tmle_nodes)
npsem <- task$npsem
glm_fast <- make_learner(Lrnr_glm_fast)
learner_list <- list(Y=glm_fast, A=glm_fast)
likelihood <- make_learner(Likelihood, learner_list)
likelihood_fit <- likelihood$train(task)
likelihood_preds <- likelihood_fit$base_predict()
likelihood_preds
