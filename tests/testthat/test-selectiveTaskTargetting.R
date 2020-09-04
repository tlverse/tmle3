context("Target specific nodes only for tasks")

library(simcausal)
library(sl3)
library(tmle3)

# Generate long format data
D <- DAG.empty()
D <- D +
  node("W1", distr = "runif", min = 0, max = 100) +
  node("W", distr = "rconst", const = round(W1)) +
  node("A", distr = "rbinom", size = 1, prob = ( 0.4 + W/300) ) +
  node("Y",  distr =  "rbinom", size = 1, prob = (0.4 + A/4 - W/300 ) )


setD <- set.DAG(D)
dat <- sim(setD, n = 3000)
dat$id <- dat$ID
dat$ID <- NULL
dat <- as.data.table(dat)

npsem <- list(define_node("W", "W", variable_type = variable_type("continuous")),
              define_node("A", "A", "W"),
              define_node("Y", "Y", c("A", "W")))

suppressWarnings(task <- tmle3_Task$new(dat, npsem, id = "id"))
setattr(task, "target_nodes", c("Y"))
attr(task, "target_nodes")



factor_list <- list(LF_emp$new("W"),
                    LF_fit$new("A", make_learner(Lrnr_glm)),
                    LF_fit$new("Y", make_learner(Lrnr_glm), type = "mean"))

lik <- Likelihood$new(factor_list)
lik<-lik$train(task)
tlik <- Targeted_Likelihood$new(lik)
#cache tasks and get initial likelihoods
v1 = tlik$get_likelihood(task, "Y")
param <- Param_ATE$new(tlik, intervention_list_control = list(LF_static$new("A", value = 0)), intervention_list_treatment= list(LF_static$new("A", value = 1)))

tlik$updater$update_step(tlik, task)
tlik$updater$update_step(tlik, task)
tlik$updater$update_step(tlik, task)


v2 = tlik$get_likelihood(task, "Y")

if(all(v1==v2)){
  stop("Nope 1")
}



suppressWarnings(task <- tmle3_Task$new(dat, npsem, id = "id"))
setattr(task, "target_nodes", c("None"))
attr(task, "target_nodes")

lik <- Likelihood$new(factor_list)
lik<-lik$train(task)
tlik <- Targeted_Likelihood$new(lik)
#cache tasks and get initial likelihoods
v1 = tlik$get_likelihood(task, "Y")
param <- Param_ATE$new(tlik, intervention_list_control = list(LF_static$new("A", value = 0)), intervention_list_treatment= list(LF_static$new("A", value = 1)))

tlik$updater$update_step(tlik, task)



suppressWarnings(v2 <- tlik$get_likelihood(task, "Y", check_sync = F))

if(any(v1!=v2)){
  stop("Nope 2")
}












