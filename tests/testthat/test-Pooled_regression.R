context("Perform pooled regression of longitudinal data structure")

library(simcausal)
# Generate long format data
D <- DAG.empty()
D <- D +
  node("L0", distr = "rbinom", size =1, prob = 0.5) +
  node("A0", distr = "rbinom", size = 1, prob = 0.5 + L0/4 ) +
  node("L1",  distr =  "rbinom", size = 1, prob = 0.5 + A0/4 ) +
  node("A1", distr = "rbinom", size = 1, prob = 0.5 + L1/4 ) +
  node("L2",  distr =  "rbinom", size = 1, prob = 0.5 + A1/4 ) +
  node("A2", distr = "rbinom", size = 1, prob = 0.5 + L2/4 ) +
  node("Y", distr = "rbinom", size = 1, prob = 0.5 + (A2 + L2)/4 )

setD <- set.DAG(D)
dat <- sim(setD, n = 1e3)
dat$id <- dat$ID
dat$ID <- NULL



long_data <- rbindlist(lapply(0:2, function(i){
  dat <- copy(dat[, c(grep(paste(i), colnames(dat), value = T), "id", "Y")])
  setnames(dat, c("L", "A", "id", "Y"))
  dat$t <- i
  return(dat)
}))
setkey(long_data, id , t)


# Longitudinal npsem with growing past
npsem <- list(define_node("L0", c("L"),c(), time = 0), define_node("A0", c("A"), c("L0"), time = 0) ,
              define_node("L1", c("L"),c("A0", "L1"), time = 1), define_node("A1", c("A"), c("L0", "L1", "A0"), time = 1) ,
              define_node("L2", c("L"),c("L1", "L0", "A1", "A0"), time = 2), define_node("A2", c("A"), c("L0", "L1", "L2", "A0", "A1"), time = 2), define_node("Y", c("Y"), c("L0", "L1", "L2", "A0", "A1", "A2"), time = 2))

task <- tmle3_Task$new(long_data, npsem, time = "t", id = "id")


# Longitudinal npsem with fixed past to allow pooling
npsem <- list(define_node("L0", c("L"),c(), time = 0), define_node("A0", c("A"), c("L0"), time = 0) ,
              define_node("L1", c("L"),c("A0", "L0"), time = 1), define_node("A1", c("A"), c("L1", "A0"), time = 1) ,
              define_node("L2", c("L"),c("L1", "A1"), time = 2), define_node("A2", c("A"), c( "L2","A1"), time = 2), define_node("Y", c("Y"), c("L0", "L1", "L2", "A0", "A1", "A2"), time = 2))

task <- tmle3_Task$new(long_data, npsem, time = "t", id = "id")


#Check names are equal up to order for pooling (though they should be the same including order)
assertthat::assert_that(identical(sort(names(task$get_regression_task("L2")$data)), sort(names(task$get_regression_task("L1")$data))))
assertthat::assert_that(identical(sort(names(task$get_regression_task("A2")$data)), sort(names(task$get_regression_task("A1")$data))))

factor_list <- list(LF_emp$new("L0"),
                    LF_fit$new("A0", make_learner(Lrnr_glm_fast)),
                    LF_fit$new(c("L1", "L2"), make_learner(Lrnr_glm_fast)),
                    LF_fit$new(c("A1", "A2"), make_learner(Lrnr_glm_fast)),
                    LF_fit$new("Y", make_learner(Lrnr_glm_fast))
)

# Do pooled training
lik <- Likelihood$new(factor_list)
assertthat::assert_that(length(lik$factor_list_pooled) == length(factor_list))
assertthat::assert_that(length(lik$factor_list) == length(task$npsem))
lik_trained <- lik$train(task)

task <- tmle3_Task$new(long_data, npsem, time = "t", id = "id")

# Check pooled regression tasks coincide (up to order of rows due to id and time)
pooled = task$get_regression_task(c("L1", "L2"), expand = T)$data
stacked = rbindlist(list(task$get_regression_task(c("L1"), expand = T, is_time_variant =  T)$data,
                         task$get_regression_task(c("L2"), expand = T, is_time_variant = T)$data), use.names = F)
setkey(stacked, id , t)
setkey(pooled, id , t)
assertthat::assert_that(all(stacked == pooled, use.names = T))


# Check likelihoods
l1 <- lik$get_likelihood(task, "L1", drop_id = F, drop_time = T, to_wide = F)
l2 <- lik$get_likelihood(task, "L2", drop_id = F, drop_time = T, to_wide = F)
l12 <- lik$get_likelihood(task, c("L1", "L2"), drop_id = F, drop_time = F, to_wide = T)
merged_l12 = merge(l1, l2, by = c("id"))
# Check that pooled and unpooled predictions coincide
assertthat::assert_that(all(l12[, c("L1", "L2")] ==  merged_l12[, c("L1", "L2")]))
