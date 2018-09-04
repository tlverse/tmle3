context("CV-Likelihood: Fold-specific Likelihood estimates")

set.seed(1234)
library(sl3)
# library(tmle3)
library(uuid)
library(assertthat)
library(data.table)
library(future)
# setup data for test
data(cpp)
data <- as.data.table(cpp)
data$parity01 <- as.numeric(data$parity > 0)
data$parity01_fac <- factor(data$parity01)
data$haz01 <- as.numeric(data$haz > 0)
data[is.na(data)] <- 0
node_list <- list(
  W = c(
    "apgar1", "apgar5", "gagebrth", "mage",
    "meducyrs", "sexn"
  ),
  A = "parity01",
  Y = "haz01"
)

qlib <- make_learner_stack(
  "Lrnr_mean",
  "Lrnr_glm_fast"
)

glib <- make_learner_stack(
  "Lrnr_mean",
  "Lrnr_glm_fast"
)

metalearner <- make_learner(Lrnr_nnls)
Q_learner <- make_learner(Lrnr_sl, qlib, metalearner)
g_learner <- make_learner(Lrnr_sl, glib, metalearner)
learner_list <- list(Y = Q_learner, A = g_learner)
tmle_spec <- tmle_TSM_all()

# define data
tmle_task <- tmle_spec$make_tmle_task(data, node_list)

# define likelihood
initial_likelihood <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)

# debugonce(initial_likelihood$get_likelihood)
EY_val <- initial_likelihood$get_likelihoods(tmle_task, "Y", 0)
EY_direct_val <- initial_likelihood$factor_list[["Y"]]$get_likelihood(tmle_task, 0)
EY_direct_val2 <- initial_likelihood$factor_list[["Y"]]$get_likelihood(tmle_task, 0)
# debugonce(initial_likelihood$factor_list[["Y"]]$get_mean)
EY_direct_val2 <- initial_likelihood$factor_list[["Y"]]$get_likelihood(tmle_task, 0)
EY_resub <- initial_likelihood$get_likelihoods(tmle_task, "Y", -1)
EY_direct_resub <- initial_likelihood$factor_list[["Y"]]$get_likelihood(tmle_task, -1)

Q_task <- tmle_task$get_regression_task("Y")
Q_task2 <- tmle_task$get_regression_task("Y")
Q_fit <- initial_likelihood$factor_list[["Y"]]$learner

preds1 <- Q_fit$predict_fold(Q_task, 0 )
preds12 <- Q_fit$predict_fold(Q_task, 0 )

# debugonce(Q_fit$fit_object$cv_meta_fit$base_predict)
preds2 <- Q_fit$predict_fold(Q_task2, 0 )

meta_task1 <- Q_fit$fit_object$cv_fit$chain_fold(Q_task, 0)
meta_task12 <- Q_fit$fit_object$cv_fit$chain_fold(Q_task, 0)
meta_task2 <- Q_fit$fit_object$cv_fit$chain_fold(Q_task2, 0)

preds1_direct <- Q_fit$fit_object$cv_meta_fit$base_predict(meta_task1)
all.equal(preds1, preds1_direct)

preds12_direct <- Q_fit$fit_object$cv_meta_fit$predict(meta_task12)
all.equal(preds12, preds12_direct)

preds2_direct <- Q_fit$fit_object$cv_meta_fit$predict(meta_task2)
all.equal(preds2, preds2_direct)

# debugonce(Q_fit$predict_fold)
preds13 <- Q_fit$predict_fold(Q_task, 0 )
preds_mat <- as.vector(as.matrix(meta_task2$X)%*%Q_fit$fit_object$cv_meta_fit$coefficients)
head(cbind(preds1, preds12, preds2, preds1_direct, preds12_direct, preds2_direct, preds13, preds_mat))


all.equal(preds1, preds2)
all.equal(preds1, preds12)

all.equal(preds1_direct, preds2_direct)
all.equal(preds1_direct, preds12_direct)

all.equal(meta_task1$data, meta_task12$data)
all.equal(meta_task1$data, meta_task2$data)

head(cbind(EY_val, EY_direct_val, EY_direct_val2, EY_resub, EY_direct_resub))

all.equal(EY,EY_direct)
