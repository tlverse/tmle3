context("Direct ATE for single node interventions")

library(sl3)
library(tmle3)
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
  W = c("sexn"),
  A = "parity01",
  Y = "haz01"
)

# generate diffeernt data
data2 <- copy(data)
data2$haz01 <- sample(data2$haz01)
tmle_spec <- tmle_ATE(1, 0)

# generate another copy of the same data
data3 <- copy(data)

# define data
tmle_task <- tmle_spec$make_tmle_task(data, node_list)
tmle_task2 <- tmle_spec$make_tmle_task(data2, node_list)
tmle_task3 <- tmle_spec$make_tmle_task(data3, node_list)
cf_task <- tmle_task$generate_counterfactual_task(uuid = "cf", new_data = data.table(A = 0))

test_that("Distinct data yields distinct uuids", expect_true(tmle_task$uuid != tmle_task2$uuid))
test_that("Identical data yields idetical uuids", expect_true(tmle_task$uuid == tmle_task3$uuid))
test_that("CF data yields distinct uuids", expect_true(tmle_task$uuid != cf_task$uuid))
