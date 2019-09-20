context("Missingness processing")
library(sl3)
data(cpp)
data <- as.data.table(cpp)
data$parity01 <- as.numeric(data$parity > 0)
data$parity01_fac <- factor(data$parity01)
data$haz01 <- as.numeric(data$haz > 0)


node_list <- list(
  W = c(
    "apgar1", "apgar5", "gagebrth", "mage",
    "meducyrs", "sexn"
  ),
  A = "parity01",
  Y = "haz01"
)


processed <- process_missing(data, node_list)
# debugonce(process_missing)
processed2 <- process_missing(processed$data, processed$node_list)
test_that("process_missing is idempotent (ignoring column ordering)", {
  expect_equal(processed$data, processed2$data[, names(processed$data), with = FALSE])
  expect_equal(processed$node_list, processed2$node_list)
})
