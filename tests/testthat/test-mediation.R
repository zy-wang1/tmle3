context("Mediation target parameters, v2")

library(R6)  # R6class
library(data.table)  # setDT
library(sl3)
library(digest)
library(uuid)  # UUIDgenerate
library(delayed)  # bundle_delayed
library(assertthat)  # assert_that
# library(speedglm)  # avoid speedglm; sometimes mess up column ordering
# library(methods)  # is
library(dplyr)  # dplyr::select, to mask other packages
library(purrr)  # map_ functions

code_list <- list.files("./R", full.names = T)
for (code in code_list) source(code)
# source("./temp_code/generate_data.R")  # included in mediationp_utils.R

set.seed(1234)

data_sim <- generate_Zheng_data(B = 1000, tau = 2, if_LY_misspec = F)
data_wide <- data.frame(data_sim)

node_list <- list(L_0 = c("L1_0", "L2_0"),
                  A_1 = "A_1",
                  R_1 = "R_1",
                  Z_1 = "Z_1",
                  L_1 = "L1_1",
                  Y_1 = "Y_1"
                  ,
                  A_2 = "A_2",
                  R_2 = "R_2",
                  Z_2 = "Z_2",
                  L_2 = "L1_2",
                  Y_2 = "Y_2"
                  # ,
                  # A_3 = "A_3",
                  # R_3 = "R_3",
                  # Z_3 = "Z_3",
                  # L_3 = "L1_3",
                  # Y_3 = "Y_3"
)

mediation_spec <- tmle_mediation(
  treatment_level = 1,
  control_level = 0
)
tmle_task <- mediation_spec$make_tmle_task(data_wide, node_list)

# choose base learners
lrnr_mean <- make_learner(Lrnr_mean)
lrnr_glm <- make_learner(Lrnr_glm, outcome_type = "binomial")
# lrnr_glm_fast <- Lrnr_glm_fast$new(outcome_type = "binomial")
# lrnr_ranger50 <- make_learner(Lrnr_ranger, num.trees = 50)
lrnr_hal_simple <- make_learner(Lrnr_hal9001, max_degree = 3, n_folds = 10)
lrnr_lasso <- make_learner(Lrnr_glmnet) # al  pha default is 1
lrnr_ridge <- make_learner(Lrnr_glmnet, alpha = 0)
lrnr_elasticnet <- make_learner(Lrnr_glmnet, alpha = .5)
learner_list <- lapply(1:length(tmle_task$npsem), function(s) Lrnr_sl$new(
  learners = list(
    # lrnr_mean,
    lrnr_glm
    # ,
    # lrnr_ranger50,
    # lrnr_hal_simple,
    # lrnr_lasso, lrnr_ridge, lrnr_elasticnet
  )
))  # simplest learner
names(learner_list) <- names(tmle_task$npsem)  # the first will be ignored; empirical dist. will be used for covariates

# define obs data initial likelihood
initial_likelihood <- mediation_spec$make_initial_likelihood(
  tmle_task,
  learner_list
)

