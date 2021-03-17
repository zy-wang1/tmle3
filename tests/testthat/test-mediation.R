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

data_sim <- generate_Zheng_data(B = 1000, tau = 1, if_LY_misspec = F)
data_wide <- data.frame(data_sim)

node_list <- list(L_0 = c("L1_0", "L2_0"),
                  A_1 = "A_1",
                  R_1 = "R_1",
                  Z_1 = "Z_1",
                  L_1 = "L1_1",
                  Y_1 = "Y_1"
                  # ,
                  # A_2 = "A_2",
                  # R_2 = "R_2",
                  # Z_2 = "Z_2",
                  # L_2 = "L1_2",
                  # Y_2 = "Y_2"
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

# define RI-based mediation parameter
tmle_params <- mediation_spec$make_params(tmle_task, initial_likelihood, options = "tc")

# test estimate
nontargeting <- tmle_params[[1]]$estimates(tmle_task)
temp_lmed3_nontargeting <- nontargeting$psi
temp_IC <- nontargeting$IC
var_D <- var(temp_IC)
n <- length(temp_IC)
se <- sqrt(var_D / n)
CI2 <- temp_lmed3_nontargeting + 1.96 * se
CI1 <- temp_lmed3_nontargeting - 1.96 * se


# test update
n_subject <- nrow(tmle_task$data)
tlik <- Targeted_Likelihood$new(initial_likelihood,
                                submodel_type_by_node = "EIC" ,
                                updater = list(convergence_type = "scaled_var",
                                               constrain_step = T,
                                               optim_delta_epsilon = F,  # fixed small step_size
                                               one_dimensional=T,
                                               delta_epsilon=function(x) {
                                                 ifelse(abs(mean(x %>% as.vector)) < sqrt(var(x %>% as.vector)/n_subject)/log(n_subject),
                                                        0,
                                                        ifelse(mean(x %>% as.vector) > 0, 0.01, -0.01)
                                                 )
                                               },
                                               maxit=100
                                               ,
                                               cvtmle=F
                                ))
tmle_params <- mediation_spec$make_params(tmle_task, tlik)
# tmle_params[[1]]$estimates()$psi
tmle_params[[1]]$clever_covariates()$IC %>% colMeans()
suppressMessages(
  tlik$updater$update(tlik, tmle_task)
)
suppressWarnings(suppressMessages(
  new_est <- tmle_params[[1]]$estimates()$psi
))
tmle_params[[1]]$clever_covariates()$IC %>% colMeans()
new_est

onestep_test <- tmle_params[[1]]$estimates()
onestep_test_est <- onestep_test$psi
temp_IC <- onestep_test$IC
var_D <- var(temp_IC)
n <- length(temp_IC)
se <- sqrt(var_D / n)
CI2_onestep_test <- onestep_test_est + 1.96 * se
CI1_onestep_test <- onestep_test_est - 1.96 * se

rbind(c(temp_lmed3_nontargeting, CI1, CI2),
      c(onestep_test_est, CI1_onestep_test, CI2_onestep_test)
)
