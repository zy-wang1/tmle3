library(R6)  # R6class
library(data.table)  # setDT
library(sl3)
library(digest)
library(uuid)  # UUIDgenerate
library(delayed)  # bundle_delayed
library(assertthat)  # assert_that
library(speedglm)  # speedglm
# library(methods)  # is

library(dplyr)  # dplyr::select, to mask other packages

code_list <- list.files("./R", full.names = T)
for (code in code_list) source(code)
source("./temp_code/generate_data.R")

# record <- record_cor <- list()

set.seed(1234)

data_sim <- generate_Zheng_data(B = 1000, tau = 2, if_LY_misspec = T)
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
)

middle_spec <- tmle_middle(
  treatment_level = 1,
  control_level = 0
)

tmle_task <- middle_spec$make_tmle_task(data_wide, node_list)

# choose base learners
lrnr_glm_fast <- Lrnr_glm_fast$new(outcome_type = "binomial")
learner_list <- lapply(1:length(tmle_task$npsem), function(s) Lrnr_sl$new(
  learners = list(
    lrnr_glm_fast
  )
))
names(learner_list) <- names(tmle_task$npsem)  # the first will be ignored; empirical dist. will be used for covariates

initial_likelihood <- middle_spec$make_initial_likelihood(
  tmle_task,
  learner_list
)



###
# nontargeting
###

tmle_params <- middle_spec$make_params(tmle_task, initial_likelihood, static_likelihood = NULL)
suppressMessages(
  nontargeting <- tmle_params[[1]]$estimates(tmle_task)
)

temp_lmed3_nontargeting <- nontargeting$psi
temp_IC <- nontargeting$IC
var_D <- var(temp_IC)
n <- length(temp_IC)
se <- sqrt(var_D / n)
# temp_var <- var(temp_IC)/length(temp_IC)
CI2 <- temp_lmed3_nontargeting + 1.96 * se
CI1 <- temp_lmed3_nontargeting - 1.96 * se




###
# one-step grid search, each part of EIC & each step can have different directions
###

updater <- tmle3_Update_middle$new(maxit = 100, convergence_type = "scaled_var",
                                   fluctuation_type = "standard", submodel_type = "onestep", d_epsilon = 0.0001, if_direction = T,
                                   cvtmle = F)
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)
tmle_params <- middle_spec$make_params(tmle_task, targeted_likelihood)
updater$tmle_params <- tmle_params

suppressMessages(
  test <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)
)
onestep_targeting_dir <- test$estimates[[1]]
temp_lmed3_est <- onestep_targeting_dir$psi
temp_IC <- onestep_targeting_dir$IC

var_D <- var(temp_IC)
n <- length(temp_IC)
se <- sqrt(var_D / n)

# temp_var <- var(temp_IC)/length(temp_IC)
CI2_onestep_dir <- onestep_targeting_dir$psi + 1.96 * se
CI1_onestep_dir <- onestep_targeting_dir$psi - 1.96 * se

onestep_steps <- test$steps
onestep_ED <- ED_from_estimates(tmle_params %>% lapply(function(x) x$estimates()))
onestep_if_local_max <- updater$if_local_max

se_Dstar <- sqrt(var(temp_IC)
                 / tmle_task$nrow
)
ED_threshold <- se_Dstar / min(log(tmle_task$nrow), 10)
onestep_threshold <- ED_threshold

temp_result <- list(nontargeting_est = nontargeting$psi,
                    nontargeting_ci1 = CI1,
                    nontargeting_ci2 = CI2,
                    onestep_dir_est = onestep_targeting_dir$psi,
                    onestep_dir_ci1 = CI1_onestep_dir,
                    onestep_dir_ci2 = CI2_onestep_dir,
                    onestep_dir_steps = onestep_steps,
                    onestep_dir_ED = onestep_ED,
                    onestep_dir_threshold = onestep_threshold,
                    onestep_dir_if_local = onestep_if_local_max
)

if (!dir.exists("./temp_output")) dir.create("./temp_output")

saveRDS(temp_result, "./temp_output/1.RDS")

