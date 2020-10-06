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
# first-step logistic
###

updater <- tmle3_Update_middle$new(maxit = 1, convergence_type = "scaled_var",
                                   fluctuation_type = "standard", submodel_type = "logistic"
                                   # , cvtmle = T
)
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)
tmle_params <- middle_spec$make_params(tmle_task, targeted_likelihood)
updater$tmle_params <- tmle_params
suppressMessages(
  test <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)
)

firsttargeting <- test$estimates[[1]]

temp_lmed3_est <- firsttargeting$psi
temp_IC <- firsttargeting$IC
var_D <- var(temp_IC)
n <- length(temp_IC)
se <- sqrt(var_D / n)
# temp_var <- var(temp_IC)/length(temp_IC)
CI2_first <- firsttargeting$psi + 1.96 * se
CI1_first <- firsttargeting$psi - 1.96 * se




###
# iterative logistic, maxit = 10
###

updater <- tmle3_Update_middle$new(maxit = 10, convergence_type = "scaled_var",
                                   fluctuation_type = "standard", submodel_type = "logistic"
                                   # , cvtmle = T
)
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)
tmle_params <- middle_spec$make_params(tmle_task, targeted_likelihood)
updater$tmle_params <- tmle_params

suppressMessages(
  test <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)
)
tentargeting <- test$estimates[[1]]

temp_lmed3_est <- tentargeting$psi
temp_IC <- tentargeting$IC
var_D <- var(temp_IC)
n <- length(temp_IC)
se <- sqrt(var_D / n)
# temp_var <- var(temp_IC)/length(temp_IC)
CI2_tentargeting <- tentargeting$psi + 1.96 * se
CI1_tentargeting <- tentargeting$psi - 1.96 * se


###
# one-step grid search, fixed direction
###

updater <- tmle3_Update_middle$new(maxit = 100, convergence_type = "scaled_var",
                                   fluctuation_type = "standard", submodel_type = "onestep", d_epsilon = 0.01)
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)
tmle_params <- middle_spec$make_params(tmle_task, targeted_likelihood)
updater$tmle_params <- tmle_params

suppressMessages(
  test <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)
)
onestep_targeting <- test$estimates[[1]]
temp_lmed3_est <- onestep_targeting$psi
temp_IC <- onestep_targeting$IC

var_D <- var(temp_IC)
n <- length(temp_IC)
se <- sqrt(var_D / n)


# temp_var <- var(temp_IC)/length(temp_IC)
CI2_onestep <- onestep_targeting$psi + 1.96 * se
CI1_onestep <- onestep_targeting$psi - 1.96 * se





###
# one-step grid search, each part of EIC & each step can have different directions
###

updater <- tmle3_Update_middle$new(maxit = 100, convergence_type = "scaled_var",
                                   fluctuation_type = "standard", submodel_type = "onestep", d_epsilon = 0.01, if_direction = T)
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




###
# one-step grid search, projection-based; direction should be flexible
###

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
node_list$L_1 <- "L_1"
names(data_sim[[2]])[grep("L1_1", names(data_sim[[2]]))] <- "L_1"
names(data_wide)[grep("L1_1", names(data_wide))] <- "L_1"
node_list$L_2 <- "L_2"
names(data_sim[[3]])[grep("L1_2", names(data_sim[[3]]))] <- "L_2"
names(data_wide)[grep("L1_2", names(data_wide))] <- "L_2"



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

tlik <- Targeted_Likelihood$new(initial_likelihood,
                                submodel_type_by_node = "EIC" ,
                                updater = list(convergence_type = "scaled_var",
                                               constrain_step = T,
                                               optim_delta_epsilon = F,
                                               one_dimensional=TRUE,
                                               delta_epsilon = 0.01,
                                               maxit = 100))
tmle_params <- middle_spec$make_params(tmle_task, tlik, if_projection = T, initial_likelihood)
tmle_params[[1]]$estimates()$psi

# aaa <- tlik$cache$tasks %>% names
# tlik$cache$tasks %>% length
# tlik$cache$tasks %>% lapply(function(x) x$data %>% dim)
# tlik$updater$epsilons

suppressMessages(
  tlik$updater$update(tlik, tmle_task)
)

# tlik$updater$update_step(tlik, tmle_task)
#
# tlik$cache$tasks %>% length
#
# bbb <- tlik$cache$tasks %>% names
# setdiff(bbb, aaa)
# tlik$cache$tasks$`5bcaf7aa013cc544cdfe854fc5be0b83`$data$A_2 %>% table

onestep_projected <- tmle_params[[1]]$estimates()

onestep_projected_est <- onestep_projected$psi

temp_IC <- onestep_projected$IC
var_D <- var(temp_IC)
n <- length(temp_IC)
se <- sqrt(var_D / n)

# temp_var <- var(temp_IC)/length(temp_IC)
CI2_onestep_projected <- onestep_projected_est + 1.96 * se
CI1_onestep_projected <- onestep_projected_est - 1.96 * se



temp_result <- list(nontargeting_est = nontargeting$psi,
                    nontargeting_ci1 = CI1,
                    nontargeting_ci2 = CI2,
                    first_est = firsttargeting$psi,
                    first_ci1 = CI1_first,
                    first_ci2 = CI2_first,
                    ten_est = tentargeting$psi,
                    ten_ci1 = CI1_tentargeting,
                    ten_ci2 = CI2_tentargeting,
                    onestep_fixed_est = onestep_targeting$psi,
                    onestep_fixed_ci1 = CI1_onestep,
                    onestep_fixed_ci2 = CI2_onestep,
                    onestep_dir_est = onestep_targeting_dir$psi,
                    onestep_dir_ci1 = CI1_onestep_dir,
                    onestep_dir_ci2 = CI2_onestep_dir,
                    onestep_pro_est = onestep_projected_est,
                    onestep_pro_ci1 = CI1_onestep_projected,
                    onestep_pro_ci2 = CI2_onestep_projected
)

if (!dir.exists("./temp_output")) dir.create("./temp_output")

saveRDS(temp_result, "./temp_output/1.RDS")

