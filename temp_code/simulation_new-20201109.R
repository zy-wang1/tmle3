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
                  # ,
                  # A_3 = "A_3",
                  # R_3 = "R_3",
                  # Z_3 = "Z_3",
                  # L_3 = "L1_3",
                  # Y_3 = "Y_3"
)
node_list$L_1 <- "L_1"
names(data_sim[[2]])[grep("L1_1", names(data_sim[[2]]))] <- "L_1"
names(data_wide)[grep("L1_1", names(data_wide))] <- "L_1"
node_list$L_2 <- "L_2"
names(data_sim[[3]])[grep("L1_2", names(data_sim[[3]]))] <- "L_2"
names(data_wide)[grep("L1_2", names(data_wide))] <- "L_2"
# node_list$L_3 <- "L_3"
# names(data_sim[[4]])[grep("L1_3", names(data_sim[[4]]))] <- "L_3"
# names(data_wide)[grep("L1_3", names(data_wide))] <- "L_3"

med_spec <- tmle_med(
  treatment_level = 1,
  control_level = 0
)

tmle_task <- med_spec$make_tmle_task(data_wide, node_list)

# choose base learners
lrnr_mean <- make_learner(Lrnr_mean)
# lrnr_glm <- make_learner(Lrnr_glm)
lrnr_glm_fast <- Lrnr_glm_fast$new(outcome_type = "binomial")
# lrnr_ranger50 <- make_learner(Lrnr_ranger, num.trees = 50)
lrnr_hal_simple <- make_learner(Lrnr_hal9001, max_degree = 3, n_folds = 10)
lrnr_lasso <- make_learner(Lrnr_glmnet) # al  pha default is 1
lrnr_ridge <- make_learner(Lrnr_glmnet, alpha = 0)
lrnr_elasticnet <- make_learner(Lrnr_glmnet, alpha = .5)
learner_list <- lapply(1:length(tmle_task$npsem), function(s) Lrnr_sl$new(
  learners = list(
    # lrnr_mean,
    # lrnr_glm,
    lrnr_glm_fast
    # ,
    # lrnr_ranger50,
    # lrnr_hal_simple
    # ,
    # lrnr_lasso, lrnr_ridge, lrnr_elasticnet
  )
))
# learner_list <- lapply(1:length(tmle_task$npsem), function(s) lrnr_glm_fast)  # simplest learner
names(learner_list) <- names(tmle_task$npsem)  # the first will be ignored; empirical dist. will be used for covariates

initial_likelihood <- med_spec$make_initial_likelihood(
  tmle_task,
  learner_list
)

record_time <- c()
time_count <- 1


###
# nontargeting
###
start.time <- Sys.time()
tmle_params <- med_spec$make_params(tmle_task, initial_likelihood, static_likelihood = NULL)
suppressMessages(
  nontargeting <- tmle_params[[1]]$estimates(tmle_task)
)
temp_lmed3_nontargeting <- nontargeting$psi
temp_IC <- nontargeting$IC
var_D <- var(temp_IC)
n <- length(temp_IC)
se <- sqrt(var_D / n)
CI2 <- temp_lmed3_nontargeting + 1.96 * se
CI1 <- temp_lmed3_nontargeting - 1.96 * se
end.time <- Sys.time()
record_time[time_count] <- end.time - start.time
time_count <- time_count + 1

###
# logistic, maxit=1
###
start.time <- Sys.time()
updater <- tmle3_Update$new(maxit = 1,
                            convergence_type = "scaled_var",
                            fluctuation_type = "standard",
                            cvtmle = F
)
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater,
                                               submodel_type_by_node = "logistic")
tmle_params <- med_spec$make_params(tmle_task, targeted_likelihood)
updater$tmle_params <- tmle_params
suppressMessages(
  fit_logistic_1 <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)
)
end.time <- Sys.time()
record_time[time_count] <- end.time - start.time
time_count <- time_count + 1

###
# logistic, maxit=1, cvtmle
###
start.time <- Sys.time()
updater <- tmle3_Update$new(maxit = 1,
                            convergence_type = "scaled_var",
                            fluctuation_type = "standard",
                            cvtmle = T
)
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater,
                                               submodel_type_by_node = "logistic")
tmle_params <- med_spec$make_params(tmle_task, targeted_likelihood)
updater$tmle_params <- tmle_params
suppressMessages(
  fit_logistic_1_cvtmle <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)
)
end.time <- Sys.time()
record_time[time_count] <- end.time - start.time
time_count <- time_count + 1

###
# logistic, maxit=10
###
start.time <- Sys.time()
updater <- tmle3_Update$new(maxit = 10,
                            convergence_type = "scaled_var",
                            fluctuation_type = "standard",
                            cvtmle = F
)
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater,
                                               submodel_type_by_node = "logistic")
tmle_params <- med_spec$make_params(tmle_task, targeted_likelihood)
updater$tmle_params <- tmle_params
suppressMessages(
  fit_logistic_10 <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)
)
end.time <- Sys.time()
record_time[time_count] <- end.time - start.time
time_count <- time_count + 1

###
# logistic, maxit=10, cvtmle
###
start.time <- Sys.time()
updater <- tmle3_Update$new(maxit = 10,
                            convergence_type = "scaled_var",
                            fluctuation_type = "standard",
                            cvtmle = T
)
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater,
                                               submodel_type_by_node = "logistic")
tmle_params <- med_spec$make_params(tmle_task, targeted_likelihood)
updater$tmle_params <- tmle_params
suppressMessages(
  fit_logistic_10_cvtmle <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)
)
end.time <- Sys.time()
record_time[time_count] <- end.time - start.time
time_count <- time_count + 1

###
# one-step, maxit=100, de=0.01
###
start.time <- Sys.time()
updater <- tmle3_Update$new(convergence_type = "scaled_var",
                            constrain_step = T,
                            optim_delta_epsilon = T,
                            one_dimensional=F,
                            delta_epsilon=function(x) {
                              ifelse(abs(mean(x %>% as.vector)) < sqrt(var(x %>% as.vector)/tmle_task$nrow)/log(tmle_task$nrow),
                                     0.00000001,
                                     ifelse(mean(x %>% as.vector) > 0, 0.01, -0.01)
                              )
                            },
                            maxit=100
                            ,
                            cvtmle=F)
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood,
                                submodel_type_by_node = "EIC" ,
                                updater = updater)
tmle_params <- med_spec$make_params(tmle_task, targeted_likelihood)
updater$tmle_params <- tmle_params
suppressMessages(
  fit_onestep_100 <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)
)
end.time <- Sys.time()
record_time[time_count] <- end.time - start.time
time_count <- time_count + 1

###
# one-step, maxit=100, de=0.01, cvtmle
###
start.time <- Sys.time()
updater <- tmle3_Update$new(convergence_type = "scaled_var",
                            constrain_step = T,
                            optim_delta_epsilon = T,
                            one_dimensional=F,
                            delta_epsilon=function(x) {
                              ifelse(abs(mean(x %>% as.vector)) < sqrt(var(x %>% as.vector)/tmle_task$nrow)/log(tmle_task$nrow),
                                     0.00000001,
                                     ifelse(mean(x %>% as.vector) > 0, 0.01, -0.01)
                              )
                            },
                            maxit=100
                            ,
                            cvtmle=T)
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood,
                                               submodel_type_by_node = "EIC" ,
                                               updater = updater)
tmle_params <- med_spec$make_params(tmle_task, targeted_likelihood)
updater$tmle_params <- tmle_params
suppressMessages(
  fit_onestep_100_cvtmle <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)
)
end.time <- Sys.time()
record_time[time_count] <- end.time - start.time
time_count <- time_count + 1

# ###
# # one-step, maxit=250, de=0.004
# ###
# start.time <- Sys.time()
# updater <- tmle3_Update_middle$new(maxit = 250, convergence_type = "scaled_var",
#                                    fluctuation_type = "standard", submodel_type = "onestep", d_epsilon = 0.004, if_direction = T,
#                                    cvtmle = F)
# targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)
# tmle_params <- middle_spec$make_params(tmle_task, targeted_likelihood)
# updater$tmle_params <- tmle_params
# suppressMessages(
#   fit_onestep_1000 <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)
# )
# end.time <- Sys.time()
# record_time[time_count] <- end.time - start.time
# time_count <- time_count + 1

# ###
# # one-step, maxit=250, de=0.004, cvtmle
# ###
# start.time <- Sys.time()
# updater <- tmle3_Update_middle$new(maxit = 250, convergence_type = "scaled_var",
#                                    fluctuation_type = "standard", submodel_type = "onestep", d_epsilon = 0.004, if_direction = T,
#                                    cvtmle = T)
# targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)
# tmle_params <- middle_spec$make_params(tmle_task, targeted_likelihood)
# updater$tmle_params <- tmle_params
# suppressMessages(
#   fit_onestep_1000_cvtmle <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)
# )
# end.time <- Sys.time()
# record_time[time_count] <- end.time - start.time
# time_count <- time_count + 1

###
# projection one-step, maxit=100, de=0.01
###
start.time <- Sys.time()
n_subject <- nrow(tmle_task$data)
tlik <- Targeted_Likelihood$new(initial_likelihood,
                                submodel_type_by_node = "EIC" ,
                                updater = list(convergence_type = "scaled_var",
                                               constrain_step = T,
                                               optim_delta_epsilon = T,
                                               one_dimensional=T,
                                               delta_epsilon=function(x) {
                                                 ifelse(abs(mean(x %>% as.vector)) < sqrt(var(x %>% as.vector)/n_subject)/log(n_subject),
                                                        0.00000001,
                                                        ifelse(mean(x %>% as.vector) > 0, 0.01, -0.01)
                                                 )
                                               },
                                               maxit=100
                                               ,
                                               cvtmle=F
                                ))
tmle_params <- med_spec$make_params(tmle_task, tlik, if_projection = T, static_likelihood = initial_likelihood)
tmle_params[[1]]$estimates()$psi
capture.output(
  tlik$updater$update(tlik, tmle_task)
)
onestep_projected <- tmle_params[[1]]$estimates()
onestep_projected_est <- onestep_projected$psi
temp_IC <- onestep_projected$IC
var_D <- var(temp_IC)
n <- length(temp_IC)
se <- sqrt(var_D / n)
CI2_onestep_projected <- onestep_projected_est + 1.96 * se
CI1_onestep_projected <- onestep_projected_est - 1.96 * se
end.time <- Sys.time()
record_time[time_count] <- end.time - start.time
time_count <- time_count + 1


###
# projection one-step, maxit=100, de=0.01, cvtmle
###
start.time <- Sys.time()
tlik <- Targeted_Likelihood$new(initial_likelihood,
                                submodel_type_by_node = "EIC" ,
                                updater = list(convergence_type = "scaled_var",
                                               constrain_step = T,
                                               optim_delta_epsilon = T,
                                               one_dimensional=T,
                                               delta_epsilon=function(x) {
                                                 ifelse(abs(mean(x %>% as.vector)) < sqrt(var(x %>% as.vector)/n_subject)/log(n_subject),
                                                        0.00000001,
                                                        ifelse(mean(x %>% as.vector) > 0, 0.01, -0.01)
                                                 )
                                               },
                                               maxit=100
                                               ,
                                               cvtmle=T
                                ))
tmle_params <- med_spec$make_params(tmle_task, tlik, if_projection = T, static_likelihood = initial_likelihood)
tmle_params[[1]]$estimates()$psi
capture.output(
  tlik$updater$update(tlik, tmle_task)
)
onestep_projected <- tmle_params[[1]]$estimates()
onestep_projected_est_cvtmle <- onestep_projected$psi
temp_IC <- onestep_projected$IC
var_D <- var(temp_IC)
n <- length(temp_IC)
se <- sqrt(var_D / n)
CI2_onestep_projected_cvtmle <- onestep_projected$psi + 1.96 * se
CI1_onestep_projected_cvtmle <- onestep_projected$psi - 1.96 * se
end.time <- Sys.time()
record_time[time_count] <- end.time - start.time
time_count <- time_count + 1


###
# projection one-step, maxit=100, de=0.01
###
start.time <- Sys.time()
n_subject <- nrow(tmle_task$data)
tlik <- Targeted_Likelihood$new(initial_likelihood,
                                submodel_type_by_node = "EIC" ,
                                updater = list(convergence_type = "scaled_var",
                                               constrain_step = T,
                                               optim_delta_epsilon = T,
                                               one_dimensional=T,
                                               delta_epsilon=function(x) {
                                                 ifelse(abs(mean(x %>% as.vector)) < sqrt(var(x %>% as.vector)/n_subject)/log(n_subject),
                                                        0.00000001,
                                                        ifelse(mean(x %>% as.vector) > 0, 0.01, -0.01)
                                                 )
                                               },
                                               maxit=100
                                               ,
                                               cvtmle=F
                                ))
tmle_params <- med_spec$make_params(tmle_task, tlik, if_projection = T, static_likelihood = initial_likelihood, n_resampling = 10000)
tmle_params[[1]]$estimates()$psi
capture.output(
  tlik$updater$update(tlik, tmle_task)
)
onestep_re <- tmle_params[[1]]$estimates()
onestep_re_est <- onestep_projected$psi
temp_IC <- onestep_re$IC
var_D <- var(temp_IC)
n <- length(temp_IC)
se <- sqrt(var_D / n)
CI2_onestep_re <- onestep_re_est + 1.96 * se
CI1_onestep_re <- onestep_re_est - 1.96 * se
end.time <- Sys.time()
record_time[time_count] <- end.time - start.time
time_count <- time_count + 1




temp_result <- list(nontargeting_est = nontargeting$psi,
                    nontargeting_ci1 = CI1,
                    nontargeting_ci2 = CI2,
                    fit_logistic_1_est = fit_logistic_1$summary$tmle_est,
                    fit_logistic_1_ci1 = fit_logistic_1$summary$lower,
                    fit_logistic_1_ci2 = fit_logistic_1$summary$upper,
                    fit_logistic_1_cvtmle_est = fit_logistic_1_cvtmle$summary$tmle_est,
                    fit_logistic_1_cvtmle_ci1 = fit_logistic_1_cvtmle$summary$lower,
                    fit_logistic_1_cvtmle_ci2 = fit_logistic_1_cvtmle$summary$upper,
                    fit_logistic_10_est = fit_logistic_10$summary$tmle_est,
                    fit_logistic_10_ci1 = fit_logistic_10$summary$lower,
                    fit_logistic_10_ci2 = fit_logistic_10$summary$upper,
                    fit_logistic_10_cvtmle_est = fit_logistic_10_cvtmle$summary$tmle_est,
                    fit_logistic_10_cvtmle_ci1 = fit_logistic_10_cvtmle$summary$lower,
                    fit_logistic_10_cvtmle_ci2 = fit_logistic_10_cvtmle$summary$upper,
                    fit_onestep_100_est = fit_onestep_100$summary$tmle_est,
                    fit_onestep_100_ci1 = fit_onestep_100$summary$lower,
                    fit_onestep_100_ci2 = fit_onestep_100$summary$upper,
                    fit_onestep_100_cvtmle_est = fit_onestep_100_cvtmle$summary$tmle_est,
                    fit_onestep_100_cvtmle_ci1 = fit_onestep_100_cvtmle$summary$lower,
                    fit_onestep_100_cvtmle_ci2 = fit_onestep_100_cvtmle$summary$upper,
                    # fit_onestep_1000_est = fit_onestep_1000$summary$tmle_est,
                    # fit_onestep_1000_ci1 = fit_onestep_1000$summary$lower,
                    # fit_onestep_1000_ci2 = fit_onestep_1000$summary$upper,
                    # fit_onestep_1000_cvtmle_est = fit_onestep_1000_cvtmle$summary$tmle_est,
                    # fit_onestep_1000_cvtmle_ci1 = fit_onestep_1000_cvtmle$summary$lower,
                    # fit_onestep_1000_cvtmle_ci2 = fit_onestep_1000_cvtmle$summary$upper,
                    onestep_projected_est = onestep_projected_est,
                    CI1_onestep_projected = CI1_onestep_projected,
                    CI2_onestep_projected = CI2_onestep_projected,
                    onestep_projected_est_cvtmle = onestep_projected_est_cvtmle,
                    CI1_onestep_projected_cvtmle = CI1_onestep_projected_cvtmle,
                    CI2_onestep_projected_cvtmle = CI2_onestep_projected_cvtmle
                    ,
                    onestep_re_est = onestep_re_est,
                    CI1_onestep_re = CI1_onestep_re,
                    CI2_onestep_re = CI2_onestep_re
)
temp_result[[length(temp_result) + 1]] <- record_time
names(temp_result)[length(temp_result)] <- "time"

if (!dir.exists("./temp_output")) dir.create("./temp_output")

saveRDS(temp_result, "./temp_output/1.RDS")

