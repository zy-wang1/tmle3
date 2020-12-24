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

if_misspec <- T

data_sim <- generate_Zheng_data(B = 1000, tau = 2, if_LY_misspec = if_misspec)
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

# record_time <- c()
# time_count <- 1


###
# nontargeting
###
# start.time <- Sys.time()
tmle_params <- med_spec$make_params(tmle_task, initial_likelihood, static_likelihood = NULL)
suppressMessages(
  nontargeting1 <- tmle_params[[1]]$estimates(tmle_task)
)
temp_lmed3_nontargeting <- nontargeting1$psi
temp_IC <- nontargeting1$IC
var_D <- var(temp_IC)
n <- length(temp_IC)
se <- sqrt(var_D / n)
CI2_1 <- temp_lmed3_nontargeting + 1.96 * se
CI1_1 <- temp_lmed3_nontargeting - 1.96 * se
# end.time <- Sys.time()
# record_time[time_count] <- end.time - start.time
# time_count <- time_count + 1

###
# projected
###

# projected-EIC, resampling
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
tlik <- Targeted_Likelihood$new(initial_likelihood,
                                submodel_type_by_node = "EIC" ,
                                updater = updater)

# projection param
tmle_params <- med_spec$make_params(tmle_task, tlik, options = list("tc"), static_likelihood = initial_likelihood,
                                    if_projection = T
                                    , n_resampling = 50000
)
updater$tmle_params <- tmle_params

suppressMessages(
  nontargeting_proj_1 <- tmle_params[[1]]$estimates(tmle_task)
)
temp_IC <- nontargeting_proj_1$IC
var_D <- var(temp_IC)
n <- length(temp_IC)
se <- sqrt(var_D / n)
CI2_proj_1 <- nontargeting_proj_1$psi + 1.96 * se
CI1_proj_1 <- nontargeting_proj_1$psi - 1.96 * se

capture.output(
  suppressMessages(
    new_fit <- fit_tmle3(tmle_task, tlik, tmle_params, updater)
  )
)
new_fit
# tlik$cache$tasks %>% lapply(function(x) x$data %>% dim)
new_psi <- new_fit$summary$tmle_est
CI1_new <- new_fit$summary$lower
CI2_new <- new_fit$summary$upper
updater$step_number


###
# one-step, maxit=100, de=0.01
###
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






temp_result <- list(nontargeting1$psi,
                    CI1_1,
                    CI2_1,
                    nontargeting_proj_1$psi,
                    CI1_proj_1,
                    CI2_proj_1,
                    new_psi,
                    CI1_new,
                    CI2_new,
                    est_onestep = fit_onestep_100$summary$tmle_est,
                    ci1_onestep = fit_onestep_100$summary$lower,
                    ci2_onestep = fit_onestep_100$summary$upper
)


# temp_result[[length(temp_result) + 1]] <- record_time
# names(temp_result)[length(temp_result)] <- "time"

if (!dir.exists("./temp_output")) dir.create("./temp_output")

saveRDS(temp_result, "./temp_output/1.RDS")

