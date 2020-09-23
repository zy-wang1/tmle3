library(R6)  # R6class
library(data.table)  # setDT
library(sl3)
library(digest)
library(uuid)  # UUIDgenerate
library(delayed)  # bundle_delayed
library(assertthat)  # assert_that
# library(methods)  # is


code_list <- list.files("./R", full.names = T)
for (code in code_list) source(code)
source("./temp_code/generate_data.R")

data_sim <- generate_Zheng_data(B = 400, tau = 2, if_LY_misspec = F)
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

tmle_params <- middle_spec$make_params(tmle_task, initial_likelihood)

nontargeting <- tmle_params[[1]]$estimates(tmle_task)
temp_lmed3_nontargeting <- nontargeting$psi
temp_IC <- nontargeting$IC

var_D <- var(temp_IC)
n <- length(temp_IC)
se <- sqrt(var_D / n)
# temp_var <- var(temp_IC)/length(temp_IC)
CI2 <- temp_lmed3_nontargeting + 1.96 * se
CI1 <- temp_lmed3_nontargeting - 1.96 * se

CI1
CI2
temp_lmed3_nontargeting






updater <- tmle3_Update_middle$new(maxit = 1, convergence_type = "scaled_var",
                            fluctuation_type = "standard", submodel_type = "logistic"
                            # , cvtmle = T
)
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)
tmle_params <- middle_spec$make_params(tmle_task, targeted_likelihood)
updater$tmle_params <- tmle_params
test <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)
firsttargeting <- test$estimates[[1]]
temp_lmed3_est <- firsttargeting$psi
temp_IC <- firsttargeting$IC

var_D <- var(temp_IC)
n <- length(temp_IC)
se <- sqrt(var_D / n)


# temp_var <- var(temp_IC)/length(temp_IC)
CI2_first <- firsttargeting$psi + 1.96 * se
CI1_first <- firsttargeting$psi - 1.96 * se

temp_lmed3_est
CI1_first
CI2_first
