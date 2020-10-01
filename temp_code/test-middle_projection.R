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

set.seed(123)

data_sim <- generate_Zheng_data(B = 500, tau = 1, if_LY_misspec = F)
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
)

node_list$L_1 <- "L_1"
names(data_sim[[2]])[grep("L1_1", names(data_sim[[2]]))] <- "L_1"
names(data_wide)[grep("L1_1", names(data_wide))] <- "L_1"

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

# tmle_params <- middle_spec$make_params(tmle_task, initial_likelihood, if_projection = T)
# tmle_params[[1]]$clever_covariates()


# tlik <- Targeted_Likelihood$new(initial_likelihood,
#                                 submodel_type_by_node = "EIC" ,
#                                 updater = list(one_dimensional = T,
#                                                constrain_step = T,
#                                                convergence_type = "scaled_var",
#                                                delta_epsilon = list("Y" = 1e-3, "A" = function(x) {
#                                                  res = 0.75/max(abs(x))
#                                                  res <- min(res, 0.1)
#                                                  return(res)
#                                                },
#                                                "Z" = function(x) {
#                                                  res = 0.75/max(abs(x))
#                                                  res <- min(res, 0.1)
#                                                  return(res)
#                                                })))
tlik <- Targeted_Likelihood$new(initial_likelihood,
                                submodel_type_by_node = "EIC" ,
                                updater = list(convergence_type = "sample_size",
                                               constrain_step = T,
                                               optim_delta_epsilon = F,
                                               one_dimensional=TRUE,
                                               delta_epsilon = 0.01))
tlik$cache$tasks %>% length

tmle_params <- middle_spec$make_params(tmle_task, tlik, if_projection = T, initial_likelihood)

# tlik$cache$tasks %>% length
#
#
# tmle_params[[1]]$clever_covariates()
# tmle_params[[1]]$estimates()
# tmle_params[[1]]$gradient$compute_component(tmle_task, "L_1", fold_number = "validation")
#
# tlik$updater$update_step(tlik, tmle_task)
# tlik$cache$tasks %>% length
#
#
# tmle_params[[1]]$gradient$compute_component(tmle_task, "L_1", fold_number = "full")


tlik$updater$update(tlik, tmle_task)


