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

set.seed(1234)

data_sim <- generate_Zheng_data(B = 2000, tau = 2, if_LY_misspec = F)
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

initial_likelihood <- middle_spec$make_initial_likelihood(
  tmle_task,
  learner_list
)

tlik <- Targeted_Likelihood$new(initial_likelihood,
                                submodel_type_by_node = "EIC" ,
                                updater = list(convergence_type = "sample_size", constrain_step = T, delta_epsilon = 0.01))


tmle_params <- middle_spec$make_params(tmle_task, tlik, if_projection = T)


tlik$updater$update(tlik, tmle_task)
