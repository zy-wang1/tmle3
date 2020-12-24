timepoint <- 2
if_misspec <- T
# data_truth <- generate_Zheng_data(B = 100000, tau = timepoint, seed = 202008, setAM = c(1, 0), if_LY_misspec = if_misspec)
# truth <- data_truth[[timepoint + 1]]$Y %>% mean
# truth

# set.seed(1234)

data_sim <- generate_Zheng_data(B = 1000, tau = timepoint, if_LY_misspec = if_misspec)
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




# updater <- tmle3_Update_middle$new(maxit = 1, convergence_type = "sample_size",
#                                    fluctuation_type = "standard", submodel_type = "logistic"
#                                    # , cvtmle = T
# )
# targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)
# tmle_params <- middle_spec$make_params(tmle_task, targeted_likelihood)
# updater$tmle_params <- tmle_params
#
# tmle_params[[1]]$estimates()
# tmle_params[[1]]$list_D %>% lapply(mean)
# suppressMessages(
#   test <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)
# )
# test
# tmle_params[[1]]$list_D %>% lapply(mean)




updater <- tmle3_Update_middle$new(maxit = 20, convergence_type = "scaled_var",
                                   fluctuation_type = "standard", submodel_type = "onestep", d_epsilon = 0.001
                                   , if_direction = T
)
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)
tmle_params <- middle_spec$make_params(tmle_task, targeted_likelihood)
updater$tmle_params <- tmle_params




intervention_nodes <- union(names(tmle_params[[1]]$intervention_list_treatment), names(tmle_params[[1]]$intervention_list_control))


temp_node_names <- names(tmle_task$npsem)
loc_A <- grep("A", temp_node_names)
loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)

# get list of all possible predicted lkds
obs_data <- tmle_task$data %>% as.data.frame %>% dplyr::select(-c(id, t))
obs_variable_names <- colnames(obs_data)

cf_task_control <- tmle_params[[1]]$cf_likelihood_control$enumerate_cf_tasks(tmle_task)[[1]]
cf_task_treatment <- tmle_params[[1]]$cf_likelihood_treatment$enumerate_cf_tasks(tmle_task)[[1]]

intervention_variables <- map_chr(tmle_task$npsem[intervention_nodes], ~.x$variables)
intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
intervention_levels_treat <- map_dbl(tmle_params[[1]]$intervention_list_treatment, ~.x$value %>% as.character %>% as.numeric)
intervention_levels_control <- map_dbl(tmle_params[[1]]$intervention_list_control, ~.x$value %>% as.character %>% as.numeric)
names(intervention_levels_treat) <- names(tmle_params[[1]]$intervention_list_treatment)
names(intervention_levels_control) <- names(tmle_params[[1]]$intervention_list_control)

list_H <- get_obs_H(tmle_task, obs_data, current_likelihood = targeted_likelihood,
                    cf_task_treatment, cf_task_control,
                    intervention_variables, intervention_levels_treat, intervention_levels_control)

