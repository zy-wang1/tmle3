set.seed(1234)

data_sim <- generate_Zheng_data(B = 400, tau = timepoint, if_LY_misspec = if_misspec)
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

updater <- tmle3_Update_middle$new(maxit = 20, convergence_type = "scaled_var",
                                   fluctuation_type = "standard", submodel_type = "onestep", d_epsilon = 0.01, if_direction = T)
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)
tmle_params <- middle_spec$make_params(tmle_task, targeted_likelihood)
updater$tmle_params <- tmle_params

tmle_params[[1]]$estimates()$psi

record <- list()

updater$update_step(targeted_likelihood, tmle_task, update_fold)
record[[updater$step_number + 1]] <- ED_from_estimates(tmle_params %>% lapply(function(x) x$estimates()))

record



list_all_predicted_lkd <- tmle_params[[1]]$list_all_predicted_lkd
obs_data <- tmle_task$data %>% select(-c(id, t))
intervention_nodes <- union(names(tmle_params[[1]]$intervention_list_treatment), names(tmle_params$intervention_list_control))
intervention_variables <- map_chr(tmle_task$npsem[intervention_nodes], ~.x$variables)
intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
intervention_levels_treat <- map_dbl(tmle_params[[1]]$intervention_list_treatment, ~.x$value %>% as.character %>% as.numeric)
intervention_levels_control <- map_dbl(tmle_params[[1]]$intervention_list_control, ~.x$value %>% as.character %>% as.numeric)

loc_node <- 3

current_newH <- get_current_newH(loc_node,
                                 tmle_task, obs_data,
                                 intervention_variables, intervention_levels_treat, intervention_levels_control,
                                 list_all_predicted_lkd
)


temp <- list_all_predicted_lkd$R_1
temp[["newH"]] <- current_newH

left_join(obs_data[, 1:4], temp)$newH %>% head
tmle_params[[1]]$list_newH$fold_number$R_1 %>% head %>% as.vector()

tmle_params[[1]]$clever_covariates(update = T)
tmle_params[[1]]$estimates(update = T)
tmle_params[[1]]$list_newH$fold_number$R_1 %>% head %>% as.vector()
left_join(obs_data[, 1:4], temp)$newH %>% head


{
  ind_var <- loc_current_var <- which(obs_variable_names == tmle_task$npsem[[loc_node]]$variables)
  intervention_variables_loc_needed <- intervention_variables_loc[intervention_variables_loc < ind_var]

  # only update t!=0 and non-A nodes
  if (loc_node %in% c(loc_Z, loc_RLY)) {
    temp_current <- list_all_predicted_lkd[[loc_node]]
    Xt_input <- temp_current[ind_var]
    # loc_to_update <- Xt_input == 1
    # df_to_update <- temp_current[loc_to_update, ]  # the possible input where L_t = 1
    df_to_update <- df_to_update_0 <- temp_current
    df_to_update[[ind_var]] <- 1
    df_to_update_0[[ind_var]] <- 0

    # df_to_update_0 <- df_to_update
    # df_to_update_0[ind_var] <- 0  # replacing the above input with Lt = 0


    # get Q1 current
    {
      data_temp <- df_to_update[1:ind_var]
      # generate all possible 0, 1 valued rlz; set A to 0 first; set y_tau to always 1;
      all_possible_rlz_1 <- expand_values(obs_variable_names, to_drop = c(1:(ind_var) ),
                                          A = 1,
                                          rule_variables = c(last(obs_variable_names)), rule_values = c(1))
      all_possible_rlz_0 <- expand_values(obs_variable_names, to_drop = c(1:(ind_var) ),
                                          A = 0,
                                          rule_variables = c(last(obs_variable_names)), rule_values = c(1))
      # for each observed L_0 vector, generate all needed combinations, one version for A = 1, one version for A = 0
      unique_input <- data_temp[1:(ind_var)] %>% unique
      library_output <- data.frame(unique_input, output =
                                     map_dbl(1:nrow(unique_input), function(which_row) {
                                       # probs in the integrals, A=1 or A=0 is inserted
                                       temp_all_comb_0 <- temp_all_comb_1 <- unique_input[which_row, ]
                                       if (length(all_possible_rlz_0) > 0) {
                                         temp_all_comb_0 <- suppressWarnings(cbind(temp_all_comb_0, all_possible_rlz_0))
                                         temp_all_comb_1 <- suppressWarnings(cbind(temp_all_comb_1, all_possible_rlz_1))
                                       }
                                       if (length(intervention_variables_loc_needed) > 0) {
                                         for (i in 1:length(intervention_variables_loc_needed)) {
                                           temp_all_comb_0[, intervention_variables_loc[i]] <- intervention_levels_control[i]
                                           temp_all_comb_1[, intervention_variables_loc[i]] <- intervention_levels_treat[i]
                                         }
                                       }
                                       # for all non-A, non-0 variables, calculate the variable by rule
                                       # for Z's, use A = 0 values; outputs are predicted probs at each possible comb
                                       loc_Z_need <- loc_Z[loc_Z > loc_node]  # only integrate previous variables; debug: not including current
                                       temp_list_0 <- lapply(loc_Z_need,
                                                             function(each_t) {
                                                               left_join(temp_all_comb_0, list_all_predicted_lkd[[each_t]])$output
                                                             })
                                       loc_other_needed <- loc_RLY[loc_RLY > loc_node]
                                       temp_list_1 <- lapply(loc_other_needed,
                                                             function(each_t) {
                                                               left_join(temp_all_comb_1, list_all_predicted_lkd[[each_t]])$output
                                                             })
                                       temp_list <- c(temp_list_0, temp_list_1)
                                       pmap_dbl(temp_list, prod) %>% sum %>% return
                                     })
      )
      Q1_current <- left_join(data_temp, library_output)$output
      if (loc_node == length(temp_node_names)) Q1_current <- rep(1, length(Q1_current))
    }

    # get Q0 current
    {
      data_temp <- df_to_update_0[1:ind_var]
      # generate all possible 0, 1 valued rlz; set A to 0 first; set y_tau to always 1;
      all_possible_rlz_1 <- expand_values(obs_variable_names, to_drop = c(1:(ind_var) ),
                                          A = 1,
                                          rule_variables = c(last(obs_variable_names)), rule_values = c(1))
      all_possible_rlz_0 <- expand_values(obs_variable_names, to_drop = c(1:(ind_var) ),
                                          A = 0,
                                          rule_variables = c(last(obs_variable_names)), rule_values = c(1))
      # for each observed L_0 vector, generate all needed combinations, one version for A = 1, one version for A = 0
      unique_input <- data_temp[1:(ind_var)] %>% unique
      library_output <- data.frame(unique_input, output =
                                     map_dbl(1:nrow(unique_input), function(which_row) {
                                       # probs in the integrals, A=1 or A=0 is inserted
                                       temp_all_comb_0 <- temp_all_comb_1 <- unique_input[which_row, ]
                                       if (length(all_possible_rlz_0) > 0) {
                                         temp_all_comb_0 <- suppressWarnings(cbind(temp_all_comb_0, all_possible_rlz_0))
                                         temp_all_comb_1 <- suppressWarnings(cbind(temp_all_comb_1, all_possible_rlz_1))
                                       }
                                       if (length(intervention_variables_loc_needed) > 0) {
                                         for (i in 1:length(intervention_variables_loc_needed)) {
                                           temp_all_comb_0[, intervention_variables_loc[i]] <- intervention_levels_control[i]
                                           temp_all_comb_1[, intervention_variables_loc[i]] <- intervention_levels_treat[i]
                                         }
                                       }
                                       # for all non-A, non-0 variables, calculate the variable by rule
                                       # for Z's, use A = 0 values; outputs are predicted probs at each possible comb
                                       loc_Z_need <- loc_Z[loc_Z > loc_node]  # only integrate previous variables; debug: not including current
                                       temp_list_0 <- lapply(loc_Z_need,
                                                             function(each_t) {
                                                               left_join(temp_all_comb_0, list_all_predicted_lkd[[each_t]])$output
                                                             })
                                       loc_other_needed <- loc_RLY[loc_RLY > loc_node]
                                       temp_list_1 <- lapply(loc_other_needed,
                                                             function(each_t) {
                                                               left_join(temp_all_comb_1, list_all_predicted_lkd[[each_t]])$output
                                                             })
                                       temp_list <- c(temp_list_0, temp_list_1)
                                       pmap_dbl(temp_list, prod) %>% sum %>% return
                                     })
      )
      Q0_current <- left_join(data_temp, library_output)$output
      if (loc_node == length(temp_node_names)) Q0_current <- rep(0, length(Q0_current))
    }

    # get H current
    {
      data_temp <- temp_current[1:ind_var]  # take =1 probs
      temp_ind <- ind_var

      all_observed_1 <- all_observed_0 <- data_temp
      if (length(intervention_variables_loc_needed) > 0) {
        for (i in 1:length(intervention_variables_loc_needed)) {
          all_observed_0[, intervention_variables_loc[i]] <- intervention_levels_control[i]
          all_observed_1[, intervention_variables_loc[i]] <- intervention_levels_treat[i]
        }
      }

      # R and L (including Y before t = tau) (or not Z not A not t=0, not Y_tau) can be calculated in the same way
      if (loc_node %in% loc_RLY) {
        loc_A_needed <- loc_A[loc_A < loc_node]  # all needed A nodes
        loc_Z_needed <- loc_Z[loc_Z < loc_node]  # all needed Z nodes

        # all predicted probs at these locs are needed
        # needed inputs are in either all_observed_1 or all_observed_0

        # these A probs will be taken as product
        part_A <- lapply(loc_A_needed, function(k) {
          left_join(all_observed_1, list_all_predicted_lkd[[k]])$output
        }) %>% pmap_dbl(prod)

        # these ratios Z probs will be taken as product
        part_Z <- ifelse_vec(length(loc_Z_needed) == 0, rep(1, length(part_A)),
                             lapply(loc_Z_needed, function(k) {
                               left_join(all_observed_0, list_all_predicted_lkd[[k]])$output /
                                 left_join(all_observed_1, list_all_predicted_lkd[[k]])$output
                             }) %>% pmap_dbl(prod))

        H_current <- ifelse(
          # data_temp[last(intervention_variables_loc_needed)] == 1
          apply(data_temp[intervention_variables_loc_needed] == 1, 1, prod) == 1
          , 1/part_A*part_Z, 0) %>% as.vector
      }
      # Z nodes
      if (loc_node %in% loc_Z) {
        loc_A_needed <- loc_A[loc_A < loc_node]  # all needed A nodes
        loc_RLY_needed <- loc_RLY[loc_RLY < loc_node]

        # these A probs will be taken as product
        part_A <- lapply(loc_A_needed, function(k) {
          left_join(all_observed_0, list_all_predicted_lkd[[k]])$output
        }) %>% pmap_dbl(prod)

        # these ratios Z probs will be taken as product
        part_LR <- lapply(loc_RLY_needed, function(k) {
          left_join(all_observed_1, list_all_predicted_lkd[[k]])$output /
            left_join(all_observed_0, list_all_predicted_lkd[[k]])$output
        }) %>% pmap_dbl(prod)

        H_current <- ifelse(
          # ZW todo: stochastic intervention
          # data_temp[last(intervention_variables_loc_needed)] == 0
          apply(data_temp[intervention_variables_loc_needed] == 0, 1, prod) == 1
          , 1/part_A*part_LR, 0)
      }
    }

    # get newH current
    newH_current <- H_current * (Q1_current - Q0_current)
    # return(newH_current)
  } else {
    # return(NULL)
  }
}


temp[["Q1"]] <- Q1_current
temp[["Q0"]] <- Q0_current
temp[["H"]] <- H_current

cf_task_treatment <- tmle_params[[1]]$cf_likelihood_treatment$enumerate_cf_tasks(tmle_task)[[1]]
cf_task_control <- tmle_params[[1]]$cf_likelihood_control$enumerate_cf_tasks(tmle_task)[[1]]

list_H <- get_obs_H(tmle_task, obs_data, current_likelihood = targeted_likelihood,
                    cf_task_treatment, cf_task_control,
                    intervention_variables, intervention_levels_treat, intervention_levels_control)
list_Q_1 <- get_obs_Q(tmle_task, obs_data, list_H,
                      intervention_variables, intervention_levels_treat, intervention_levels_control,
                      list_all_predicted_lkd,
                      lt = 1)
list_Q_0 <- get_obs_Q(tmle_task, obs_data, list_H,
                      intervention_variables, intervention_levels_treat, intervention_levels_control,
                      list_all_predicted_lkd,
                      lt = 0)


identical(left_join(obs_data[, 1:4], temp)$Q1, list_Q_1[[3]])
identical(left_join(obs_data[, 1:4], temp)$Q0, list_Q_0[[3]])
identical(left_join(obs_data[, 1:4], temp)$H, list_H[[3]])

identical(list_H[[3]]*(list_Q_1[[3]] - list_Q_0[[3]]),
          left_join(obs_data[, 1:4], temp)$newH
)

identical(left_join(obs_data[, 1:4], temp)$Q1, list_Q_1[[3]])
identical(left_join(obs_data[, 1:4], temp)$Q0, list_Q_0[[3]])
identical(left_join(obs_data[, 1:4], temp)$H, list_H[[3]])



left_join(obs_data[, 1:4], temp)$Q0 %>% tail
list_Q_0[[3]] %>% tail
list_Q[[3]] %>% tail
lt <- 0
loc_node <- 3
left_join(obs_data[, 1:(loc_current_var-1)], library_output)$output[1:10]

left_join(obs_data[, 1:4], temp)$H %>% head
list_H[[3]] %>% head




