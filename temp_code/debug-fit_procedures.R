# timepoint <- 1
# if_misspec <- T
# data_truth <- generate_Zheng_data(B = 100000, tau = timepoint, seed = 202008, setAM = c(1, 0), if_LY_misspec = if_misspec)
# truth <- data_truth[[timepoint + 1]]$Y %>% mean
# truth
#
# set.seed(1234)

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

{
  if (is.null(tmle_task)) {
    tmle_task <- tmle_params[[1]]$observed_likelihood$training_task
  }

  intervention_nodes <- union(names(tmle_params[[1]]$intervention_list_treatment), names(tmle_params[[1]]$intervention_list_control))

  # todo: extend for stochastic
  cf_task_treatment <- tmle_params[[1]]$cf_likelihood_treatment$enumerate_cf_tasks(tmle_task)[[1]]
  cf_task_control <- tmle_params[[1]]$cf_likelihood_control$enumerate_cf_tasks(tmle_task)[[1]]

  Y <- tmle_task$get_tmle_node(tmle_params[[1]]$outcome_node)

  # all not A, not t=0 nodes
  temp_node_names <- names(tmle_task$npsem)
  loc_A <- grep("A", temp_node_names)
  loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
  loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
  if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)

  # get list of all possible predicted lkds
  obs_data <- tmle_task$data %>% as.data.frame %>% dplyr::select(-c(id, t))
  obs_variable_names <- colnames(obs_data)
  # ZW todo: to handle long format and wide format

  param.list_all_predicted_lkd <- tmle_params[[1]]$observed_likelihood$list_all_predicted_lkd
  # only calculate list of lkd here when it is null; otherwise only update it in updater$apply_update_all
  if (!is.null(param.list_all_predicted_lkd)) {
    list_all_predicted_lkd <- param.list_all_predicted_lkd
  } else {
    list_all_predicted_lkd <- lapply(1:length(temp_node_names), function(loc_node) {
      if (loc_node > 1) {
        # currently only support univariate node for t>0
        current_variable <- tmle_task$npsem[[loc_node]]$variables
        temp_input <- expand_values(variables = obs_variable_names[1:which(obs_variable_names == current_variable)])  # all possible inputs
        temp_task <- tmle3_Task$new(temp_input, tmle_task$npsem[1:loc_node])
        temp_output <- self$get_likelihood(temp_task, node = temp_node_names[loc_node], fold_number = fold_number)  # corresponding outputs
        data.frame(temp_input, output = temp_output) %>% return
      }
    })
    private$.list_all_predicted_lkd <- list_all_predicted_lkd
  }

  # recalculate if list_D or result is null, or if we force it to update
  # make sure we force update this after each updating step
  # this helps speed up updater$check_convergence
  if (!is.null(param.list_D) & !is.null(param.result) & update == F) {
    # return(param.result)
    param.result
  } else {
    intervention_variables <- map_chr(tmle_task$npsem[intervention_nodes], ~.x$variables)
    intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
    intervention_levels_treat <- map_dbl(tmle_params[[1]]$intervention_list_treatment, ~.x$value %>% as.character %>% as.numeric)
    intervention_levels_control <- map_dbl(tmle_params[[1]]$intervention_list_control, ~.x$value %>% as.character %>% as.numeric)
    # nodes to integrate out in the target identification
    # only support univaraite node for now; assume treatment level is one
    all_possible_RZLY_1 <- expand_values(obs_variable_names, to_drop = c(1:length(tmle_task$npsem[[1]]$variables) ),
                                         rule_variables = c(intervention_variables,
                                                            last(obs_variable_names)),
                                         rule_values = c(intervention_levels_treat, 1))
    all_possible_RZLY_0 <- expand_values(obs_variable_names, to_drop = c(1:length(tmle_task$npsem[[1]]$variables) ),
                                         rule_variables = c(intervention_variables,
                                                            last(obs_variable_names)),
                                         rule_values = c(intervention_levels_control, 1))

    # for each observed L_0 vector, generate all needed combinations, one version for A = 1, one version for A = 0
    unique_L0 <- obs_data[, tmle_task$npsem[[1]]$variables] %>% unique
    library_L0 <- data.frame(unique_L0, output =
                               map_dbl(1:nrow(unique_L0), function(which_row) {
                                 temp_all_comb_0 <- cbind(unique_L0[which_row, ], all_possible_RZLY_0)
                                 temp_all_comb_1 <- cbind(unique_L0[which_row, ], all_possible_RZLY_1)
                                 # for all non-A, non-0 variables, calculate the variable by rule
                                 # for Z's, use A = 0 values; outputs are predicted probs at each possible comb
                                 # note that list_all_predicted_lkd is ordered by node
                                 temp_list_0 <- lapply(loc_Z,
                                                       function(each_t) {
                                                         left_join(temp_all_comb_0, list_all_predicted_lkd[[each_t]])$output
                                                       })
                                 temp_list_1 <- lapply(loc_RLY,
                                                       function(each_t) {
                                                         left_join(temp_all_comb_1, list_all_predicted_lkd[[each_t]])$output
                                                       })
                                 temp_list <- c(temp_list_0, temp_list_1)
                                 pmap_dbl(temp_list, prod) %>% sum %>% return
                               })
    )
    # substitution estimator
    vec_est <- left_join(obs_data[, tmle_task$npsem[[1]]$variables], library_L0)$output
    psi <- mean(vec_est)

    # get true IC

    # # getting H's
    # all_observed_1 <- all_observed_0 <- obs_data
    # for (temp_A in intervention_variables) {
    #   all_observed_1 <- all_observed_1 %>% mutate(!!temp_A := 1)
    #   all_observed_0 <- all_observed_0 %>% mutate(!!temp_A := 0)
    # }
    #
    # list_H <- get_obs_H(tmle_task, obs_data, current_likelihood = self$observed_likelihood,
    #                     cf_task_treatment, cf_task_control,
    #                     intervention_variables, intervention_levels_treat, intervention_levels_control)
    # # get a list of needed deltaQ
    # list_Q_1 <- get_obs_Q(tmle_task, obs_data, list_H,
    #                       intervention_variables, intervention_levels_treat, intervention_levels_control,
    #                       list_all_predicted_lkd,
    #                       lt = 1)
    # list_Q_0 <- get_obs_Q(tmle_task, obs_data, list_H,
    #                       intervention_variables, intervention_levels_treat, intervention_levels_control,
    #                       list_all_predicted_lkd,
    #                       lt = 0)

    list_newH <- tmle_params[[1]]$clever_covariates(tmle_task, fold_number)
    list_newH_raw <- tmle_params[[1]]$list_newH_raw$fold_number

    list_D <- list_D_trt <- list_D_ctrl <- list()
    for (ind_var in 1:length(list_newH)) {
      if(!is.null(list_newH[[ind_var]])) {
        # ZW todo: for discretized variables
        current_ind <- (obs_data[[tmle_task$npsem[[ind_var]]$variables]] == 1)*1
        if (ind_var %in% loc_Z) temp_p <- tmle_params[[1]]$observed_likelihood$get_likelihoods(cf_task_control, temp_node_names[ind_var]) else
          temp_p <- tmle_params[[1]]$observed_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[ind_var])
        temp_p <- ifelse(current_ind == 1, temp_p, 1 - temp_p)
        list_D[[ind_var]] <- (current_ind - temp_p) *list_newH[[ind_var]]

        # e.g. L nodes, ind_A all 1's for treat cf task
        if (ind_var %in% loc_Z) {
          list_D_trt[[ind_var]] <- rep(0, tmle_task$nrow)
          list_D_ctrl[[ind_var]] <- (current_ind - temp_p) *list_newH_raw[[ind_var]]
        } else {
          list_D_trt[[ind_var]] <- (current_ind - temp_p) *list_newH_raw[[ind_var]]
          list_D_ctrl[[ind_var]] <- rep(0, tmle_task$nrow)
        }
      }
    }
    list_D[[1]] <- list_D_trt[[1]] <- list_D_ctrl[[1]] <- vec_est - psi
    # # for debug the last node
    # list_D[[length(list_D)]] <- length(list_D) <- length(list_D) <- rep(0, vec_est)
    names(list_D) <- names(list_D_trt) <- names(list_D_ctrl) <- names(list_newH)

    vec_D <- list_D %>% compact %>% pmap_dbl(sum)
    IC <- vec_D

    result <- list(psi = psi, IC = IC
                   # , full_IC = list_D
    )

    # these are cached; unless likelihood is updated, or we force it to update, they shouldn't be changed
    param.list_D <- list_D
    param.list_D_trt <- list_D_trt
    param.list_D_ctrl <- list_D_ctrl
    param.result <- result
  }
}


# suppressMessages(
#   test <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)
# )

{
  suppressMessages(
    updater$update_step(targeted_likelihood, tmle_task, fold_number = "full")
  )
  updater$check_convergence(tmle_task, fold_number = "full") %>% print
  ED_from_estimates(tmle_params %>% lapply(function(x) x$estimates())) %>% print
  tmle_params[[1]]$estimates()$psi
}
