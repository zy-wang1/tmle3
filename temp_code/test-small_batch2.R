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
# source("./temp_code/generate_data.R")

library(parallel)

n_sim <- 8*2
nCores <- detectCores()

timepoint <- 1
if_LY_misspec <- T
if_A_misspec <- F

data_truth <- generate_Zheng_data(B = 1000000, tau = timepoint, seed = 202008, setAM = c(1, 0),
                                  if_LY_misspec = if_LY_misspec)
truth <- data_truth[[timepoint + 1]]$Y %>% mean
truth

results <- mclapply(X = 1:n_sim, mc.cores = nCores, FUN = function(i) {

  # set.seed(1234)

  data_sim <- generate_Zheng_data(B = 20000, tau = timepoint, if_LY_misspec = if_LY_misspec, if_A_misspec = if_A_misspec)
  data_wide <- data.frame(data_sim)

  if (timepoint == 2) {
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
  }  else {
    # data_wide <- data_wide %>% select(-L1_1)
    node_list <- list(L_0 = c("L1_0", "L2_0"),
                      A_1 = "A_1",
                      R_1 = "R_1",
                      Z_1 = "Z_1",
                      # L_1 = "L1_1",
                      Y_1 = "Y_1"
    )
    # node_list$L_1 <- "L_1"
    # names(data_sim[[2]])[grep("L1_1", names(data_sim[[2]]))] <- "L_1"
    # names(data_wide)[grep("L1_1", names(data_wide))] <- "L_1"
  }

  med_spec <- tmle_med(
    treatment_level = 1,
    control_level = 0
  )

  tmle_task <- med_spec$make_tmle_task(data_wide, node_list)

  # choose base learners
  lrnr_glm <- Lrnr_glm$new(outcome_type = "binomial")
  # choose base learners
  lrnr_mean <- make_learner(Lrnr_mean)
  # lrnr_glm <- make_learner(Lrnr_glm)
  lrnr_glm_fast <- Lrnr_glm_fast$new(outcome_type = "binomial")
  # lrnr_ranger50 <- make_learner(Lrnr_ranger, num.trees = 50)
  lrnr_hal_simple <- make_learner(Lrnr_hal9001, max_degree = 5, n_folds = 10)
  lrnr_lasso <- make_learner(Lrnr_glmnet) # al  pha default is 1
  lrnr_ridge <- make_learner(Lrnr_glmnet, alpha = 0)
  lrnr_elasticnet <- make_learner(Lrnr_glmnet, alpha = .5)

  learner_list <- lapply(1:length(tmle_task$npsem), function(s) Lrnr_sl$new(
    learners = list(
      lrnr_glm
      # lrnr_glm_fast
      # ,
      # lrnr_mean,
      # lrnr_elasticnet
      # lrnr_hal_simple
    )
  ))
  names(learner_list) <- names(tmle_task$npsem)  # the first will be ignored; empirical dist. will be used for covariates

  initial_likelihood <- med_spec$make_initial_likelihood(
    tmle_task,
    learner_list
  )

  # # analytic EIC
  # updater <- tmle3_Update$new(convergence_type = "sample_size",
  #                             constrain_step = T,
  #                             optim_delta_epsilon = T,
  #                             one_dimensional=F,
  #                             delta_epsilon=function(x) {
  #                               ifelse(abs(mean(x %>% as.vector)) < sqrt(var(x %>% as.vector)/tmle_task$nrow)/log(tmle_task$nrow),
  #                                      0.00000001,
  #                                      ifelse(mean(x %>% as.vector) > 0, 0.001, -0.001)
  #                               )
  #                             },
  #                             maxit=100
  #                             ,
  #                             cvtmle=F)
  # tlik <- Targeted_Likelihood$new(initial_likelihood,
  #                                 submodel_type_by_node = "EIC" ,
  #                                 updater = updater)

  # logistic
  updater <- tmle3_Update$new(convergence_type = "scaled_var",
                              # fluctuation_type = "standard",
                              # constrain_step = T,
                              # optim_delta_epsilon = T,
                              # one_dimensional=F,
                              # delta_epsilon=function(x) {
                              #   ifelse(abs(mean(x %>% as.vector)) < sqrt(var(x %>% as.vector)/tmle_task$nrow)/log(tmle_task$nrow),
                              #          0.00000001,
                              #          ifelse(mean(x %>% as.vector) > 0, 1, -1)
                              #   )
                              # },
                              maxit=10
                              ,
                              cvtmle=F)
  tlik <- Targeted_Likelihood$new(initial_likelihood,
                                  submodel_type_by_node = "logistic" ,
                                  updater = updater)

  # # projected-EIC, resampling
  # updater <- tmle3_Update$new(convergence_type = "scaled_var",
  #                             constrain_step = T,
  #                             optim_delta_epsilon = T,
  #                             one_dimensional=F,
  #                             delta_epsilon=function(x) {
  #                               ifelse(abs(mean(x %>% as.vector)) < sqrt(var(x %>% as.vector)/tmle_task$nrow)/log(tmle_task$nrow),
  #                                      0.00000001,
  #                                      ifelse(mean(x %>% as.vector) > 0, 0.01, -0.01)
  #                               )
  #                             },
  #                             maxit=100
  #                             ,
  #                             cvtmle=F)
  # tlik <- Targeted_Likelihood$new(initial_likelihood,
  #                                 submodel_type_by_node = "EIC" ,
  #                                 updater = updater)


  tmle_params_no <- med_spec$make_params(tmle_task, initial_likelihood, options = list("tc")
                                         # , static_likelihood = initial_likelihood,
                                         # if_projection = T
                                         # , n_resampling = 0
  )

  # projection param
  tmle_params <- med_spec$make_params(tmle_task, tlik, options = list("tc")
                                      # , static_likelihood = initial_likelihood,
                                      # if_projection = T
                                      # , n_resampling = 50000
  )  # new est

  # tmle_params_raw <- med_spec$make_params(tmle_task, initial_likelihood, options = list("tc")
  #                                     # ,
  #                                     # if_projection = T,
  #                                     # if_resampling = T
  #                                     )

  # tmle_params <- med_spec$make_params(tmle_task, tlik, options = list("tc")
  #                                     # ,
  #                                     # if_projection = T,
  #                                     # if_resampling = T
  #                                     )

  updater$tmle_params <- tmle_params

  suppressMessages(
    nontargeting <- tmle_params_no[[1]]$estimates(tmle_task)
  )
  temp_IC <- nontargeting$IC
  var_D <- var(temp_IC)
  n <- length(temp_IC)
  se <- sqrt(var_D / n)
  CI2 <- nontargeting$psi + 1.96 * se
  CI1 <- nontargeting$psi - 1.96 * se

  est_dr <- nontargeting$psi + mean(temp_IC)
  CI1_dr <- est_dr - 1.96*se
  CI2_dr <- est_dr + 1.96*se

  # suppressMessages(
  #   new <- tmle_params[[1]]$estimates(tmle_task)
  # )
  # temp_IC <- new$IC
  # var_D <- var(temp_IC)
  # n <- length(temp_IC)
  # se <- sqrt(var_D / n)
  # CI2_new <- new$psi + 1.96 * se
  # CI1_new <- new$psi - 1.96 * se
  # new_psi <- new$psi

  new_fit <- fit_tmle3(tmle_task, tlik, tmle_params, updater)
  new_fit
  tlik$cache$tasks %>% lapply(function(x) x$data %>% dim)
  new_psi <- new_fit$summary$tmle_est
  CI1_new <- new_fit$summary$lower
  CI2_new <- new_fit$summary$upper
  updater$step_number

  # # load full_p list first
  # full_task <- tmle_task
  # full_node_names <- names(full_task$npsem)
  # full_data <- full_task$data %>% as.data.frame %>% dplyr::select(-c(id, t))  # exactly the obs data
  # full_variable_names <- colnames(full_data)
  # list_all_predicted_lkd <- lapply(1:length(full_node_names), function(loc_node) {
  #   if (loc_node > 1) {
  #     # currently only support univariate node for t>0
  #     current_variable <- full_task$npsem[[loc_node]]$variables
  #     temp_input <- expand_values(variables = full_variable_names[1:which(full_variable_names == current_variable)])  # all possible inputs
  #     temp_task <- tmle3_Task$new(temp_input, full_task$npsem[1:loc_node])
  #     temp_target_node <- intersect(tmle_params[[1]]$update_nodes, full_node_names[loc_node])
  #     if (length(temp_target_node) == 1) {
  #       setattr(temp_task, "target_nodes", full_node_names[loc_node])
  #       temp_output <- tlik$get_likelihood(temp_task, node = full_node_names[loc_node], fold_number)  # corresponding outputs
  #     } else {
  #       setattr(temp_task, "target_nodes", "no_update")
  #       temp_output <- tlik$get_likelihood(temp_task, node = full_node_names[loc_node], fold_number)  # corresponding outputs
  #     }
  #     data.frame(temp_input, output = temp_output) %>% return
  #   }
  # })
  # names(list_all_predicted_lkd) <- full_node_names

  # n_resampling <- 1000000
  # temp_node_names <- tmle_task$npsem %>% names
  # loc_A <- grep("A", temp_node_names)
  # loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
  # loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
  #
  # temp_input <- tmle_task$get_tmle_node(temp_node_names[1])
  # temp_input <- temp_input[c(
  #   sample(nrow(temp_input), abs(round(n_resampling)), replace = T)
  # ), ]
  # for (i in 2:length(temp_node_names)) {
  #   if (i %in% loc_A) {
  #     temp_input <- cbind(temp_input, 1) %>% as.data.frame()
  #     names(temp_input)[ncol(temp_input)] <- temp_node_names[i]
  #   } else if (i %in% loc_RLY) {
  #     temp_input <- cbind(temp_input, 1) %>% as.data.frame()
  #     names(temp_input)[ncol(temp_input)] <- temp_node_names[i]
  #     temp_input[,ncol(temp_input)] <- rbinom(nrow(temp_input), 1,
  #                                             left_join(temp_input, list_all_predicted_lkd[[i]])$output)
  #   } else if (i %in% loc_Z) {
  #     temp_input_copy <- cbind(temp_input, 1) %>% as.data.frame()
  #     names(temp_input_copy)[ncol(temp_input_copy)] <- temp_node_names[i]
  #     temp_input_copy[, grep("A", colnames(temp_input_copy))] <- 0
  #     temp_input <- cbind(temp_input,
  #                         rbinom(nrow(temp_input), 1,
  #                                left_join(temp_input_copy, list_all_predicted_lkd[[i]])$output)
  #                         ) %>% as.data.frame()
  #     names(temp_input)[ncol(temp_input)] <- temp_node_names[i]
  #   }
  # }
  # MC_est <- temp_input$Y_1 %>% mean
  # new_fit$summary$tmle_est
  MC_est <- 0


  # ipw
  temp_node_names <- tmle_task$npsem %>% names()
  cf_task_treatment <- tmle_params[[1]]$cf_likelihood_treatment$enumerate_cf_tasks(tmle_task)[[1]]
  cf_task_control <- tmle_params[[1]]$cf_likelihood_control$enumerate_cf_tasks(tmle_task)[[1]]
  fold_number = "full"
  intervention_nodes <- union(names(tmle_params[[1]]$intervention_list_treatment), names(tmle_params[[1]]$intervention_list_control))
  obs_data <- tmle_task$data %>% select(-c(id, t))
  obs_variable_names <- colnames(obs_data)
  intervention_variables <- map_chr(tmle_task$npsem[intervention_nodes], ~.x$variables)
  intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
  intervention_levels_treat <- map_dbl(tmle_params[[1]]$intervention_list_treatment, ~.x$value %>% as.character %>% as.numeric)
  intervention_levels_control <- map_dbl(tmle_params[[1]]$intervention_list_control, ~.x$value %>% as.character %>% as.numeric)
  names(intervention_levels_treat) <- names(tmle_params[[1]]$intervention_list_treatment)
  names(intervention_levels_control) <- names(tmle_params[[1]]$intervention_list_control)

  # load full_p list first
  full_task <- tmle_task
  full_node_names <- names(full_task$npsem)
  full_data <- full_task$data %>% as.data.frame %>% dplyr::select(-c(id, t))  # exactly the obs data
  full_variable_names <- colnames(full_data)
  list_all_predicted_lkd <- lapply(1:length(full_node_names), function(loc_node) {
    if (loc_node > 1) {
      # currently only support univariate node for t>0
      current_variable <- full_task$npsem[[loc_node]]$variables
      temp_input <- expand_values(variables = full_variable_names[1:which(full_variable_names == current_variable)])  # all possible inputs
      temp_task <- tmle3_Task$new(temp_input, full_task$npsem[1:loc_node])
      temp_target_node <- intersect(tmle_params[[1]]$update_nodes, full_node_names[loc_node])
      if (length(temp_target_node) == 1) {
        setattr(temp_task, "target_nodes", full_node_names[loc_node])
        temp_output <- tlik$get_likelihood(temp_task, node = full_node_names[loc_node], fold_number)  # corresponding outputs
      } else {
        setattr(temp_task, "target_nodes", "no_update")
        temp_output <- tlik$get_likelihood(temp_task, node = full_node_names[loc_node], fold_number)  # corresponding outputs
      }
      data.frame(temp_input, output = temp_output) %>% return
    }
  })
  names(list_all_predicted_lkd) <- full_node_names
  list_H <- get_obs_H_full(tmle_task, obs_data, current_likelihood = initial_likelihood,
                      cf_task_treatment, cf_task_control,
                      intervention_variables, intervention_levels_treat, intervention_levels_control,
                      fold_number = "full")
  est_ipw <- mean(last(list_H) * data_wide$Y_1)
  est_ipw2 <- weighted.mean(data_wide$Y_1[data_wide$A_1 == 1], last(list_H)[data_wide$A_1 == 1])


  # updater$update_step(likelihood = tlik, tmle_task, fold_number = "full")
  # tmle_params[[1]]$clever_covariates()$IC %>% rowSums %>% mean
  # sd(tmle_params[[1]]$clever_covariates()$IC %>% rowSums) / sqrt(tmle_task$nrow) / log(tmle_task$nrow)
  # var(tmle_params[[1]]$clever_covariates()$IC %>% rowSums)
  # var_D
  # tmle_params[[1]]$estimates()$psi
  # nontargeting$psi
  # tmle_params[[1]]$estimates()$IC %>% hist
  # tmle_params_no[[1]]$estimates()$IC %>% hist
  updater$epsilons

  tmle_params_no[[1]]$list_D %>% lapply(mean)
  tmle_params[[1]]$list_D %>% lapply(mean)



  # capture.output(
  #   updater$update(tlik, tmle_task)
  # )
  # new_est <- updater$tmle_params[[1]]$estimates()
  # new_psi <- new_est$psi
  # nontargeting$psi
  # new_psi
  # temp_IC <- new_est$IC
  # var_D <- var(temp_IC)
  # n <- length(temp_IC)
  # se <- sqrt(var_D / n)
  # CI2_new <- new_est$psi + 1.96 * se
  # CI1_new <- new_est$psi - 1.96 * se

  # est_seq <- est_seq_reg(data_sim)

  return(
    list(non = nontargeting$psi,
         CI1 = CI1,
         CI2 = CI2,
         new = new_psi,
         CI1_new = CI1_new,
         CI2_new = CI2_new,
         step = updater$step_number,
         MC_est = MC_est,
         est_dr = est_dr,
         CI1_dr = CI1_dr,
         CI2_dr = CI2_dr,
         est_ipw = est_ipw,
         est_ipw2 = est_ipw2
         # ,
         # ed = test$ED,
         # threshold = threshold1
         # ,
         # one_cvtmle = test2$estimates[[1]]$psi,
         # step_cvtmle = test2$steps,
         # ed_cvtmle = test2$ED,
         # threshold_cvtmle = sd(tmle_params[[1]]$estimates()$IC) / sqrt(1000) / log(1000)
         # ,
         # last_dir = updater$record_direction %>% last %>% unlist,
         # dir10 = ifelse_vec(length(updater$record_direction) >= 10, updater$record_direction[[10]] %>% unlist, 0)
    )
  )
})



(do.call(rbind, results) %>% as.matrix)
report <- sapply(1:2, function(i) {
  c(mean(( unlist((do.call(rbind, results) %>% as.matrix)[, i * 3 - 2]) - truth)^2),
    mean(( unlist((do.call(rbind, results) %>% as.matrix)[, i * 3 - 2]) - truth)),
    sd( unlist((do.call(rbind, results) %>% as.matrix)[, i * 3 - 2])),
    mean(( unlist((do.call(rbind, results) %>% as.matrix)[, i * 3 - 1]) <= truth & unlist((do.call(rbind, results) %>% as.matrix)[, i * 3]) >= truth)),
    mean(( unlist((do.call(rbind, results) %>% as.matrix)[, i * 3]) - unlist((do.call(rbind, results) %>% as.matrix)[, i * 3 - 1]) ))
  )
})
rownames(report) <- c("MSE", "bias", "sd", "coverage", "width")
report
(do.call(rbind, results) %>% as.matrix)[, 7] %>% unlist

# mean(((do.call(rbind, results) %>% as.matrix)[, 8] %>% unlist) - truth)
mean(((do.call(rbind, results) %>% as.matrix)[, 9] %>% unlist) - truth)

((((do.call(rbind, results) %>% as.matrix)[, 10] %>% unlist) < truth) & (((do.call(rbind, results) %>% as.matrix)[, 11] %>% unlist) > truth)) %>% mean
mean(((do.call(rbind, results) %>% as.matrix)[, 12] %>% unlist) - truth)  # ipw raw
sd(((do.call(rbind, results) %>% as.matrix)[, 12] %>% unlist))
mean(((do.call(rbind, results) %>% as.matrix)[, 13] %>% unlist) - truth)  # ipw 2
sd(((do.call(rbind, results) %>% as.matrix)[, 13] %>% unlist))
